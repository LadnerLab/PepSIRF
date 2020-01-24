#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "omp_opt.h"
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range.hpp>
#include <algorithm>

#include "fs_tools.h"
#include "module_deconv.h"
#include "kmer_tools.h"
#include "time_keep.h"
#include "fasta_parser.h"
#include "setops.h"
#include "distance_tools.h"
#include "overlap_data.h"
#include "distance_matrix.h"

module_deconv::module_deconv() = default;

std::string module_deconv::get_name()
{
    return "Deconv";
}

void module_deconv::run( options *opts )
{
    options_deconv *d_opts = ( options_deconv * ) opts;
    time_keep::timer timer;

    timer.start();

    fs_tools::path input_base{ d_opts->enriched_fname };

    if( fs_tools::is_directory( input_base ) )
        {
            if( d_opts->output_fname == "deconv_output.tsv" )
                {
                    d_opts->output_fname = "deconv_output";
                }

            fs_tools::path output_base{ d_opts->output_fname };
            bool output_existed = !fs_tools::create_directory( output_base );

            if( output_existed )
                {
                    std::cout << "WARNING: Output directory '"
                              << output_base.string()
                              << "' already exists. Any "
                              << "files with colliding names will be "
                              << "overwritten.\n";
                }

            auto in_dir_iter = fs_tools::directory_iterator( input_base );

            for( auto& input_f : boost::make_iterator_range( in_dir_iter, {} ) )
                {
                    std::string file_name = input_f.path().filename().string();
                    fs_tools::path in_path = input_f;
                    fs_tools::path out_path = output_base/( file_name + d_opts->outfile_suffix );

                    d_opts->enriched_fname = in_path.string();
                    d_opts->output_fname   = out_path.string();

                    choose_kmers( d_opts );
                }
        }
    else
        {
            choose_kmers( d_opts );
        }

    timer.stop();

    std::cout << "Took " << time_keep::get_elapsed( timer ) << " second(s).\n";
}

void module_deconv::choose_kmers( options_deconv *opts )
{
    options_deconv *d_opts = opts;

    auto enriched_species = parse_enriched_file( d_opts->enriched_fname );
    auto pep_species_vec  = parse_linked_file( d_opts->linked_fname );

    evaluation_strategy::score_strategy score_strat   = get_evaluation_strategy( d_opts );
    evaluation_strategy::filter_strategy filter_strat = get_filter_method( d_opts );
    evaluation_strategy::tie_eval_strategy tie_eval_strat = get_tie_eval_strategy( d_opts );

    // filter out the peptides that are not enriched
    auto it = std::remove_if( pep_species_vec.begin(), pep_species_vec.end(),
                              [&]( std::pair<std::string,std::vector<std::pair<std::string,double>>>& i) -> bool 
                              { return enriched_species.find( std::get<0>( i ) )
                                      == enriched_species.end();
                              }
                            );
    pep_species_vec.erase( it, pep_species_vec.end() );

    std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>
        peptide_assignment_global;

    std::unordered_map<species_id<std::string>,
                       scored_entity<peptide, double>
                       >
        species_highest_peptide;

    // add what species were originally shown to "hit" which peptides
    for( const auto& x : pep_species_vec )
        {
            auto ref = peptide_assignment_global
                .emplace( x.first, std::vector<std::pair<std::string,double>>() ).first;
            for( const auto& i : x.second )
                {
                    ref->second.emplace_back( i );
                }
        }

    std::size_t thresh = d_opts->threshold;

    omp_set_num_threads( d_opts->single_threaded ? 1 : 2 );

    std::unordered_map<std::string, std::vector<std::string>> id_pep_map;
    std::unordered_map<std::string, std::vector<std::pair<std::string,double>>>
        pep_id_map;

    // vector holding the species ids with the highest count
    std::vector<std::pair<std::string, double>> species_scores;
    std::unordered_map<std::string, double> species_peptide_counts;

    std::vector<std::pair<species_data, bool>> output_counts;

    std::string id_name_map_fname;

    std::map<std::string,std::string> name_id_map;
    std::map<std::string,std::string>* name_id_map_ptr = nullptr;

    if( !util::empty( d_opts->id_name_map_fname ) )
        {
            parse_name_map( d_opts->id_name_map_fname, name_id_map );
            name_id_map_ptr = &name_id_map;
        }

    auto filter = [&]( evaluation_strategy::filter_strategy filter_strat )
        {
            if( filter_strat
                == evaluation_strategy::filter_strategy::SCORE_FILTER )
                {
                    filter_counts<std::string,double>
                    ( species_scores, thresh );


                }
            else if( filter_strat
                      == evaluation_strategy::filter_strategy::COUNT_FILTER
                   )
                {
                    
                    // recreate species_scores
                    filter_counts<std::string,double>
                    ( species_peptide_counts, thresh );
                }
        };

    auto make_map = [&]()
        {
        #pragma omp parallel
            {
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        id_to_pep( id_pep_map, pep_species_vec );
                    }

                    #pragma omp section
                    {
                        pep_to_id( pep_id_map, pep_species_vec );
                    }
                }

                // implicit barrier
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        score_species( species_scores, id_pep_map, pep_id_map,
                                       score_strat
                                     );
                    }

                    #pragma omp section
                    {
                        get_species_counts_per_peptide( id_pep_map,
                                                        species_peptide_counts
                                                      );
                    }
                }
            }
        };

    auto make_map_and_filter = [&]( evaluation_strategy::filter_strategy filter_strat )
        {
            make_map();
            filter( filter_strat );
        };

    make_map();


    auto write_round_scores = [&]( std::size_t round_no )
        {
            std::unordered_map<std::string,std::pair<double,double>> round_scores;

            std::ofstream out_f;
            std::string fname;
            fs_tools::create_fname( fname,
                                    d_opts->orig_scores_dname,
                                    "round_", round_no
                                  );

            // populate the unfiltered scores with counts, scores for all species
            combine_count_and_score( round_scores,
                                     species_peptide_counts,
                                     species_scores
                                   );

            out_f.open( fname );

            write_scores( out_f,
                          name_id_map_ptr,
                          round_scores
                        );

            out_f.close();
        };

    filter( filter_strat );

    // used to track the original scores for each species
    std::unordered_map<std::string,std::pair<double,double>> original_scores;
    std::size_t round_no = 0;

    if( !util::empty( d_opts->orig_scores_dname ) )
        {
            write_round_scores( round_no );
        }

    combine_count_and_score( original_scores,
                             species_peptide_counts,
                             species_scores
                           );

    std::unordered_map<std::string,std::vector<std::string>> peptide_assignment_map;

    std::unordered_map<std::string,std::unordered_map<std::string,double>>
        pep_spec_map_w_counts;

    if( tie_eval_strat
        == evaluation_strategy
           ::tie_eval_strategy
             ::PERCENT_TIE_EVAL
        || tie_eval_strat
          == evaluation_strategy
        ::tie_eval_strategy
        ::SUMMATION_SCORING_TIE_EVAL
        )
        {
            for( const auto& pep : pep_id_map )
                {

                    std::unordered_map<std::string,double>
                        val_map;

                    util::pairs_to_map
                        ( val_map,
                          pep.second.begin(),
                          pep.second.end()
                        );

                    pep_spec_map_w_counts
                        .insert( std::make_pair( pep.first,
                                                 val_map
                                               )
                               );
                }
        }

    std::unordered_map<std::string, std::vector<scored_peptide<double>>>
        species_peptide_scores;

    std::unordered_map<std::string, scored_peptide<double>>
        species_with_highest_peptide;

    while( species_peptide_counts.size()
           && species_scores[ 0 ].second > thresh )
        {
            pep_species_vec.clear();
            std::vector<std::pair<std::string,double>>
                tie_candidates;

            std::vector<std::pair<std::string, double>> tied_species;

            tie_data::tie_type tie = tie_data::tie_type::SINGLE_WAY_TIE;

            if( !use_ratio_score_tie_thresh( d_opts->score_tie_threshold ) )
                {
                        tie = get_tie_candidates( tie_candidates,
                                                  species_scores,
                                                  thresh,
                                                  d_opts->score_tie_threshold,
                                                  util::difference<double>()
                                                );
                }
            else
                {
                        tie = get_tie_candidates( tie_candidates,
                                                  species_scores,
                                                  thresh,
                                                  d_opts->score_tie_threshold,
                                                  util::ratio<double>()
                                                );
                }

            handle_ties( tied_species,
                         id_pep_map,
                         pep_spec_map_w_counts,
                         tie_candidates,
                         tie_eval_strat,
                         tie,
                         d_opts->score_overlap_threshold
                       );

            int is_tie = tied_species.size() - 1;

            // calculate the score each peptide attributes to
            // its species score
            score_species_peptides( species_peptide_scores,
                                    id_pep_map,
                                    pep_id_map,
                                    score_strat
                                    );

            // get the highest of these scores for each species
            get_highest_score_per_species( species_with_highest_peptide,
                                           species_peptide_scores
                                           );


            for( auto& tied_peptide : tied_species )
                {
                    std::string id = tied_peptide.first;
                    double score   = tied_peptide.second;

                    auto new_data = species_data( species_id<std::string>( id ),
                                                  score,
                                                  species_peptide_counts
                                                  .find( id )->second,
                                                  species_with_highest_peptide[ id ]
                                                );

                    output_counts.emplace_back( std::make_pair(
                                                                new_data,
                                                                is_tie > 0
                                                                )
                                                );

                    // reduce the number of ties
                    --is_tie;

                    auto current_peptides = id_pep_map.find( id )->second;

                    for( auto& pep_id : current_peptides )
                        {
                            peptide_assignment_map
                                .emplace( pep_id, std::vector<std::string>() );

                            peptide_assignment_map
                                .find( pep_id )->second.push_back( id );

                            pep_id_map.erase( pep_id );
                        }
                }

            for( auto it = pep_id_map.begin(); it != pep_id_map.end(); ++it )
                {
                    pep_species_vec.emplace_back( it->first, it->second );
                }

            if( !util::empty( d_opts->orig_scores_dname ) )
                {
                    write_round_scores( round_no );
                }

            tied_species.clear();
            tie_candidates.clear();
            id_pep_map.clear();
            pep_id_map.clear();
            species_scores.clear();
            species_peptide_counts.clear();
            species_peptide_scores.clear();
            species_with_highest_peptide.clear();

            make_map_and_filter( filter_strat );

            ++round_no;
        }

    write_outputs( d_opts->output_fname,
                   name_id_map_ptr,
                   output_counts,
                   original_scores
                 );

    if( d_opts->species_peptides_out.compare( "" ) )
        {
            write_species_assign_map( d_opts->species_peptides_out,
                                      peptide_assignment_global,
                                      peptide_assignment_map
                                    );
        }
}


std::vector<std::pair<std::string,std::vector<std::pair<std::string,double>>>>
module_deconv::parse_linked_file( std::string fname )
{
    std::ifstream in_stream( fname );

    std::vector<std::pair<std::string,std::vector<std::pair<std::string,double>>>>
        ret_vec;

    std::string line;

    int lineno = 0;
    const boost::regex id_count_re{ "([^\\s:]+):{0,1}(\\d*\\.{0,1}\\d*)" };

    while( std::getline( in_stream, line ) )
        {
            if( lineno )
                {
                    std::vector<std::string> split_line;
                    boost::trim_right( line );
                    boost::split( split_line, line,
                                  boost::is_any_of( "\t" )
                                );


                    std::vector<std::string> comma_delimited;
                    if( split_line.size() > 1 )
                        {
                            boost::split( comma_delimited, split_line[ 1 ],
                                          boost::is_any_of( "," )
                                        );

                            std::vector<std::pair<std::string,double>> id_ints;

                            for( auto& item : comma_delimited )
                                {
                                    boost::smatch match;
                                    boost::regex_search( item, match, id_count_re );

                                    // matched 'id', no count found
                                    if( match[ 1 ] != ""
                                        && match[ 2 ] == "" )
                                        {
                                            id_ints.emplace_back(
                                                                 std::make_pair(
                                                                 ( match[ 1 ] ),
                                                                 boost::lexical_cast<double>
                                                                 ( 1 )
                                                                                )
                                                                );
                                        }
                                    // matched 'id:count'
                                    else if( match[ 1 ] != ""
                                             && match[ 2 ] != "" ) 
                                        {
                                            id_ints.emplace_back(
                                                                 std::make_pair(
                                                                 ( match[ 1 ] ),
                                                                 boost::lexical_cast<double>
                                                                 ( match[ 2 ] )
                                                                                )
                                                                );
                                        }
                                    else
                                        {
                                            std::cout << "Failed to line: \n";
                                            std::cout << item << "\n";
                                        }

                                }
                            ret_vec.emplace_back( split_line[ 0 ], id_ints );

                        }
                }
            ++lineno;
        }

    return ret_vec;
}

void module_deconv::id_to_pep( std::unordered_map<std::string, std::vector<std::string>>&
                               id_pep_map,
                               std::vector<std::pair<std::string, std::vector<std::pair<std::string,double>>>>&
                               pep_species_vec )
{
    for( auto it = pep_species_vec.begin(); it != pep_species_vec.end(); ++it )
        {
            std::string& pep = std::get<0>( *it );

            for( auto inner = std::get<1>( *it ).begin();
                     inner != std::get<1>( *it ).end();
                 ++inner
               )
                {
                    auto find = id_pep_map.find( std::get<0>( *inner ) );
                    if( find == id_pep_map.end() )
                        {
                            find = std::get<0>(
                                               id_pep_map.emplace( std::get<0>( *inner ), std::vector<std::string>() )
                                              );
                        }
                    find->second.push_back( pep );
                }
        }
}
double module_deconv::score_peptide_for_species( const peptide& peptide,
                                std::unordered_map
                                <std::string,std::vector<std::pair<std::string,double>>>&
                                spec_count_map,
                                std::string id,
                                evaluation_strategy::score_strategy score_strat
                                )
{
    double score = 0;

    const std::string& peptide_str = peptide.get_sequence();

    if( score_strat == evaluation_strategy::score_strategy::INTEGER_SCORING )
        {
            score =  1.0;
        }
    else if( score_strat == evaluation_strategy::score_strategy::FRACTIONAL_SCORING )
        {
            score = 1.0 / (double) spec_count_map[ peptide_str ].size();
        }
    else if( score_strat == evaluation_strategy::score_strategy::SUMMATION_SCORING )
        {
                    for( auto& iter : spec_count_map[ peptide_str ] )
                        {
                            if( std::get<0>( iter ) == id )
                                {
                                    score += std::get<1>( iter );
                                }
                        }
        }
    return score;
}


double module_deconv::get_score( std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                                 spec_count_map,
                                 std::string id,
                                 std::vector<std::string>& peptides,
                                 evaluation_strategy::score_strategy strat
                               )
{
    double score = 0;

    if( strat == evaluation_strategy::score_strategy::INTEGER_SCORING )
        {
            return (double) peptides.size();
        }

    else if( strat == evaluation_strategy::score_strategy::FRACTIONAL_SCORING )
        {
            std::size_t index = 0;

            #pragma omp parallel for private( index )           \
             shared( spec_count_map ) reduction( +:score )
            for( index = 0; index < peptides.size(); ++index )
                {
                    score += 1.0 / (double) spec_count_map[ peptides[ index ] ].size();
                }
        }

    else if( strat == evaluation_strategy::score_strategy::SUMMATION_SCORING )
        {
            std::size_t index = 0;

            #pragma omp parallel for private( index )           \
             shared( spec_count_map ) reduction( +:score )
            for( index = 0; index < peptides.size(); ++index )
                {
                    for( auto& iter : spec_count_map[ peptides[ index ] ] )
                        {
                            if( std::get<0>( iter ) == id )
                                {
                                    score += std::get<1>( iter );
                                }
                        }
                }
        }
    
    return score;
}

void module_deconv::pep_to_id( std::unordered_map<std::string, std::vector<std::pair<std::string,double>>>&
                               pep_id_map,
                               std::vector<std::pair<std::string, std::vector<std::pair<std::string,double>>>>&
                               pep_species_vec
                             )
{
    for( auto it = pep_species_vec.begin(); it != pep_species_vec.end(); ++it )
        {
            pep_id_map.emplace( std::get<0>( *it ),
                                std::get<1>( *it )
                              );
        }
}

void module_deconv::score_species( std::vector<std::pair<std::string, double>>&
                                   id_counts,
                                   std::unordered_map<std::string,std::vector<std::string>>&
                                   id_count_map,
                                   std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                                   spec_count_map,
                                   evaluation_strategy::score_strategy strat
                                )

{
    id_counts.reserve( id_count_map.size() );

    for( auto it = id_count_map.begin(); it != id_count_map.end(); ++it )
        {
            double score = get_score( spec_count_map,
                                      it->first,
                                      it->second,
                                      strat
                                    );
            id_counts.emplace_back( it->first, score );
        }
    std::sort( id_counts.begin(), id_counts.end(),
               util::compare_pair_non_increasing<std::string, double>()
             );

}

void module_deconv::score_species_peptides(
                   std::unordered_map<std::string,
                   std::vector<scored_peptide<double>>
                   >& dest,
                   std::unordered_map<std::string,std::vector<std::string>>&
                   id_count_map,
                   std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                   spec_count_map,
                   evaluation_strategy::score_strategy strat
                                )
    {

        dest.reserve( id_count_map.size() );

        for( const auto& species : id_count_map )
            {
                const auto& species_id = species.first;
                const auto& peptides   = species.second;

                dest[ species_id ] = std::vector<scored_peptide<double>>();
                dest[ species_id ].reserve( peptides.size() );

                for( const auto& pep : peptides )
                    {
                        double score = score_peptide_for_species
                            ( peptide( pep ), spec_count_map, species_id, strat );

                        auto new_pep = scored_peptide<double>( pep, score );
                        dest[ species_id ].push_back( new_pep );
                    }
            }
    }

void module_deconv::get_highest_score_per_species( std::unordered_map<std::string,
                                                   scored_peptide<double>>& dest,
                                                   const std::unordered_map<
                                                   std::string,
                                                   std::vector<scored_peptide<double>>
                                                   >&
                                                   species_peptide_scores
                                                 )
{
    for( const auto& species_w_peptides : species_peptide_scores )
        {
            const auto& species_id = species_w_peptides.first;
            const auto& peptides = species_w_peptides.second;

            const auto max_peptide = std::max_element( peptides.begin(),
                                                       peptides.end()
                                                     );

            dest.emplace( species_id, *max_peptide );
        }

}



void module_deconv::write_outputs( std::string out_name,
                                   std::map<std::string,std::string>*
                                   id_name_map,
                                   std::vector<
                                   std::pair<species_data, bool>
                                   >&
                                   out_counts,
                                   std::unordered_map<std::string,std::pair<double,double>>&
                                   original_scores
                                 )

{
    std::ofstream out_file( out_name );

    if( id_name_map != nullptr )
        {
            out_file << "Species Name\t";
        }


    out_file << "Species ID\tCount\tScore\tOriginal Count\tOriginal Score\tMax Probe Score\n";

    bool tied = false;
    std::vector<std::pair<species_data,bool>> tied_items;

    for( auto it = out_counts.begin();
         it != out_counts.end();
         ++it
       )
        {
            auto tied_item = it;
            while( tied_item->second ) // *it and the next species are tied, report together
                {
                    tied_item = std::next( tied_item, 1 );
                    tied_items.push_back( *tied_item );
                    tied = true;
                }

            if( id_name_map != nullptr )
                {
                    for( auto tied_i : tied_items )
                        {
                            to_stream_if( out_file, tied, 
                                          get_map_value( id_name_map,
                                                         tied_i.first.get_id(),
                                                         tied_i.first.get_id()
                                                       ),
                                          ","
                                );
                        }
                    out_file << get_map_value( id_name_map,
                                               it->first.get_id(),
                                               it->first.get_id()
                                             ) << "\t";
                }

            else
                {
                    out_file << "\t";
                }

            auto orig_id = it->first.get_id();

            // species id for both (both are only written if tied is true)
            for( auto tied_i : tied_items )
                {

                    auto tied_id = tied_i.first;
                    to_stream_if( out_file, tied, tied_id.get_id(), "," );
                }
            out_file << orig_id << "\t";

            // count for both
            for( auto tied_i : tied_items )
                {
                    to_stream_if( out_file, tied, tied_i.first.get_count(), "," );
                }

            out_file << it->first.get_count() << "\t";

            // score for both 
            for( auto tied_i : tied_items )
                {
                    to_stream_if( out_file, tied, tied_i.first.get_score(), "," );
                }

            out_file << it->first.get_score() << "\t";

            // original count for both
            for( auto tied_i : tied_items )
                {
                    auto tied_id = std::get<0>( tied_i );
                    to_stream_if( out_file, tied,
                                  original_scores
                                  .find( tied_id.get_id() )->second.first,
                                  ","
                                  );
                }

            out_file << original_scores.find( orig_id )->second.first << "\t";

            // original score for both
            for( auto tied_i : tied_items )
                {
                    auto tied_id = std::get<0>( tied_i );
                    to_stream_if( out_file, tied,
                                  original_scores
                                  .find( tied_id.get_id() )->second.second,
                                  ","
                                );
                }

            out_file << original_scores.find( orig_id )->second.second << "\t";

            for( auto tied_i : tied_items )
                {
                    out_file << tied_i.first.get_highest_scoring_peptide().get_score() << ",";
                }
            out_file << it->first.get_highest_scoring_peptide().get_score() << "\n";

            if( tied )
                {
                    it = tied_item;
                    tied = false;
                    tied_items.clear();
                }
        }
}

sequential_set<std::string>
module_deconv::parse_enriched_file( std::string f_name )
{
    sequential_set<std::string>
        return_val;

    std::ifstream in_file( f_name );
    std::string line;

    while( std::getline( in_file, line ) )
        {
            boost::trim_right( line );
            return_val.insert( line );
        }
    return return_val;
}

evaluation_strategy::filter_strategy
module_deconv::get_filter_method( options_deconv *opts )
{
    if( opts->score_filtering )
        {
            return evaluation_strategy::filter_strategy::SCORE_FILTER;
        }
    return evaluation_strategy::filter_strategy::COUNT_FILTER;
}

evaluation_strategy::score_strategy
module_deconv::get_evaluation_strategy( options_deconv *opts )
{
    if( !( opts->fractional_scoring
           || opts->summation_scoring
          )
      )
        {
            return evaluation_strategy::score_strategy::INTEGER_SCORING;
        }

    return  opts->fractional_scoring ?
             evaluation_strategy::score_strategy::FRACTIONAL_SCORING :
             evaluation_strategy::score_strategy::SUMMATION_SCORING;
}

evaluation_strategy::tie_eval_strategy
module_deconv::get_tie_eval_strategy( options_deconv *opts )
{
    if( opts->summation_scoring )
        {
            return evaluation_strategy::tie_eval_strategy::SUMMATION_SCORING_TIE_EVAL;
        }

    return use_ratio_overlap_threshold( opts->score_overlap_threshold ) ?
        evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL        :
        evaluation_strategy::tie_eval_strategy::INTEGER_TIE_EVAL;
}


void
module_deconv::parse_name_map( std::string fname,
                               std::map<std::string,std::string>& name_map
                             )
{
    std::ifstream in_stream( fname );
    std::size_t lineno = 0;

    std::string line;

    while( std::getline( in_stream, line ) )
        {
            std::vector<std::string> split_line;
            if( lineno )
                {

                    boost::trim_right( line );
                    boost::split( split_line, line,
                                  boost::is_any_of( "|" )
                                );

                    // if this is a species-level ID,
                    // then this will provide the mapping.
                    // Otherwise some space is wasted, but that's not
                    // a major concern.
                    boost::trim( split_line[ 0 ] );
                    boost::trim( split_line[ 1 ] );

                    std::string id   = ( split_line[ 0 ] );
                    std::string name = split_line[ 1 ];
                    
                    name_map.insert( std::make_pair( id, name ) );
                }

            ++lineno;
        }
}


void module_deconv::
get_species_counts_per_peptide( std::unordered_map<std::string, std::vector<std::string>>&
                                id_pep_map,
                                std::unordered_map<std::string,double>& pep_counts
                              )
{
    for( auto count = id_pep_map.begin();
         count != id_pep_map.end();
         ++count
       )
        {
            pep_counts.emplace( count->first,
                                count->second.size()
                              );
        }
}

void
module_deconv::handle_ties( std::vector<std::pair<std::string,double>>&
                            dest_vec,
                            std::unordered_map<std::string, std::vector<std::string>>&
                            id_pep_map,
                            std::unordered_map<std::string,std::unordered_map<std::string,double>>&
                            pep_species_map_wcounts,
                            std::vector<std::pair<std::string,double>>&
                            tie_candidates,
                            evaluation_strategy::tie_eval_strategy
                            tie_evaluation_strategy,
                            tie_data::tie_type tie_type,
                            double overlap_threshold
                          )

{

    // '1-way tie': there is a single candidate
    if( tie_type == tie_data::tie_type::SINGLE_WAY_TIE )
        {
            dest_vec.push_back( tie_candidates[ 0 ] );
        }

    // two more cases: 2-way tie and k-way tie
    else if( tie_type == tie_data::tie_type::TWO_WAY_TIE )
        {
            // if the overlap between the two is high,
            // report them together. Otherwise, report the first
            // item. High is defined by the overlap_threshold
            if( calculate_overlap(  
                                   id_pep_map,
                                   pep_species_map_wcounts,
                                   tie_candidates[ 0 ].first,
                                   tie_candidates[ 1 ].first,
                                   tie_evaluation_strategy
                                 ).sufficient( overlap_threshold )
              )
                {
                    dest_vec.insert( dest_vec.end(),
                                     tie_candidates.begin(),
                                     tie_candidates.end()
                                   );
                }
            else
                {
                    // if there is not significant overlap,
                    // we only consider the first item.
                    dest_vec.emplace_back( tie_candidates.at( 0 ) );
                }
        }
    else // k-way tie
        {
            handle_kway_tie( dest_vec,
                             id_pep_map,
                             pep_species_map_wcounts,
                             tie_candidates,
                             tie_evaluation_strategy,
                             overlap_threshold
                           );
        }
}

tie_data::tie_type
module_deconv::get_tie_type( std::size_t to_convert )
{
    tie_data::tie_type ret_val;

    switch( to_convert )
        {
        case 0:
            throw std::runtime_error( "to_convert must be > 0" );
            break; 
        case 1:
            ret_val =  tie_data::tie_type::SINGLE_WAY_TIE;
            break;
        case 2:
            ret_val = tie_data::tie_type::TWO_WAY_TIE;
            break;
        default:
            ret_val = tie_data::tie_type::K_WAY_TIE;
            break;
        }
    return ret_val;
}

void
module_deconv
::write_species_assign_map( std::string fname,
                            std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                            peptide_assign_original,
                            std::unordered_map<std::string,std::vector<std::string>>&
                            out_map
                          )

{
                                          
    std::ofstream ofs( fname, std::ofstream::out );
    ofs << "Peptide\tAssigned Ids\tAll IDs\n";
    
    for( auto curr = out_map.begin(); curr != out_map.end(); ++curr )
        {
            ofs << curr->first << "\t";
            ofs << boost::algorithm::join( curr->second, "," );
            ofs << "\t";

            std::vector<std::string> id_count_pairs;
            const auto& species_wcounts = peptide_assign_original.find( curr->first )->second;

            for( const auto& x : species_wcounts )
                {
                    std::string curr_str;

                    curr_str.append( x.first );
                    curr_str.append( ":" );
                    curr_str.append( boost::lexical_cast<std::string>( x.second ) ); 

                    id_count_pairs.push_back( curr_str );
                }


            ofs << boost::algorithm::join(
                                          id_count_pairs,
                                          ","
                                         );
            ofs << "\n";
        }
    ofs << "\n";
}

void
module_deconv
::handle_kway_tie( 
                   std::vector<std::pair<std::string,double>>& tie_outputs,
                   std::unordered_map<std::string, std::vector<std::string>>& id_pep_map,
                   std::unordered_map<std::string,std::unordered_map<std::string,double>>&
                   pep_species_map_wcounts,
                   std::vector<std::pair<std::string,double>>& tie_candidates,
                   evaluation_strategy::tie_eval_strategy eval_strat,
                   const double ovlp_threshold
                 )
{
    distance_matrix<overlap_data<double>> pairwise_distances( tie_candidates.size() );

    auto distance = [&]( const std::pair<std::string, double>& first,
                         const std::pair<std::string, double>& second
                       )
        {
            return calculate_overlap( id_pep_map,
                                      pep_species_map_wcounts,
                                      first.first,
                                      second.first,
                                      eval_strat
                                    );
        };

    // compute pairwise distances for each of the species
    distance_tools
        ::pairwise_distances( pairwise_distances,
                              tie_candidates.begin(),
                              tie_candidates.end(),
                              distance
                            );

    // index, index (inner), a_to_b, b_to_a
    std::tuple<std::size_t, 
               std::size_t,
               overlap_data<double>
               > max_match = std::make_tuple( 0, 0, overlap_data<double>{0,0} );

    std::size_t max_index_outer = 0;
    std::size_t max_index_inner = 0;
    overlap_data<double> max_overlap = std::get<2>( max_match );

    for( std::size_t index = 0; index < pairwise_distances.size(); ++index )
        {
            std::size_t inner_index = 0;
            for( inner_index = index + 1;
                 inner_index < pairwise_distances[ index ].size();
                 ++inner_index
               )
                {
                    auto current = pairwise_distances[ index ][ inner_index ];

                    // does this pair have greater overlap than the highest we've found?
                    bool greater_ovlp = max_overlap < current;

                    // does this pair have a higher score than the highest we've found?
                    bool higher_score = tie_candidates[ max_index_outer ].second
                                        > tie_candidates[ index ].second
                                          && tie_candidates[ max_index_inner ].second
                                             > tie_candidates[ inner_index ].second;

                    if( greater_ovlp && higher_score )
                        {
                            max_match = std::make_tuple(
                                                        index,
                                                        inner_index,
                                                        current
                                                        );

                            max_index_outer = index;
                            max_index_inner = inner_index;
                        }
                }

        }

    tie_outputs.emplace_back( tie_candidates[ max_index_outer ] );

    if( std::get<2>( max_match ).sufficient( ovlp_threshold ) )
        {
            // get the species that are in the tie
            std::size_t max_outer = std::get<0>( max_match );
            std::size_t max_inner = std::get<1>( max_match );

            tie_outputs.emplace_back( tie_candidates[ max_inner ] );

            
            for( std::size_t index = 0;
                 index < pairwise_distances[ max_outer ].size();
                 ++index
               )
                {
                    if( !( index == max_inner || index == max_outer )
                          && pairwise_distances[ max_outer ][ index ]
                             .sufficient( ovlp_threshold )
                            && pairwise_distances[ max_inner ][ index ].sufficient( ovlp_threshold )
                      )
                        {
                            tie_outputs.emplace_back( tie_candidates[ index ] );
                        }
                }
        }
}

bool module_deconv              
::use_ratio_score_tie_thresh( double threshold )
{
    return !util::is_integer( threshold );
}


bool module_deconv
::use_ratio_overlap_threshold( double threshold )
{
    return !util::is_integer( threshold );
}
