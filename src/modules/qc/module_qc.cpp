#include "module_qc.h"
#include "options_qc.h"
#include "options.h"
#include "fif_parser.h"
#include <sstream>
#include <stdexcept>
#include <tuple>
#include "et_search.h"
#include "time_keep.h"
#include "maps.h"
#include "sample.h"
#include <unordered_map>
#include "fastq_sequence.h"
#include "fastq_score.h"
#include "samplelist_parser.h"
#include <vector>

std::string module_qc::get_name( )
{
    return "qc";
}

void module_qc::run(options *opts)
{
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::size_t>>> matched_map;
    options_qc *q_opts = (options_qc*)opts;
    std::vector<flex_idx> flexible_idx_data;

    if(q_opts->idx_fname.empty() )
        {
            flexible_idx_data.emplace_back( q_opts->sample_indexes[0],
                                        "r1",
                                        std::get<0>( q_opts->index1_data ),
                                        std::get<1>( q_opts->index1_data ),
                                        std::get<2>( q_opts->index1_data ) );
            flexible_idx_data.emplace_back( q_opts->sample_indexes[1],
                                        "r2",
                                        std::get<0>( q_opts->index2_data ),
                                        std::get<1>( q_opts->index2_data ),
                                        std::get<2>( q_opts->index2_data ) );
        }

    else
        {
            fif_parser flex_parser;
            flexible_idx_data = flex_parser.parse( q_opts->idx_fname );
        }

    

    //std::size_t read_index = 0;
    //struct time_keep::timer total_time;
    sequential_map<sequence, sample> index_map;
    samplelist_parser samplelist_p;

    if( !q_opts->idx_fname.empty() )
        {
            for( std::size_t curr_row = 0; curr_row < flexible_idx_data.size(); ++curr_row )
                {
                    q_opts->sample_indexes.emplace_back( flexible_idx_data[curr_row].idx_name );
                }
        }

    std::vector<sample> samplelist = samplelist_p.parse( q_opts ); //samplelist_p;

    std::vector<std::pair<std::string,size_t>> index_match_totals;
    for(const auto& curr_index : flexible_idx_data)
        {
            index_match_totals.emplace_back( std::make_pair( curr_index.idx_name, 0));
        }
    
    fastq_parser fastq_p;
    std::vector<fastq_sequence> index_seqs = fastq_p.parse( q_opts->input_r1_fname );

    for( const auto seq : index_seqs )
        {
            std::cout << seq.seq << std::endl;
            
        }

    //set up all index maps
    std::map<std::string, std::size_t> index_loc;
    //std::vector<std::map<std::string, std::size_t>> index_maps;
    std::map<std::string, std::size_t> index_maps[q_opts->sample_indexes.size()];

    for(std::size_t i = 0; i < q_opts->sample_indexes.size(); ++i)
    {
        index_loc[q_opts->sample_indexes.at(i)] = i;
        std::map<std::string, std::size_t> map;
        index_maps[i] = map;
    }

    for( const auto sample : samplelist )
        {
            std::vector<std::tuple<std::string, std::string, std::size_t>> barcodes;
            std::string barcode = "";
            for( std::size_t curr_barcode = 0; curr_barcode < sample.string_ids.size(); ++curr_barcode )
                {
                    barcodes.emplace_back( std::make_tuple( index_match_totals[curr_barcode].first, sample.string_ids[curr_barcode], 0 ) );

                    index_maps[ index_loc[ index_match_totals[ curr_barcode ].first ]][sample.string_ids[ curr_barcode ]]++;

                }
        }

    std::vector<std::vector<std::string>> ordered_indexes;

    for(std::size_t i = 0; i < index_loc.size(); ++i)
        {
            std::vector<std::string> vec;
            ordered_indexes.emplace_back(vec);
            for( const auto row : index_maps[i] )
                {
                    ordered_indexes.at(i).emplace_back(row.first);
                }
        }

    std::vector<std::string> out_dumb;
    std::vector<std::string> all_possible = get_combined((std::size_t)0, ordered_indexes, out_dumb);

    for( const auto id_comb : all_possible)
        {
            std::cout << id_comb << std::endl;
        }

    //std::ofstream outfile (q_opts->output_fname, std::ios::out);

    //outfile.close();

}

std::vector<std::string> module_qc::get_combined( std::size_t iteration, std::vector<std::vector<std::string>> indata, std::vector<std::string> outdata )
{
    if( iteration == indata.size( ) )
        {
            return outdata;
        }
    else if( iteration == 0 )
        {
            std::vector<std::string> out;

            for( const auto id : indata.at( iteration ) )
                {
                    out.emplace_back( id );
                }

            return get_combined( ++iteration, indata, out );
        }
    else
        {
            std::vector<std::string> out;

            for( const auto init : outdata )
                {
                    for( const auto ad : indata.at( iteration ) )
                        {
                            out.emplace_back( init + "\t" + ad );
                        }
                }

            return get_combined( ++ iteration, indata, out );
        }
}
