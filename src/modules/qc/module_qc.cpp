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
#include "fastq_parser.h"
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

    

    std::size_t read_index = 0;
    struct time_keep::timer total_time;
    sequential_map<sequence, sample> index_map;
    samplelist_parser samplelist_p;
    std::vector<sample> samplelist = samplelist_p.parse( q_opts ); //samplelist_p;

    
    std::vector<std::pair<std::string,size_t>> index_match_totals;
    for(const auto& curr_index : flexible_idx_data)
        {
            index_match_totals.emplace_back( std::make_pair( curr_index.idx_name, 0));
        }

    fastq_parser fastq_p;

    for( const auto sample : samplelist )
        {
            std::vector<std::tuple<std::string, std::string, std::size_t>> barcodes;
            for( std::size_t curr_barcode = 0; curr_barcode < sample.string_ids.size(); ++curr_barcode )
                {
                    barcodes.emplace_back( std::make_tuple( index_match_totals[curr_barcode].first, sample.string_ids[curr_barcode], 0 ) );
                }
            matched_map.emplace( sample.name, barcodes);
        }

    std::ofstream outfile (q_opts->output_fname, std::ios::out);

    std::string header = "Sample name\t";
    std::size_t index_right_offset = 1;
    std::size_t num_indexes = matched_map.begin()->second.size();

    for( std::size_t curr_index = 0; curr_index < num_indexes; curr_index++)
    {
        std::string header_col = "";
        for( std::size_t index_left_offset = 0; index_left_offset < index_right_offset; index_left_offset++)
        {
            if(index_left_offset != 0)
            {
                header_col.append( " + " );
            }
            std::string index_name = std::get<0>( matched_map.begin()->second[index_left_offset] );
            header_col.append( index_name );
        }

        header.append(header_col);
        header.append( "\t" );
        ++index_right_offset;
    }

    header.append( "\n" );
    outfile << header;
    outfile.close();
    

}
