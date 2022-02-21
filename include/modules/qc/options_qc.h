#ifndef OPTIONS_QC_HH_INCLUDED
#define OPTIONS_QC_HH_INCLUDED

#include <string>
#include <tuple>
#include <sstream>
#include <vector>
#include "options.h"

class options_qc : public options
{
    public:
        options_qc(): DEFAULT_READ_PER_LOOP( 100000 ), DEFAULT_OUTPUT_FNAME( "output.csv" ){}
        std::string idx_fname;
        std::string input_r1_fname;
        std::string input_r2_fname;
        std::tuple<std::size_t, std::size_t, std::size_t> index1_data;
        std::tuple<std::size_t, std::size_t, std::size_t> index2_data;
        std::string samplelist_fname;
        std::vector<std::string> sample_indexes;
        std::string samplename;
        std::string output_fname;

        long int read_per_loop;
        const long int DEFAULT_READ_PER_LOOP;
        const std::string DEFAULT_OUTPUT_FNAME;

        std::string get_arguments();

        std::string tup_to_string( std::tuple<std::size_t,
                               std::size_t,
                               std::size_t>& data
                             );

        void set_info( std::tuple<std::size_t, std::size_t, std::size_t>
                              options_qc:: *member, std::string info
                            );
};

#endif