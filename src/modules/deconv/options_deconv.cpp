#include <algorithm>
#include <sstream>

#include "options_deconv.h"

std::string options_deconv::get_arguments()
{
    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

    if( !create_linkage )
        {
            str_stream << "--linked              " << linked_fname  << "\n" << 
                          " --threshold           " << threshold <<  "\n" << 
                          " --enriched            " << enriched_fname <<  "\n" << 
                          " --single_threaded     " <<  bool_str( single_threaded ) <<  "\n" << 
                          " --fractional_scoring  " <<  bool_str( fractional_scoring ) <<  "\n" << 
                          " --summation_scoring   " <<  bool_str( summation_scoring ) <<  "\n" << 
                          " --id_name_map         " << id_name_map_fname << "\n" << 
                          " --output              " << output_fname <<  "\n"
                  ;
        }

    else
        {
            str_stream << "--create_linkage      " << bool_str( create_linkage )  << "\n" << 
                          " --protein_file        " << prot_file_fname <<  "\n"
                          " --peptide_file        " << peptide_file_fname <<  "\n"
                          " --output              " << output_fname <<  "\n"
                          " --kmer_size           " << k <<  "\n"
                ;
        }
    str_stream << "\n";


    return str_stream.str();
}
