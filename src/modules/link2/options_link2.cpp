#include <sstream>

#include "options_link2.h"

options_link2::options_link2() = default;

std::string options_link2::get_arguments()
{
    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

            str_stream <<
                          " --protein_file                   " << prot_file_fname <<  "\n" <<
                          " --peptide_file                   " << peptide_file_fname <<  "\n" <<
                          " --meta                           " << metadata_info << "\n" <<
                          " --kmer_size                      " << kmer_size <<  "\n" <<
                          " --span                           " << span <<  "\n" <<
                          " --num_threads                    " << num_threads << "\n"
                          " --output                         " << output_fname <<  "\n" <<
                          " --logfile                        " << logfile << "\n"
                                                               << "\n"
                              ;
    return str_stream.str();
}