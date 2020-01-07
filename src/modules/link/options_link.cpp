#include <sstream>

#include "options_link.h"

options_link::options_link() = default;

std::string options_link::get_arguments()
{
    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

            str_stream << 
                          "--protein_file                   " << prot_file_fname <<  "\n" <<
                          " --peptide_file                   " << peptide_file_fname <<  "\n" <<
                          " --tax_id_index                   " << id_index << "\n" <<
                          " --kmer_redundancy_control        " << bool_str( penalize_kmers )  << "\n" <<
                          " --output                         " << output_fname <<  "\n" <<
                          " --kmer_size                      " << k <<  "\n"
                                                               << "\n"
                              ;
    return str_stream.str();
}
