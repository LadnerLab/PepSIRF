#include "module_demux.h"

module_demux::module_demux() 
{
    name = "Demux";
}

void module_demux::run( options *opts )
{
    options_demux *d_opts = (options_demux*) opts;
    std::vector<sequence> library_seqs;

    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;

    library_seqs = fasta_p.parse( d_opts->library_fname );

    std::cout << "Number of input seqs: " << library_seqs.size() << "\n";


    // 

}


std::string module_demux::get_name()
{
    return name;
}
