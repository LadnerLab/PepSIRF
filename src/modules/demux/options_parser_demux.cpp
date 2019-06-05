#include "options_parser_demux.h"

options_parser_demux::options_parser_demux() = default;

void options_parser_demux::parse( int argc, char ***argv, options *opts )
{
    options_demux *opts_demux = (options_demux*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "input_r1", po::value<std::string>( &opts_demux->input_r1_fname )->required(), "Input forward reads fastq file to parse.\n")
        ( "input_r2", po::value<std::string>( &opts_demux->input_r2_fname ), "Input reverse reads fastq file to parse. Note that if this argument is "
          "not supplied only forward indices will be used to identify samples.\n"
        )
        ( "index", po::value<std::string>( &opts_demux->index_fname )->required(), "Name of fasta file containing forward and (if included )reverse index sequences.\n")
        ( "library,l", po::value<std::string>( &opts_demux->library_fname )->required(), "Designed library containing nucleic acid peptides. "
                                                             "Library should be in fasta format and should contain "
                                                             "sequences that were used to design input_r1.\n"
        )
        ( "read_per_loop,r", po::value<long int>( &opts_demux->read_per_loop )->default_value( opts_demux->DEFAULT_READ_PER_LOOP ), "The number of fastq "
          "records read a time. A higher value will result in more memory usage by the program, but will also result in fewer disk accesses, "
          "increasing performance of the program.\n"
        )
        ( "samplelist,s", po::value<std::string>( &opts_demux->samplelist_fname )->required(), "A tab-delimited list of samples, one sample per line. If the samples are "
          "already indexed by I2 only the forward index (I1) and the sample name are required. The first item in each tab-delimited line is the forward (I1) index, the second "
          "(if included) is the reverse (I2) index, and the third is the samplename. \n"
        )
        ( "num_threads,t", po::value<int>( &opts_demux->num_threads )->default_value( opts_demux->DEFAULT_NUM_THREADS ), "Number of threads to use for analyses.\n" );

    po::store( po::command_line_parser( argc, *argv ).options( desc ).allow_unregistered().run(), vm);

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }
    else
        {
            po::notify( vm );
        }
}
