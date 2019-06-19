#include "options_parser_demux.h"

options_parser_demux::options_parser_demux() = default;

bool options_parser_demux::parse( int argc, char ***argv, options *opts )
{
    options_demux *opts_demux = (options_demux*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework demultiplexing module. \n"
                                  "This module takes the following parameters and outputs the counts of each reference \n"
                                  "sequence for each sample."
                                );
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
        ( "max_mismatches,m", po::value<std::size_t>( &opts_demux->max_mismatches )->default_value( 0 ),
          "The maximum number of 'mismatches' to tolerate when parsing reads. If a read is not within this value of "
          "any of the sequences within the designed library it will not be considered. Note that here we define "
          "a mismatch by the Hamming distance D between a reference sequence r and a read sequence s. If D( r, s ) "
          "<= max_mismatches we say that r and s are similar. Note that if for some read q D( r, s ) <= max_mismatches and D( r, q ) "
          "<= max_mismatches, then we say that the sequence whose distance is the minimum between D( r, s ) and D( r, q ) maps to "
          "reference r. \n"
        )
        ( "seq_start", po::value<std::size_t>( &opts_demux->seq_start )->required(), "Start index (0-based) of each read where we expect the designed peptide to "
          "begin. For each read, we start at this index and read for seq_len number of characters. Remember: this argument should be zero-based!\n"
        )
        ( "seq_len", po::value<std::size_t>( &opts_demux->seq_len )->required(), "The length of the designed peptides. Note that we assume "
          "all of the designed peptides are the same length.\n"
        )
        ( "f_index_start", po::value<std::size_t>( &opts_demux->f_index_start )->required(), "Start index (0-based) of each read where we expect the forward index sequences to be found.\n" )
        ( "f_index_len", po::value<std::size_t>( &opts_demux->f_index_len )->required(), "Length of forward index sequences. For each read we start at f_index_start and grab f_index_len "
          "nucleotides.\n"
        )
        ( "concatemer", po::value<std::string>( &opts_demux->concatemer ), "Concatenated primer sequences. If this concatemer is found within a read, we know that a potential sequence "
          "from the designed library was not included. The number of times this concatemer is recorded in the input file is reported.\n"
        )
        ( "output,o", po::value<std::string>( &opts_demux->output_fname )->default_value( opts_demux->DEFAULT_OUTPUT_FNAME ), "The name of the output file to write counts to. "
          "Each line in this file will be a tab-separated list of values, where each entry i is either the name of a sequence or the counts for this sequence in "
          "sample i. This file will have a header labelling each column, i'th tab-separated value of column i of the header will be the sample name of sample i. "
          "If we traverse this column, we will see the count of this sample for each sequence. \n"
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
            return false;
        }
    else
        {
            po::notify( vm );
            return true;
        }
}
