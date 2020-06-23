#include "options_parser_demux.h"
#include <stdexcept>

options_parser_demux::options_parser_demux() = default;

bool options_parser_demux::parse( int argc, char ***argv, options *opts )
{
    options_demux *opts_demux = (options_demux*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response Framework demultiplexing module. \n"
                                  "This module takes the following parameters and outputs counts for each reference \n"
                                  "sequence (i.e. probe/peptide) for each sample. For this module we define 'distance' "
                                  "as the Hamming distance D between a reference sequence r and a read sequence s. If D( r, s ) "
                                  "<= max_mismatches we say that s maps to r. Note that if multiple reference sequences, r1 and r2, "
                                  "are similar to s, D( r1, s ) <= max_mismatches and D( r2, s ) <= max_mismatches, "
                                  "then we say that s maps to the reference with the minimum distance. Additionally if D( r1, s ) "
                                  "== D( r2, q ) then we discard s as we cannot say whether s maps to r1 or r2. \n",
                                  line_width
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit" )
        ( "input_r1",
          po::value<std::string>( &opts_demux->input_r1_fname )->required(),
          "Fastq-formatted file containing reads with DNA tags. "
          "If PepSIRF was NOT compiled with Zlib support, this file must be uncompressed. If PepSIRF was compiled with Zlib support, "
          "then this file can be uncompressed or compressed using gzip. In this case, the file format will be automatically determined.\n"
          )
        ( "input_r2",
          po::value<std::string>( &opts_demux->input_r2_fname ),
          "Optional index-only fastq file. "
          "If PepSIRF was NOT compiled with Zlib support, this file must be uncompressed. If PepSIRF was compiled with Zlib support, "
          "then this file can be uncompressed or compressed using gzip. In this case, the file format will be automatically determined.\n"
          "Note that if this argument is not supplied, only \"index1\" will be used to identify samples.\n"
        )
        ( "index,i", po::value<std::string>( &opts_demux->index_fname )->required(), "Name of fasta-formatted file containing forward and "
                                                                                     "(potentially) reverse index sequences. Sequence names must match "
                                                                                     "exactly with those supplied in the \"samplelist\"\n")
        ( "library,l", po::value<std::string>( &opts_demux->library_fname )->default_value( "" ), "Fasta-formatted file containing reference DNA tags. "
          "If this flag is not included, reference-independent demultiplexing will be performed. "
          "In reference-independent mode, each sequence in the region specified by '--seq' will be considered its own reference, "
          "and the observed sequences will be used as the row names in the output count matrix.\n"
        )
        ( "read_per_loop,r", po::value<long int>( &opts_demux->read_per_loop )->default_value( opts_demux->DEFAULT_READ_PER_LOOP ), "The number of fastq "
          "records read a time. A higher value will result in more memory usage by the program, but will also result in fewer disk accesses, "
          "increasing performance of the program.\n"
        )
        ( "index1", po::value<std::string>()
                     ->required()
                     ->notifier( [&]( const std::string &vals ) {
                             opts_demux->set_info( &options_demux::f_index_data,
                                                    vals
                                                 );
                                                                }
                               ),
          "Positional information for index1 (i.e barcode 1). This argument must be passed as 3 comma-separated values. The first item represents the (0-based) expected "
          "start position of the first index; the second represents the length of the first index; and the third represents the number "
          "of mismatches that are tolerated for this index. An example is '--index1 12,12,1'. This says that the index starts at (0-based) position "
          "12, the index is 12 nucleotides long, and if a perfect match is not found, then up to one mismatch will be tolerated.\n."
        )
        ( "seq", po::value<std::string>()
                     ->required()
                     ->notifier( [&]( const std::string &vals ) {
                             opts_demux->set_info( &options_demux::seq_data,
                                                    vals
                                                 );
                                                                }
                               ),
          "Positional information for the DNA tags. This argument must be passed in the same format specified for \"index1\".\n"
        )
        ( "index2", po::value<std::string>()
                     ->default_value( "0,0,0" )
                     ->notifier( [&]( const std::string &vals ) {
                             opts_demux->set_info( &options_demux::r_index_data,
                                                    vals
                                                 );
                                                                }
                               ),
          "Positional information for index2, optional. This argument must be passed in the same format specified for \"--index1\". "
          "If \"--input2\" is provided, this positional information is assummed to refer to the reads contained in this second, index-only fastq file. "
          "If \"--input_r2\" is NOT provided, this positional information is assumed to refer to the reads contained in the \"--input_r1\" fastq file.\n"
        )
        ( "concatemer,c", po::value<std::string>( &opts_demux->concatemer ),
        "Concatenated adapter/primer sequences (optional). The presence of this sequence within a read indicates that the "
        "expected DNA tag is not present. If supplied, the number of times this concatemer is recorded in the input file is reported.\n"
        )
        ( "output,o", po::value<std::string>( &opts_demux->output_fname )->default_value( opts_demux->DEFAULT_OUTPUT_FNAME ), "Name for the output counts file. "
          "This output file will be tab-delimited and will contain a header row. The first column will contain probe/peptide names. Each subsequent column will "
          "contain probe/peptide counts for a sample, with one column per sample.\n"
        )
        ( "aa_counts,a", po::value<std::string>( &opts_demux->aggregate_fname )->default_value( "" ),
          "Name for an output file that will contain aggregated aa-level counts. This is relevant when peptides from a designed library "
          "have multiple different nt-level encodings. "
          "If this option is included without the '--translate_aggregates' flag, names of sequences in the file supplied "
          "by the \"--library flag\" MUST be of the form ID-NUM, where ID can contain any characters except '-', and NUM "
          "represents the id of this encoding. ID and NUM MUST be separated by a single dash '-' character. For example, suppose we have TG1_1-1 and TG1_1-2 "
          "in our library, which says that we generated two encodings for the TG1_1 peptide. The \"--aa_counts\" file will have a single TG1_1 "
          "entry, with per sample counts that are the sum of the counts from TG1_1-1 and TG1_1-2.\n"

        )
        ( "translate_aggregates", po::bool_switch( &opts_demux->translation_aggregation )->default_value( false ),
          "Include this flag to use translation-based aggregation. "
          "In this mode, counts for nt sequences will be combined if they translate into the same aa sequence. "
          "Note: When this mode is used, the name of the aggregate sequence will be the sequence that was a result of the translation. "
          "Therefore, this mode is most appropriate for use with reference-independent demultiplexing.\n"
        )
        ( "samplelist,s", po::value<std::string>( &opts_demux->samplelist_fname )->required(), "A tab-delimited list of samples, one sample per line. "
          "The first item in each tab-delimited line is \"--index1\", the second (if included) is the \"--index2\", and the third is the samplename. "
          "\"--index2\" can be omitted if only using a single index.\n"
        )
        ( "phred_base", po::value<int>( &opts_demux->phred_base )->default_value( 33 )
          ->notifier( []( const int value ){
                      if( !( value == 33 || value == 64 ) )
                          { throw std::runtime_error( "Phred values can only be 33 or 64" ); }
                                            }
                    ),
          "Phred base to use when parsing fastq quality scores. Valid options include 33 or 64.\n"
        )
        ( "phred_min_score", po::value<int>( &opts_demux->min_phred_score )->default_value( 0 ),
          "The minimum average phred-scaled quality score for the DNA tag "
          "portion of a read for it to be considered for matching. This means that if the average "
          "phred33/64 score for a read at the expected locations of the DNA tag is not at least "
          "this then the read will be discarded.\n"
        )
        ( "num_threads,t", po::value<int>( &opts_demux->num_threads )->default_value( opts_demux->DEFAULT_NUM_THREADS ), "Number of threads to use for analyses.\n" );


    po::store( po::command_line_parser( argc, *argv ).options( desc ).allow_unregistered().run(), vm);

    if( vm.count( "help" )
	    || argc == 2
	  )
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
