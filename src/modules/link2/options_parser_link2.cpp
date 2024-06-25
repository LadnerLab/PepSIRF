#include "logger.h"
#include "options_parser_link2.h"
#include "options_link2.h"

#include <boost/program_options.hpp>

options_parser_link2::options_parser_link2() = default;

bool options_parser_link2::parse( int argc, char ***argv, options *opts )
{
    options_link2 *opts_link2 = (options_link2*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune "
                                  "Response Framework link module.\n",
                                  line_width
                                  );

    desc.add_options()
        (
         "help,h", "Produce help message and exit.\n"
         "The link module is used to create the \"--linked\" input file for the deconv module. "
         "The output file from this module defines linkages between taxonomic "
         "groups (or other groups of interest) and peptides based on shared kmers.\n"
        )
        // TODO: describe format of linkage file
        ( "output,o", po::value<std::string>( &opts_link2->output_fname )
          ->default_value( "link_output.tsv" ),
          "Name of the file to which output is written. Output will be in the form of "
          "a tab-delimited file with a header. Each entry/row will be of the form: "
          "peptide_name TAB id:score,id:score, and so on. By default, \"score\" "
          "is defined as the number of shared kmers.\n"
        )
        ( "protein_file", po::value<std::string>( &opts_link2->prot_file_fname ),
          "Name of fasta file containing protein sequences of interest.\n"
        )
        ( "peptide_file", po::value<std::string>( &opts_link2->peptide_file_fname ),
          "Name of fasta file containing aa peptides of interest. These will generally "
          "be peptides that are contained in a particular assay.\n"
        )
        ( "meta", po::value<std::string>( &opts_link2->metadata_info ),
          "Taxonomic information for each protein "
          "contained in \"--protein_file\". Three comma-separated strings should be provided: "
          "1) name of tab-delimited metadata file, 2) header for column containing protein "
          "sequence name and 3) header for column containing ID to be used in creating the linkage map.\n"
        )
        ( "kmer_size,k", po::value<std::size_t>( &opts_link2->kmer_size ), "Kmer size to use when creating "
          "the linkage map.\n"
        )
        ( "span,s", po::value<std::size_t>( &opts_link2->span ), "Number of positions across which the match can span.\n"
        )
        ("logfile", po::value( &opts_link2->logfile )
         ->default_value( options_link2::set_default_log() ),
         "Designated file to which the module's processes are logged. By "
         "default, the logfile's name will include the module's name and the "
         "time the module started running.\n"
        )
        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if (
        vm.count( "help" )
        || argc == 2
    ) {
        std::ostringstream info_str;
        info_str << desc << "\n";

        Log::info(info_str.str());

        return false;
    }
    else
    {
        po::notify( vm );
        return true;
    }
}
