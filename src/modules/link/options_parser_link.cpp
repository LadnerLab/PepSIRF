#include "options_parser_link.h"
#include "options_link.h"

#include <boost/program_options.hpp>

options_parser_link::options_parser_link() = default;

bool options_parser_link::parse( int argc, char ***argv, options *opts )
{
    options_link *opts_link = (options_link*) opts;
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
        ( "output,o", po::value<std::string>( &opts_link->output_fname )
          ->default_value( "link_output.tsv" ),
          "Name of the file to which output is written. Output will be in the form of "
          "a tab-delimited file with a header. Each entry/row will be of the form: "
          "peptide_name TAB id:score,id:score, and so on. By default, \"score\" "
          "is defined as the number of shared kmers.\n"
        )
        ( "protein_file", po::value<std::string>( &opts_link->prot_file_fname ),
          "Name of fasta file containing protein sequences of interest.\n"
        )
        ( "peptide_file", po::value<std::string>( &opts_link->peptide_file_fname ),
          "Name of fasta file containing aa peptides of interest. These will generally "
          "be peptides that are contained in a particular assay.\n"
        )
        ( "meta", po::value<std::string>( &opts_link->metadata_fname ),
          "Optional method for providing taxonomic information for each protein "
          "contained in \"--protein_file\". Three comma-separated strings should be provided: "
          "1) name of tab-delimited metadata file, 2) header for column containing protein "
          "sequence name and 3) header for column containing ID to be used in creating the linkage map.\n"
        )
        ( "tax_id_index,t", po::value<std::size_t>( &opts_link->id_index )->default_value( 1 )
          ->notifier( [&]( const std::size_t val ) {
                  if( val > 3 )
                      {
                          throw boost::program_options::invalid_option_value( "tax_id_index must be one of either 0, "
                                                                             "1, 2, or 3."
                                                                           );

                      }
              } ),
          "Optional method for parsing taxonomic information from the names of the protein sequences. "
          "This is an alternative to the \"--meta\" argument and will only work if the protein names contain "
          "\"OXX\" tags of the form \"OXX=variableID,speciesID,genusID,familyID\". If used, a signle "
          "integer should be provided that corresponds to the index (0-based, valid values include 0-3) of "
          "the tax id to use for creating a linkage map. For example, if this argument is passed with the value 1, \n"
          "the species ID will be used (2 for genus, 3 for family. 0 can vary depending upon the \n"
          "method used for assigning the 0'th ID).\n"
        )
        ( "kmer_redundancy_control,r", po::bool_switch( &opts_link->penalize_kmers )->default_value( false ),
          "Optional modification to the way scores are calculated. If this flag is used, then instead "
          "of a peptide receiving one point for each kmer it shares with proteins of a given taxonomic "
          "group, it receives 1 / ( the number of times the kmer appears in the provided peptides ) points.\n"
        )
        ( "kmer_size,k", po::value<std::size_t>( &opts_link->k ), "Kmer size to use when creating "
          "the linkage map.\n"
        )
        (
         "logfile", po::value( &opts_link->logfile )
         ->default_value( "" ),
          "Designated file to which the module's processes are logged. By "
          "default, the logfile's name will include the module's name and the "
          "time the module started running.\n"
        )
        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

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
