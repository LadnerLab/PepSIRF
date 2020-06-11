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
        )
        // TODO: describe format of linkage file
        ( "output,o", po::value<std::string>( &opts_link->output_fname )
          ->default_value( "link_output.tsv" ),
          "Name of the file to write output to. Output will be in the form of "
          "a tab-delimited file with a header. Each entry will be of the form:\n"
          "peptide_name TAB id:score,id:score, and so on. \n"
        )
        ( "protein_file", po::value<std::string>( &opts_link->prot_file_fname ),
          "Name of fasta file containing protein sequences from which a design was "
          "created.\n"
        )
        ( "peptide_file", po::value<std::string>( &opts_link->peptide_file_fname ),
          "Name of fasta file containing aa peptides that have been designed as part "
          "of a library.\n"
        )
        ( "meta", po::value<std::string>( &opts_link->metadata_fname ),
          "Name of metadata file with \".metadata\" extension, protein sequence name and species identification name. "
          "The three entries should be comma delimited. The protein sequence and species identification name should be "
          "consistent with their metadata file header column names."
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
          "The index (0-based, valid values include 0-3) of the tax id to use for "
          "linking peptides and species. For example, if this argument is passed with the value 1, \n"
          "the species ID will be used. (2 for genus, 3 for family. 0 can vary depending upon the \n"
          "method used for assigning the 0'th ID.\n"
        )
        ( "kmer_redundancy_control,r", po::bool_switch( &opts_link->penalize_kmers )->default_value( false ),
          "Control for kmer redundancy when creating the peptide linkage map. Instead of a peptide receiving "
          "one point for each kmer it receives for a species, it recieves 1 / ( the number of times the kmer "
          "appears in the original design ) points.\n"
        )
        ( "kmer_size,k", po::value<std::size_t>( &opts_link->k ), "Kmer size to use when creating "
          "the linkage map. Entries in the linkage file will contain peptides and the species ids of "
          "species that share a kmer with this peptide. For example, if k is 7 and there exists a line "
          "in the linkage file of the form: \n 'peptide_1 TAB 455:12,423:10'\n then peptide_1 "
          "shares 12 7-mers with the species with id '455', "
          "and 10 7-mers with the species that has id 423.\n"
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
