#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "fs_tools.h"
#include "logger.h"
#include "options_parser_deconv.h"
#include "util.h"


options_parser_deconv::options_parser_deconv() = default;

bool options_parser_deconv::parse( int argc, char ***argv, options *opts )
{
    options_deconv *opts_deconv = (options_deconv*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response "
                                  "Framework species deconvolution module.\n",
                                  line_width
                                );
    desc.add_options()
        ( "help,h", "Produce help message\n"
          "The deconv module converts a list of enriched peptides into a parsimony-based "
          "list of likely taxa to which the assayed individual has likely been exposed. "
          "This module has two modes: batch and singular. In batch mode, the input given "
          "to '--enriched' is a directory containing files of enriched peptides for each "
          "sample (e.g., as output by enrich). In this case, '--output' "
          "is a directory where the output deconvolution reports will be written. In "
          "singular mode, both '--enriched' and '--output' are treated as files, not directories. "
          "The chosen mode is determined by the type of argument provided with '--enriched'. "
          "If a directory is specified, batch mode will be used. If a file is specified, "
          "singular mode will be used.\n"
        )
        ( "linked,l", po::value<std::string>( &opts_deconv->linked_fname ),
          "Name of linkage map to be used for deconvolution. It should be in the format output "
          "by the 'link' module.\n"
        )
        ( "threshold,t", po::value<std::size_t>( &opts_deconv->threshold ),
          "Minimum score that a taxon must obtain in order to be included in the deconvolution report.\n"
        )
        ( "outfile_suffix", po::value( &opts_deconv->outfile_suffix )
          ->default_value( "" ),
          "Used for batch mode only. When specified, the name of each file written to the output "
          " directory will have this suffix.\n"
        )
        ( "output,o", po::value<std::string>( &opts_deconv->output_fname )->default_value( "deconv_output.tsv" ),
          "Name for the output directory or file. Output will be in the form of either one file "
          "or a directory containing files. Each output file will be tab-delimited and will include "
          "a header.\n"
        )
        ( "remove_file_types,r", po::bool_switch( &opts_deconv->remove_file_types )
          ->default_value( false ),
          "Use this flag to exclude input file ('--enrich') extensions from the names of output files. "
          "Not used in singular mode.\n"
        )
        ("scores_per_round,s", po::value<std::string>( &opts_deconv->orig_scores_dname )->default_value( "" )
         ->notifier(
            [&]( const std::string& val )
            {
                // check that a directory was actually supplied
                if( val.length() )
                {
                    if( boost::filesystem::exists( val ) )
                    {
                        Log::error("Directory '" + val + "' already exists!");
                    }
                    else
                    {
                        fs_tools::create_directory(val);
                    }
                }
            }),
          "Optional. Name of directory to write counts/scores to after every round. If included, "
          "the counts and scores for all remaining taxa will be recorded after every round. "
          "Filenames will be written in the format '$dir/round_x', where x is the round number. "
          "The original scores will be written to '$dir/round_0'. A new file will be written to the "
          "directory after each subsequent round. If this flag is included "
          "and the specified directory exists, the program will exit with an error.\n"
        )
        ( "single_threaded", po::bool_switch( &opts_deconv->single_threaded )->default_value( false ),
          "By default this module uses two threads. Include this option with no arguments if you only want "
          "only one thread to be used. \n"
        )
        ( "scoring_strategy", po::value<std::string>( &opts_deconv->scoring_strategy )->default_value( "summation" ),
          "Scoring strategies \"summation\", \"integer\", or \"fraction\" can be specified. "
          "By not including this flag, summation scoring will be used by default.\n"
          "The --linked file passed must be of the form created by the link module. This means a file of "
          "tab-delimited values, one per line.\n"
          "Each line is of the form peptide_name TAB id:score,id:score, and so on. An error will occur"
          "if input is not in this format.\n"
          "For summation scoring, the score assigned to each peptide/ID pair is determined by the \":score\" "
          "portion of the --linked file.\n"
          "For example, assume a line in the --linked file looks like the following: \n"
          "peptide_1 TAB 123:4,543:8\n"
          "The IDs '123' and '543' will receive scores of 4 and 8 respectively.\n"
          "For integer scoring, each ID receives a score of 1 for every enriched peptide to which it is "
          "linked (\":score\" is ignored).\n"
          "For fractional scoring, the score is assigned to each peptide/ID pair is defined by 1/n for "
          "each peptide, where n is the number of IDs to which a peptide is linked. "
          "In this method of scoring peptides, a peptide with fewer linked IDs is worth more points.\n"
        )
        ( "enriched,e", po::value<std::string>( &opts_deconv->enriched_fname ),
          "Name of a directory containing files, or a single file containing the names of enriched peptides, "
          "one per line. Each peptide contained in this file (or files) should have a corresponding entry in the "
          "'--linked' input file.\n"
        )
        ( "score_filtering", po::bool_switch( &opts_deconv->score_filtering )->default_value( false ),
          "Include this option if you want filtering to be done by the score of each taxon, rather than the count of "
          "linked peptides. If used, any taxon with a score below '--threshold' will be removed from consideration, "
          "even if it is the highest scoring taxon. Note that for integer scoring, both score filtering and count "
          "filtering (default) are the same. If this flag is not included, then any species whose count falls below "
          "'--threshold' will be removed from consideration. Score filtering is best suited for the summation scoring "
          "method.\n"
        )
        ( "peptide_assignment_map,p", po::value<std::string>( &opts_deconv->species_peptides_out ),
          "Optional output. If specified, a map detailing which peptides were assigned to which taxa will be "
          "written. If this module is run in batch mode, this will be used as a directory name for the peptide maps to "
          "be stored. Maps will be tab-delimited files with the first column being peptide names; the second column "
          "containing a comma-separated list of taxa to which the peptide was assigned; the third column will be a list "
          "of the taxa with which the peptide originally shared a kmer. "
          "Note that the second column will only contain multiple values in the event of a tie.\n"
        )
        ( "mapfile_suffix", po::value( &opts_deconv->map_suffix )
          ->default_value( "" ),
          "Used for batch mode only. When specified, the name of each '--peptide_assignment_map' will have this suffix.\n"
        )
        ( "score_tie_threshold", po::value<double>( &opts_deconv->score_tie_threshold )->default_value( 0.00 )
          ->notifier( [&]( const double val ) {
                  if( val > 1 && !util::is_integer( val ) )
                      {
                          throw boost::program_options::invalid_option_value( "If score_tie_threshold is not an integer, "
                                                                              "it must be in (0,1)"
                                                                            );
                      }
                      } ),
          "Threshold for two species to be evaluated as a tie. Note that this value can be either an integer "
          "or a ratio that is in (0,1). When provided as an integer this value dictates the difference in score "
          "that is allowed for two taxa to be considered as potentially tied. For example, if this flag is provided with the value "
          "0, then two or more taxa must have the exact same score to be tied. If this flag is provided with "
          "the value 4, then the difference between the scores of two taxa must be no greater than 4 to be considered tied. "
          "For example, if taxon 1 has a score of 5, and taxon 2 has a score anywhere between the integer values in [1,9], "
          "then these species will be considered tied, and their tie will be evaluated as dictated by the specified "
          "'--score_overlap_threshold'. If the argument provided to this flag is in (0, 1), then the score for a taxon "
          "must be at least this proportion of the score for the highest scoring taxon, to trigger a tie. So if species "
          "1 has the highest score with a score of 9, and species 2 has a score of 5, then this flag must be provided with "
          "value >= 4/5 = 0.8 for the species 1 and 2 to be considered tied. Note that any values provided to this flag that are in "
          "the set { x: x >= 1 } - Z, where Z is the set of integers, will result in an error. So 4.45 is not a valid value, "
          "but both 4 and 0.45 are.\n"
        )
        ( "id_name_map", po::value<std::string>( &opts_deconv->id_name_map_fname )->default_value( "" ),
          "Optional file containing mappings from taxonomic ID to taxon name. This file should be formatted like the "
          "file 'rankedlineage.dmp' from NCBI. It is recommended to either use this file or a subset of this file "
          "that contains all of the taxon ids linked to peptides of interest. If included, the output will contain "
          "a column denoting the name of the species as well as the id.\n"
        )
        ("custom_id_name_map",
         po::value<std::string>(&opts_deconv->custom_id_name_map_fname)
            ->default_value(""),
         "Optional file containing mappings from taxonomic IDs to taxon names."
         " The format of this file is dictated by the user. It is recommended"
         " to either use this file or a subset of this file that contains all"
         " of the taxon ids linked to peptides of interest. If included, the"
         " output will contain a column denoting the name of the species as"
         " well as the ID.\n"
        )
        ("score_overlap_threshold",
         po::value<double>(&opts_deconv->score_overlap_threshold)->default_value(1.0)
            ->notifier(
                [&](const double val)
                {
                    if(val > 1 && !util::is_integer(val))
                    {
                        throw boost::program_options::invalid_option_value(
                            "If score_overlap_threshold is not an integer, it"
                            " must be in (0, 1). Otherwise an integer must be"
                            " provided."
                        );
                    }
                }),
         "Once two species have been determined to be tied, according to '--score_tie_threshold', they are "
         "then evaluated as a tie. To use integer tie evaluation, where species must share an integer number of "
         "peptides, not a ratio of their total peptides, provide this argument with a value in the interval [1, inf). "
         "For ratio tie evaluation, which is used when this argument is provided with a value in the "
         "interval (0,1), two taxon must reciprocally share at least the specified proportion of peptides to be "
         "reported together. For example, suppose species 1 shares half (0.5) of its peptides with species 2, but species "
         "2 only shares a tenth (0.1) of its peptides with species 1. These two will only be reported together if "
         "score_overlap_threshold' <= 0.1.\n"
        )
        ( "enriched_file_ending", po::value<std::string>( &opts_deconv->enriched_file_ending )->default_value( "_enriched.txt" ),
          "Optional flag that specifies what string is expected at the end of each file containing enriched peptides. "
          "Set to \"_enriched.txt\" by default \n"
        )
        ("logfile", po::value(&opts_deconv->logfile)
            ->default_value(options_deconv::set_default_log()),
         "Designated file to which the module's processes are logged. By"
         " default, the logfile's name will include the module's name and the"
         " time the module started running.\n"
        )
        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if (vm.count( "help" ) || argc == 2)
    {
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

