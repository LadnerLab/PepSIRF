#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "options_parser_deconv.h"


options_parser_deconv::options_parser_deconv() = default;

bool options_parser_deconv::parse( int argc, char ***argv, options *opts )
{
    options_deconv *opts_deconv = (options_deconv*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework species deconvolution module. \n"
                                  "This module has two different modes: scoring species and creating a linkage file. Each of \n"
                                  "these modes has its own set of arguments and parameters. The description of each module is \n"
                                  "followed by the mode for which this command is, in brackets. For example, if the description is \n"
                                  "followed by [scoring_species], then this argument is for the scoring species mode. Similarly, \n"
                                  "[create_linkage] is followed by linkage creation arguments. Arguments pertinent to both \n"
                                  "modes are followed by [scoring_species,create_linkage].\n"
                                );
    desc.add_options()
        ( "help,h", "Produce help message\n" )
        ( "linked,l", po::value<std::string>( &opts_deconv->linked_fname ),
          "Name of file containing peptide to species linkages. [scoring_species]\n"
        )
        ( "threshold,t", po::value<std::size_t>( &opts_deconv->threshold ),
          "Threshold number of peptides for a species to be considered. [scoring_species]\n"
        )
        ( "output,o", po::value<std::string>( &opts_deconv->output_fname )->default_value( "deconv_output.tsv" ),
          "Name of the file to write output to. Output will be in the form of "
          "a tab-delimited file with a header. Each entry will be of the form:\n"
          "species_id\\tcount\n [create_linkage,scoring_species]\n"
        )
        ( "single_threaded", po::bool_switch( &opts_deconv->single_threaded )->default_value( false ),
          "By default this module uses two threads. Include this option with no arguments if you only want "
          " one thread to be used. [create_linkage,scoring_species]\n"
        )
        ( "fractional_scoring", po::bool_switch( &opts_deconv->fractional_scoring )->default_value( false ),
          "Use fractional instead of integer scoring. For integer scoring the score of each species is "
          "defined by the number of peptides that share a 7mer with that species. For fractional scoring "
          "the score of each species is defined by 1/n for each peptide, where n is the number of species "
          "a peptide shares a 7mer with. In this method of scoring "
          "peptides with fewer species are worth more. Note that if neither this flag nor --summation_scoring "
          "are included, integer scoring will be used. In integer scoring each species is scored by the "
          "number of peptides it shares a kmer with. [scoring_species]\n" 
        )
        ( "summation_scoring", po::bool_switch( &opts_deconv->summation_scoring )->default_value( false ),
          "Include this flag (without any arguments) if you want summation scoring to be used instead of "
          "fractional or integer scoring. For summation scoring, the --linked file passed must be of the "
          "form created by --create_linkage. This means a file of tab-delimited values, one per line. \n"
          "Each line is of the form peptide_name TAB id:score,id:score, and so on. Undefined behavior "
          "will result if input is not in this format. For summation scoring, each species is scored "
          "based on the number of kmers it shares with each peptide with which it shares a kmer.\n "
          "For example, assume a line in the --linked file looks like the following: \n"
          "peptide_1 TAB 123:4,543:8\n"
          "Both species '123' and '543' will receive a score of 4 and 8 respectively."
          "Note that if neither this flag nor --summation_scoring "
          "are included, integer scoring will be used. In integer scoring each species is scored by the "
          "number of peptides it shares a kmer with. [scoring_species]\n" 
        )
        ( "enriched,e", po::value<std::string>( &opts_deconv->enriched_fname ),
          "File containing the names of enriched peptides, one per line. "
          "Each file in this file should have a corresponding entry in the "
          "file provided by the --linked option. [scoring_species]\n"
        )
        ( "score_filtering", po::bool_switch( &opts_deconv->score_filtering )->default_value( false ),
          "Include this flag if you want filtering to be done by the score of each species. Note that score is "
          "determined by the different flags specifying how a species should be scored. This means "
          "that any species whose score falls below --threshold "
          "will be removed from consideration. Note that for integer scoring, both score filtering and count filtering "
          "are the same. If this flag is not included, then any species whose count falls below --threshold will "
          "be removed from consideration. Score filtering is best suited for the summation scoring algorithm. [scoring_species]\n"
        )
        ( "peptide_assignment_map", po::value<std::string>( &opts_deconv->species_peptides_out ),
          "If specified, a map detailing which peptides were assigned to which species will be "
          "written. This map will be a tab-delimited file with the first column peptide names, "
          "and the second column is a comma-separated list of species the peptide was assigned to. "
          "Note that this comma-separated list will only contain multiple values in the event of "
          "a tie. [scoring_species]\n"
        )
        ( "score_tie_threshold", po::value<double>( &opts_deconv->score_tie_threshold )->default_value( 0.00 ),
          "Threshold for two species to be evaluated as a tie. For example, if this value is "
          "set to 4, and species_1 has a count of 5, species_2 has count 4 then their "
          "distance is '1', which is less than 4. And so these species are considered "
          "for evaluating whether tye are tied."
          "Note that this does not necessarily mean that these species are tied, "
          "in order for them to be considered tied they must also share a certain percentage or number of "
          "their peptides, as defined by the score_overlap_threshold and tie evaluation strategy "
          "specified. [scoring_species] \n"
        )
        ( "score_overlap_threshold", po::value<double>( &opts_deconv->score_overlap_threshold )->default_value( 1.0 ),
          "Once two species have been found to be within 'score_tie_threshold' number of peptides "
          "of one another, they are then evaluated as a tie. For a two-way tie where integer tie evaluation is used, "
          "if the species share "
          "more than score_overlap_threshold number of peptides, then they are both reported. An example value "
          "is 10. For ratio tie evaluation, two species must share at leat this amount of peptides with each other. "
          "For example, suppose species 1 shares 0.5 of its peptides with species 2, but species 2 only shares 0.1 "
          "of its peptides with species 1. These two will only be reported together if score_overlap_threshold "
          "<= 0.1. [scoring_species] \n" 
        )
        ( "integer_tie_eval", po::bool_switch( &opts_deconv->integer_tie_eval )->default_value( false ),
          "Include this flag if tie evaluation should be done by comparing integer overlap of peptides. "
          "For example, if this flag is included and --score_overlap_threshold is set to 10, then if two "
          "species are tied and share at least 10 peptides then they will be reported together. "
          "Important note: if summation_scoring is used then a special tie-breaking strategy is used. [scoring_species]\n"
        )
        ( "ratio_tie_eval", po::bool_switch( &opts_deconv->ratio_tie_eval )->default_value( false ),
          "Include this flag if tie evaluation should be done by comparing the ratio of a species "
          "peptides it shares with another. Note that two species must share at least --overlap_tie_threshold "
          "with eachother to be considered tied. "
          "For example, suppose species 1 shares 0.5 of its peptides with species 2, but species 2 only shares 0.1 "
          "of its peptides with species 1. These two will only be reported together if score_overlap_threshold "
          "<= 0.1. "
          "Important note: if summation_scoring is used then a special tie-breaking strategy is used. [scoring_species]\n"
        )
        ( "create_linkage", po::bool_switch( &opts_deconv->create_linkage )->default_value( false ),
          "Boolean switch to create the linkage file that is used as input for "
          "the species deconvolution process. If this option is included "
          "then 'protein_file' and 'peptide_file' must also be included. [create_linkage]\n"
        )
        ( "protein_file", po::value<std::string>( &opts_deconv->prot_file_fname ),
          "Name of fasta file containing protein sequences from which a design was "
          "created. [create_linkage]\n"
        )
        ( "peptide_file", po::value<std::string>( &opts_deconv->peptide_file_fname ),
          "Name of fasta file containing aa peptides that have been designed as part "
          "of a library. [create_linkage]\n"
        )
        ( "kmer_size,k", po::value<std::size_t>( &opts_deconv->k ), "Kmer size to use when creating "
          "the linkage map. Entries in the linkage file will contain peptides and the species ids of "
          "species that share a kmer with this peptide. For example, if k is 7 and there exists a line "
          "in the linkage file of the form: \n 'peptide_1 TAB 455:12,423:10'\n then peptide_1 "
          "shares 12 7-mers with the species with id '455', and 10 7-mers with the species that has id 423. [create_linkage]\n"
        )
        ( "id_name_map", po::value<std::string>( &opts_deconv->id_name_map_fname )->default_value( "" ),
          "File containing mappings from taxonomic id to name. This file should be formatted like the "
          "file 'rankedlineage.dmp' from NCBI. It is recommended to either use this file or a subset of this file "
          "that at least contains the species ids of the designed peptides. If included, the output will contain "
          "a column denoting the name of the species as well as the id. [create_linkage]\n"
        );

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

