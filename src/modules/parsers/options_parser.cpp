#include "options_parser.h"
#include <sys/ioctl.h>

options_parser::options_parser() 
{

    #ifdef TIOCGSIZE

    struct ttysize ts;
    ioctl( STDIN_FILENO, TIOCGSIZE, &ts );
    line_width = ts.ts_cols;

    #elif defined( TIOCGWINSZ )

    struct winsize ts;
    ioctl( STDIN_FILENO, TIOCGWINSZ, &ts );
    line_width = ts.ws_col;

    #endif 


}

options_parser *options_parser::set_version( const std::string version_no )
{
    this->version_no = version_no;
    return this;
}

options_parser::~options_parser() = default;

bool options_parser::parse( int argc, char ***argv, options *opts )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response Framework"
                                );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "num_threads,t", po::value<int>( &opts->num_threads )->default_value( opts->DEFAULT_NUM_THREADS ), "Number of threads to use for analyses." );

    po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
            return false;
        }

    po::notify( vm );
    return true;
}


std::string options_parser::format_version_string()
{
    return "(v" + version_no + ")";
}
