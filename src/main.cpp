#include <iostream>
#include <exception>
#include <stdlib.h>

#include "options_parser.h"
#include "options.h"
#include "module.h"
#include "cli_validator.h"
#include "module_initializer.h"
#include "pepsirf_version.h"

int main( int argc, char **argv )
{
    // create objects to be used
    module_initializer init;
    cli_validator cli_val;

    bool help_msg_only = false;
	int ret_code = EXIT_SUCCESS;
    const std::string version_no = PEPSIRF_VERSION;

    try
        {
            // create our option parser for parsing command-line options
            help_msg_only = cli_val.validate( argc, &argv );

            // the '-h' argument was passed without a module, help message has been displayed
            if( help_msg_only )
                {
                    return EXIT_SUCCESS;
                }
            // create the correct options object
            init.initialize( argv[ 1 ] );

            // parse the arguments, any incorrect arguments will raise an error
            // note that parser->parse returns true if arguments other than
            // help are passed
            // help_msg_only = !parser->parse( argc, &argv, opts );
            help_msg_only = !init.get_options_parser()
                             ->set_version( version_no )
                             ->parse( argc, &argv, init.get_opts() );

            if( !help_msg_only )
                {
				    // run PepSIRF with options parsed from command-line
                    std::cout << "PepSIRF (v"
                              << version_no << ") \n"
                              << "Starting module "
                              << init.get_module()->get_name()
                              << " with arguments: \n "
                              << init.get_opts()->get_arguments();

                    init.get_module()->run( init.get_opts() );
                }
        }
    catch( std::exception& e )
        {
            // command-line errors will be handled here
            std::cerr << "Error: " << e.what() << std::endl;

            ret_code = EXIT_FAILURE;
        }

	return ret_code;
}
