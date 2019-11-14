#include <iostream>
#include <exception>
#include <stdlib.h>

#include "options_parser_factory.h"
#include "options_factory.h"
#include "options_parser.h"
#include "options.h"
#include "module.h"
#include "module_factory.h"

int main( int argc, char **argv )
{
    // create objects to be used
    options_parser_factory parser_fact;
    options_parser *parser = nullptr;

    options_factory opts_fact;
    options *opts = nullptr;

    module_factory mod_fact;
    module *mod = nullptr;

    bool help_msg_only = false;

	int ret_code = EXIT_SUCCESS;

    try
        {
            // create our option parser for parsing command-line options
            parser = parser_fact.create( argc, &argv );

            // the '-h' argument was passed without a module, help message has been displayed
            if( !parser )
                {
                    return EXIT_SUCCESS;
                }
            // create the correct options object
            opts = opts_fact.create( argc, &argv );

            // parse the arguments, any incorrect arguments will raise an error
            // note that parser->parse returns true if arguments other than
            // help are passed
            help_msg_only = !parser->parse( argc, &argv, opts );

            if( !help_msg_only )
                {
				    // run PepSIRF with options parsed from command-line
                    mod = mod_fact.create( argv[ 1 ] );
                    std::cout << "Starting module " << mod->get_name() << " with arguments: \n " << opts->get_arguments();
                    mod->run( opts );
                }
        }
    catch( std::exception& e )
        {
            // command-line errors will be handled here
            std::cerr << "Error: " << e.what() << std::endl;

            ret_code = EXIT_FAILURE;
        }

    delete parser;
    delete mod;
	delete opts;
	return ret_code;
}
