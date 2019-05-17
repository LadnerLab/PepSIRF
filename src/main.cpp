#include <iostream>
#include <omp.h>
#include <exception>
#include <stdlib.h>

#include "options_parser.h"
#include "options.h"

int main( int argc, char **argv )
{
    // create objects to be used
    options_parser parser;
    options opts;

    try
        {
            // create our option parser for parsing command-line options
            parser = options_parser();
            // parse the arguments, any incorrect arguments will raise an error
            opts = parser.parse( argc, &argv );
        }
    catch( std::exception& e )
        {
            // command-line errors will be handled here
            std::cerr << "Error: " << e.what() << std::endl;

            return EXIT_FAILURE;
        }

}
