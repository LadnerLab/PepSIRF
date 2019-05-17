#include <iostream>
#include <omp.h>
#include <exception>
#include <stdlib.h>

#include "options_parser.h"
#include "options.h"

int main( int argc, char **argv )
{
    options_parser parser;
    options opts;

    try
        {
            parser = options_parser();
        }
    catch( std::exception& e )
        {
            std::cerr << "Error: " << e.what() << std::endl;

            return EXIT_FAILURE;
        }

    opts = parser.parse( argc, &argv );
}
