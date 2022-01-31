#ifndef OPTIONS_PARSER_QC_HH_INCLUDED
#define OPTIONS_PASRER_QC_HH_INCLUDED
#include "options_parser.h"
#include "options.h"

/**
 * Options parser for the qc module
 * 
 */
class options_parser_qc : public options_parser
{
    public:

        options_parser_qc();

        bool parse( int argc, char ***argv, options *opts );
};

#endif