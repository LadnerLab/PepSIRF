#ifndef MODULE_DETAILS_HH_INCLUDED
#define MODULE_DETAILS_HH_INCLUDED
#include <string>
#include <stdexcept>

#include "options.h"
#include "module.h"
#include "options_parser.h"

// options for each module
#include "options_deconv.h"
#include "options_demux.h"
#include "options_normalize.h"
#include "options_subjoin.h"

// options parsers for each module
#include "options_parser_deconv.h"
#include "options_parser_demux.h"
#include "options_parser_normalize.h"
#include "options_parser_subjoin.h"

// the modules themselves
#include "module_deconv.h"
#include "module_demux.h"
#include "module_normalize.h"
#include "module_subjoin.h"


/**
 * Contains the items necessary to instantiate a module.
 * To add a module, add the entry to module_details::initialize.
 **/
class module_details
{
 public:

    /**
     * Construct the module_details, initialize everything to null.
     **/
    module_details();

    /**
     * Destroy the initializer, deletes any allocated memory.
     **/
    ~module_details();
    void initialize( const std::string& module_name );

    options *get_opts();
    module *get_module();
    options_parser *get_options_parser();

 private:

    /**
     * A pointer to a subclass of options.
     **/
    options *opts;

    /**
     * A pointer to a subclass of module.
     **/
    module *mod;

    /**
     * A pointer to a subclass of options_parser
     **/
    options_parser *opt_parser;


};


#endif // MODULE_DETAILS_HH_INCLUDED
