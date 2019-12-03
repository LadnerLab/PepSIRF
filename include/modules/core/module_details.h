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

    /**
     * Initialize a module based upon a module name. 
     * If a valid module name is supplied, this will be populated 
     * with the module itself, the module's options, and options_parser.
     * Each of these items can be retrieved after creation.
     * @param module_name The name of the module to initialize
     * @throws std::runtime_error if an invalid module name is provided.
     **/
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
