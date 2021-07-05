#ifndef MODULE_INITIALIZER_HH_INCLUDED
#define MODULE_INITIALIZER_HH_INCLUDED
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
#include "options_zscore.h"
#include "options_bin.h"
#include "options_link.h"
#include "options_enrich.h"
#include "options_info.h"

// options parsers for each module
#include "options_parser_deconv.h"
#include "options_parser_demux.h"
#include "options_parser_normalize.h"
#include "options_parser_subjoin.h"
#include "options_parser_zscore.h"
#include "options_parser_bin.h"
#include "options_parser_link.h"
#include "options_parser_enrich.h"
#include "options_parser_info.h"

// the modules themselves
#include "module_deconv.h"
#include "module_demux.h"
#include "module_normalize.h"
#include "module_subjoin.h"
#include "module_zscore.h"
#include "module_bin.h"
#include "module_link.h"
#include "module_enrich.h"
#include "module_info.h"


/**
 * Contains the items necessary to instantiate a module.
 * to add a module, add the entry to module_initializer::initialize.
 **/
class module_initializer
{
 public:

    /**
     * Construct the module_initializer, initialize everything to null.
     **/
    module_initializer();

    /**
     * Destroy the initializer, deletes any allocated memory.
     **/
    ~module_initializer();

    /**
     * Initialize a module based upon a module name. 
     * If a valid module name is supplied, this will be populated 
     * with the module itself, the module's options, and options_parser.
     * Each of these items can be retrieved after creation.
     * @param module_name The name of the module to initialize
     * @throws std::runtime_error if an invalid module name is provided.
     **/
    void initialize( const std::string& module_name );

    /**
     * Get the initialized options for the 
     * specified module. 
     * @pre This has been successfully initialized
     * @returns A pointer to the options initialized by this module
     * @throws std::runtime_error if this has not yet been initialized.
     **/
    options *get_opts();

    /**
     * Get the initialized module object for the 
     * specified module. 
     * @pre This has been successfully initialized
     * @returns A pointer to the module object initialized by this module
     * @throws std::runtime_error if this has not yet been initialized.
     **/
    module *get_module();

    /**
     * Get the initialized options_parser for the 
     * specified module. 
     * @pre This has been successfully initialized
     * @returns A pointer to the options_parser initialized by this module
     * @throws std::runtime_error if this has not yet been initialized.
     **/
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

    /**
     * Throws std::runtime_error if ptr == nullptr
     **/
    template<typename T>
    void throw_on_null( const T *ptr ) const
        {
            if( ptr == nullptr )
                {
                    throw std::runtime_error(
                                             "The provided pointer is NULL, "
                                             "operations cannot be performed on it."
                                            );
                }
        }

};


#endif // MODULE_INITIALIZER_HH_INCLUDED
