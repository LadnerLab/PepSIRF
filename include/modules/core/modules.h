#ifndef MODULES_HH_DEFINED
#define MODULES_HH_DEFINED
#include <unordered_set>
#include <iostream>

/**
 * Used to track information for module names, 
 * and tools for operating upon those names.
 **/
namespace modules
{
    /**
     * A set of module names, each with a corresponding module.
     **/
    const std::unordered_set<std::string> module_names = { "demux", "deconv",
                                                           "norm", "subjoin",
                                                           "zscore"
                                                         };

    // a 'StringLike' is something that is string-constructible, i.e.
    // std::string( StringLike a ) works.
    template<typename StringLike>
    bool is_module( const StringLike& name )
    {
        return module_names.find( std::string( name ) ) != module_names.end();
    }

    /**
     * Join the current set of modules on a string.
     * @param on The string to join each of the modules on,
     *        so 'on' is placed in between concatenated module names.
     * @note The last module name will also have the suffix 'on'
     * @returns a string of the form:
     *          mod1onmod2onmod3on
     **/
    std::string join_names( const std::string& on );

    /**
     * Join names with both a prefix and a suffix.
     * @param prefix The string to prefix each module name with.
     * @param on The string to use as a suffix for each joined name. 
     * @returns a string containing the names with prefix prepended, separated
     *          by on.
     **/
    std::string join_names( const std::string& prefix,
                            const std::string& on
                          );
    /**
     * Join all of the names found in an iterable 'on' a string.
     * @param names An iterable that is iterable by a range-based 
     *        for loop. This iterable will be used as the names that are
     *        joined on 'on'.
     * @param on The string to join each of the modules on,
     *        so 'on' is placed in between concatenated module names.
     * @note The last module name will also have the suffix 'on'
     * @returns a string of the form:
     *          mod1onmod2onmod3on
     **/
    template<typename Iterable>
        std::string join_names_from_iterable( const Iterable& names,
                                              const std::string& on
                                            )
        {
            std::string out( "" );
            for( const auto& x: names )
                {
                    out += x + on;
                }
            return out;
        }



};

#endif // MODULES_HH_DEFINED
