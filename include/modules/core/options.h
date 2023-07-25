#ifndef OPTIONS_HH_INCLUDED
#define OPTIONS_HH_INCLUDED
#include <string>

/*! Data class to contain and handle 
 * arguments that will be passed in from the 
 * command-line.
*/
class options
{
public:
    options(): DEFAULT_NUM_THREADS( 2 ){} //!< Default constructor.
    virtual ~options();

    // name of file to which runtime logs will be written
    static std::string logfile;

    int num_threads; //!< The number of threads to use for computation.
    const int DEFAULT_NUM_THREADS; //!< The default number of threads to use

    /**
     * Get the arguments that were supplied to the options 
     * class by the command-line.
     * @returns String containing arguments, one per line. Arguments are formatted in 
     *          '--arg_name argument' format.
     */
    virtual std::string get_arguments();

    /**
     * Sets the logfile name to the default dynamic name
     * @return dynamic name which includes the module name and time of run
     */
    virtual std::string set_default_log();

    /**
     * Gets the module's log file
     * @return module's log file name
     */
    static std::string get_logfile();
};


#endif /* OPTIONS_HH_INCLUDED */

