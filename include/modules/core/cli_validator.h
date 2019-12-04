#ifndef CLI_VALIDATOR_HH_INCLUDED
#define CLI_VALIDATOR_HH_INCLUDED

/**
 * Ensures that a module has been passed 
 * by the cli.
 **/
class cli_validator
{
 public:
    /**
     * Default constructor
     **/
    cli_validator();

    /**
     * Validate a set of commandline args.
     * Defines a description and a usage message if 
     * '-h' or '--help' were included as arguments.
     * @param argc The number of args supplied at the commandline.
     * @param argv A pointer to an array of strings containing the arguments
     *        and their values.
     * @returns true if a valid module name was supplied, 
     *          false if either '-h', or '--help' were supplied.
     * @throws std::runtime_error if an invalid module is specified.          
     **/
    bool validate( int argc, char ***argv );
};



#endif // CLI_VALIDATION_HH_INCLUDED
