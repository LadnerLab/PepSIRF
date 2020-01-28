#ifndef NAN_REPORT_HH_INCLUDED
#define NAN_REPORT_HH_INCLUDED
#include <string>
#include <cstdint>
#include <iostream>

/**
 * A 'nan_report' contains information about a 
 * peptide whose zscore is nan. Each report tracks 
 * the name of the probe with a 'nan' zscore, 
 * the name of the sample that probe is in, and the 
 * number of the bin that probe is in.
 **/
struct nan_report
{
public:

    /**
     * The name of the probe with a 
     * zscore of 'nan'
     **/
    std::string probe_name;

    /**
     * The name of the sample the probe 
     * is in.
     **/
    std::string sample_name;

    /**
     * The bin number of the probe.
     **/
    std::size_t bin_number;

    /**
     * Character to join on when 
     * writing representing the report 
     * as a string.
     **/
    char join_char;


    nan_report( const std::string& probe_name,
                const std::string& sample_name,
                const std::size_t bin_number,
                const char join_char
                ) : probe_name{ probe_name },
                    sample_name{ sample_name },
                    bin_number{ bin_number },
                    join_char{ join_char }
    {}

    nan_report( const std::string& probe_name,
                const std::string& sample_name,
                const std::size_t bin_number
                ) : nan_report{ probe_name, sample_name, bin_number, '\t' }
    {}

    /**
     * Return the string representation of this 
     * report, joined on this->join_char
     **/
    std::string to_string() const
    {
        return probe_name + join_char
            + sample_name + join_char
            + std::to_string( bin_number);
    }

};

/**
 * Write the nan_report to a stream.
 * @tparam Stream the type of strime to write the report to.
 * @param stream the stream to write the report to
 * @param rep The report to write to the stream.
 **/
template<typename Stream>
Stream& operator<<( Stream& stream, const nan_report& rep )
{
    stream << rep.to_string();
    return stream;
}

#endif // NAN_REPORT_HH_INCLUDED
