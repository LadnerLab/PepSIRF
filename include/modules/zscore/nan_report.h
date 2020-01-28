#ifndef NAN_REPORT_HH_INCLUDED
#define NAN_REPORT_HH_INCLUDED
#include <string>
#include <cstdint>
#include <iostream>

struct nan_report
{
public:
    std::string probe_name;
    std::string sample_name;
    std::size_t bin_number;
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

    std::string to_string() const
    {
        return probe_name + join_char
            + sample_name + join_char
            + std::to_string( bin_number);
    }

};

template<typename Stream>
Stream& operator<<( Stream& stream, const nan_report& rep )
{
    stream << rep.to_string();
    return stream;
}

#endif // NAN_REPORT_HH_INCLUDED
