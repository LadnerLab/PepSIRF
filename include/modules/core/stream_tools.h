#ifndef STREAM_TOOLS_HH_INCLUDED
#define STREAM_TOOLS_HH_INCLUDED 

#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

namespace std
{
    // used to convert string '3.4,5.6' -> pair( 3.4, 5.6 )
    template<typename Stream>
    Stream& operator>>( Stream& in, std::pair<double,double>& pair )
    {
        string temp_val;
        in >> temp_val;

        vector<string> split_data;
        boost::split( split_data, temp_val, boost::is_any_of( "," ) );

        if( split_data.size() != 2 )
            {
                throw std::runtime_error( "Values must a comma-delimited "
                                          "pair of doubles. "
                                        );
            }

        pair.first = std::strtod( split_data[ 0 ].c_str(), nullptr );
        pair.second = std::strtod( split_data[ 1 ].c_str(), nullptr );

        return in;
    }

    // used to convert pair( 3.4, 5.6 ) -> string '3.4,5.6'
    template<typename Stream>
    Stream& operator<<( Stream& out, const std::pair<double,double>& pair )
    {
        ostringstream stream;
        stream << pair.first << "," << pair.second;

        out << stream.str();
        return out;
    }

};

#endif // STREAM_TOOLS_HH_INCLUDED
