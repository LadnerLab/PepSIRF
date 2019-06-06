#include "options_demux.h"


std::string options_demux::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--input_r1          " << input_r1_fname << "\n" <<
                  " --input_r2         " << input_r2_fname << "\n" << 
                  " --output           " << output_fname << "\n" <<
                  " --index            " << index_fname << "\n" << 
                  " --samplelist       " << samplelist_fname << "\n"
                  " --library          " << library_fname << "\n" <<
                  " --read_per_loop    " << read_per_loop << "\n" <<
                  " --num_threads      " << num_threads << "\n" <<
                  "\n";

    return str_stream.str();
}
