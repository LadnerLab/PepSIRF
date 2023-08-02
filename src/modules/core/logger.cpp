#include "logger.h"

std::ofstream Log::logstream;

void Log::close()
{
    logstream.close();
}

void Log::error(const std::string err_str)
{
    std::cout << "Error: " << err_str;
    logstream << "Error: " << err_str;
}

void Log::info(const std::string info_str)
{
    std::cout << info_str;
    logstream << info_str;
}

void Log::open(const std::string logfile)
{
    logstream.open(logfile, std::ofstream::out);
}

void Log::warn(const std::string warn_str)
{
    std::cout << "Warning: " << warn_str;
    logstream << "Warning: " << warn_str;
}

