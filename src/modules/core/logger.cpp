#include "logger.h"

std::ofstream* Log::logstream;

void Log::error(std::string err_str)
{
    std::cout << "Error: " << err_str;
    *logstream << "Error: " << err_str;
}

void Log::info(std::string info_str)
{
    std::cout << info_str;
    *logstream << info_str;
}

void Log::warn(std::string warn_str)
{
    std::cout << "Warning: " << warn_str;
    *logstream << "Warning: " << warn_str;
}

