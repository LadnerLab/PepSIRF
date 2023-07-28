#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iostream>
#include <string>

#include "options.h"

void error(std::string err_str);

void info(std::string info_str);

void warn(std::string warn_str);


#endif /* LOGGER_H */

