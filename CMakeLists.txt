cmake_minimum_required(VERSION 3.9)
project(PepSIRF)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS program_options)

include_directories(include/modules/core)
include_directories(include/modules/parsers)
include_directories(include/modules/demux)
add_compile_options( "-Wall" "-std=c++11" "-Wextra" 
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable(pep_sirf src/main.cpp)

# Library common to all modules
add_library(pepsirf_common SHARED src/parsers/options_parser.cpp src/modules/core/options.cpp 
            src/modules/core/sequence.cpp 
			src/parsers/fastq_parser.cpp src/parsers/fasta_parser.cpp
			src/modules/core/options_factory.cpp)

add_library(modules_demux STATIC src/modules/demux/options_demux.cpp src/modules/demux/options_parser_demux.cpp src/parsers/options_parser.cpp 
            src/modules/core/options.cpp)

add_executable(pepsirf_test test/pepsirf_test.cpp)
target_include_directories(pepsirf_test PUBLIC lib)

target_link_libraries(pep_sirf PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux ${Boost_LIBRARIES})
target_link_libraries(pepsirf_test PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux ${Boost_LIBRARIES})
