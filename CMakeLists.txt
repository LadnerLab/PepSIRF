cmake_minimum_required(VERSION 3.9)
project(PepSIRF)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS program_options)

include_directories(include)
add_compile_options( "-Wall" "-std=c++11" "-Wextra" 
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable(pep_sirf src/main.cpp)
add_library(pepsirf_lib STATIC src/parsers/options_parser.cpp src/options.cpp src/sequence.cpp src/parsers/fastq_parser.cpp src/parsers/fasta_parser.cpp)
add_executable(pepsirf_test test/pepsirf_test.cpp)
target_include_directories(pepsirf_test PUBLIC lib)

target_link_libraries(pep_sirf PUBLIC OpenMP::OpenMP_CXX pepsirf_lib ${Boost_LIBRARIES})
target_link_libraries(pepsirf_test PUBLIC OpenMP::OpenMP_CXX pepsirf_lib ${Boost_LIBRARIES})
