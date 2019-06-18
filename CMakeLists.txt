cmake_minimum_required(VERSION 3.9)
project(PepSIRF)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS program_options)

include_directories(include/modules/core)
include_directories(include/modules/parsers)
include_directories(include/modules/demux)
include_directories(lib/parallel_hashmap)
add_compile_options( "-Wall" "-std=c++11" "-Wextra" "-fopenmp" 
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable(pep_sirf src/main.cpp)


# Library common to all modules
add_library(pepsirf_common STATIC src/modules/parsers/options_parser.cpp src/modules/core/options.cpp 
            src/modules/core/sequence.cpp 
			src/modules/parsers/fastq_parser.cpp src/modules/parsers/fasta_parser.cpp
			src/modules/core/options_factory.cpp src/modules/core/options_parser_factory.cpp
			src/modules/core/module_factory.cpp src/modules/core/module_factory.cpp
			src/modules/core/module.cpp
			src/modules/core/sequence_indexer.cpp
		   )

add_library(modules_demux STATIC src/modules/demux/options_demux.cpp src/modules/demux/options_parser_demux.cpp src/modules/parsers/options_parser.cpp 
            src/modules/core/options.cpp src/modules/core/options_parser_factory.cpp
			src/modules/core/module_factory.cpp
			src/modules/core/module.cpp src/modules/demux/module_demux.cpp
			src/modules/parsers/fasta_parser.cpp src/modules/parsers/fastq_parser.cpp
            src/modules/core/sequence.cpp 
			src/modules/demux/samplelist_parser.cpp
			src/modules/core/sequence_indexer.cpp
		   )


option( BUILD_TESTS "Determines whether to build tests." ON )
if(BUILD_TESTS)
  add_executable(pepsirf_test test/pepsirf_test.cpp)
  include_directories(pepsirf_test lib)
  enable_testing()
  target_link_libraries(pepsirf_test PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux ${Boost_LIBRARIES})
  add_test( NAME test COMMAND pepsirf_test ) 
endif()

target_link_libraries(pep_sirf PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux ${Boost_LIBRARIES})
