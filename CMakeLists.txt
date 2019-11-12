cmake_minimum_required(VERSION 3.9)
project(PepSIRF LANGUAGES CXX)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)

if(OpenMP_FOUND)
  message( "OpenMP enabled" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  add_definitions( -DENABLE_OPENMP )

else()
  message( "WARNING: OpenMP not found, parallelism disabled." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-value")
endif()

find_package(Boost REQUIRED COMPONENTS program_options regex filesystem)

include_directories(include/modules/core)
include_directories(include/modules/parsers)
include_directories(include/modules/demux)
include_directories(include/modules/deconv)
include_directories(lib/parallel_hashmap)
add_compile_options( "-Wall" "-std=c++11" "-Wextra"
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
			src/modules/core/fastq_score.cpp
			src/modules/core/fs_tools.cpp
            src/modules/core/peptide.cpp
           )

add_library(modules_demux STATIC src/modules/demux/options_demux.cpp src/modules/demux/options_parser_demux.cpp src/modules/parsers/options_parser.cpp 
            src/modules/core/options.cpp src/modules/core/options_parser_factory.cpp
			src/modules/core/module_factory.cpp
			src/modules/core/module.cpp src/modules/demux/module_demux.cpp
			src/modules/parsers/fasta_parser.cpp src/modules/parsers/fastq_parser.cpp
            src/modules/core/sequence.cpp 
			src/modules/demux/samplelist_parser.cpp
			src/modules/core/sequence_indexer.cpp
			src/modules/core/fastq_score.cpp
		   )

add_library( modules_deconv STATIC src/modules/deconv/options_deconv.cpp src/modules/deconv/options_parser_deconv.cpp 
				src/modules/core/options.cpp src/modules/core/options_parser_factory.cpp
				src/modules/core/module_factory.cpp
				src/modules/core/module.cpp src/modules/deconv/module_deconv.cpp
				src/modules/core/fs_tools.cpp
                src/modules/deconv/scored_peptide.cpp
                src/modules/core/peptide.cpp
                src/modules/deconv/species_data.cpp
		)


option( BUILD_TESTS "Determines whether to build tests." ON )
if(BUILD_TESTS)
  add_executable(pepsirf_test test/pepsirf_test.cpp)
  include_directories(pepsirf_test lib)
  enable_testing()
  target_link_libraries(pepsirf_test PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux modules_deconv ${Boost_LIBRARIES} )
  add_test( NAME test COMMAND pepsirf_test ) 
endif()

target_link_libraries(pep_sirf PUBLIC OpenMP::OpenMP_CXX pepsirf_common modules_demux modules_deconv ${Boost_LIBRARIES} )
