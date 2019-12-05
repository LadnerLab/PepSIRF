cmake_minimum_required(VERSION 3.9)
project(PepSIRF LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)

find_package(Boost REQUIRED COMPONENTS program_options regex filesystem)

include_directories(include/modules/core)
include_directories(include/modules/normalize)
include_directories(include/modules/parsers)
include_directories(include/modules/demux)
include_directories(include/modules/deconv)
include_directories(include/modules/subjoin)
include_directories(lib/parallel_hashmap)
add_compile_options( "-Wall" "-std=c++11" "-Wextra" "-fvisibility=hidden"
                     "-fvisibility-inlines-hidden"
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>"
                   )


set(PepSIRF_MODULES pepsirf_common modules_demux modules_deconv modules_normalize modules_subjoin )                

if(OpenMP_FOUND)
  message( "OpenMP enabled" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  add_definitions( -DENABLE_OPENMP )

else()
  message( "WARNING: OpenMP not found, parallelism disabled." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-value")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable(pep_sirf src/main.cpp)

# Library common to all modules
add_library(pepsirf_common STATIC src/modules/parsers/options_parser.cpp src/modules/core/options.cpp 
            src/modules/core/sequence.cpp 
			src/modules/parsers/fastq_parser.cpp src/modules/parsers/fasta_parser.cpp
			src/modules/core/module.cpp
			src/modules/core/sequence_indexer.cpp
			src/modules/core/fastq_score.cpp
			src/modules/core/fs_tools.cpp
            src/modules/core/peptide.cpp
            src/modules/core/module_initializer.cpp
            src/modules/core/cli_validator.cpp
            src/modules/core/modules.cpp
           )

add_library(modules_demux STATIC src/modules/demux/options_demux.cpp src/modules/demux/options_parser_demux.cpp src/modules/parsers/options_parser.cpp 
            src/modules/core/options.cpp 
			src/modules/core/module.cpp src/modules/demux/module_demux.cpp
			src/modules/parsers/fasta_parser.cpp src/modules/parsers/fastq_parser.cpp
            src/modules/core/sequence.cpp 
			src/modules/demux/samplelist_parser.cpp
			src/modules/core/sequence_indexer.cpp
			src/modules/core/fastq_score.cpp
		   )

add_library( modules_deconv STATIC src/modules/deconv/options_deconv.cpp src/modules/deconv/options_parser_deconv.cpp 
				src/modules/core/module.cpp src/modules/deconv/module_deconv.cpp
				src/modules/core/fs_tools.cpp
                src/modules/deconv/scored_peptide.cpp
                src/modules/core/peptide.cpp
                src/modules/deconv/species_data.cpp
		)

add_library( modules_normalize STATIC src/modules/normalize/options_normalize.cpp
             src/modules/normalize/options_parser_normalize.cpp
			 src/modules/core/module.cpp src/modules/normalize/module_normalize.cpp
		   )

add_library( modules_subjoin STATIC src/modules/subjoin/options_subjoin.cpp
             src/modules/subjoin/options_parser_subjoin.cpp
             src/modules/core/module.cpp src/modules/subjoin/module_subjoin.cpp
           )

           



option( BUILD_TESTS "Determines whether to build tests." ON )
if(BUILD_TESTS)
  add_executable(pepsirf_test test/pepsirf_test.cpp)
  add_definitions( -DMATRIX_CHECK_BOUNDS )
  include_directories(pepsirf_test lib)
  enable_testing()
if(OpenMP_FOUND)
  target_link_libraries(pepsirf_test PUBLIC OpenMP::OpenMP_CXX ${PepSIRF_MODULES} ${Boost_LIBRARIES} )

else()
  target_link_libraries(pepsirf_test PUBLIC ${PepSIRF_MODULES} ${Boost_LIBRARIES} )
endif()

  add_test( NAME test COMMAND pepsirf_test ) 
endif()

if(OpenMP_FOUND)
  target_link_libraries(pep_sirf PUBLIC OpenMP::OpenMP_CXX ${PepSIRF_MODULES} ${Boost_LIBRARIES} )
else()
  target_link_libraries(pep_sirf PUBLIC ${PepSIRF_MODULES} ${Boost_LIBRARIES} )
endif()
