cmake_minimum_required(VERSION 3.9)
project(PepSIRF LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)
find_package(ZLIB)

find_package(Boost REQUIRED COMPONENTS program_options regex filesystem iostreams)

include_directories(include/modules/core)
include_directories(include/modules/normalize)
include_directories(include/modules/parsers)
include_directories(include/modules/demux)
include_directories(include/modules/deconv)
include_directories(include/modules/subjoin)
include_directories(include/modules/zscore)
include_directories(include/modules/bin)
include_directories(include/modules/link)
include_directories(include/modules/s_enrich)
include_directories(include/modules/p_enrich)
include_directories(include/modules/info)
include_directories(lib)

add_compile_options( "-Wall" "-std=c++11" "-Wextra" "-fvisibility=hidden"
                     "-fvisibility-inlines-hidden"
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>"
                   )


set(PepSIRF_MODULES
    pepsirf_common
    modules_demux modules_deconv 
    modules_normalize modules_subjoin 
	modules_zscore modules_bin 
	modules_link
	modules_s_enrich
	modules_p_enrich
	modules_info
   )                

list( APPEND PepSIRF_LINK_LIBS
  ${PepSIRF_MODULES}
  ${Boost_LIBRARIES}
)

if(OpenMP_FOUND)
  message( "OpenMP enabled" )
  list( APPEND PepSIRF_LINK_LIBS OpenMP::OpenMP_CXX )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  add_definitions( -DENABLE_OPENMP )

else()
  message( "WARNING: OpenMP not found, parallelism disabled." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-value")
endif()

if(ZLIB_FOUND)
  message( "Gzip file manipulation enabled." )
  add_definitions( -DZLIB_ENABLED )
  list( APPEND PepSIRF_LINK_LIBS ZLIB::ZLIB )
else()
  message( "WARNING: ZLIB was not found, demux cannot process gzipped files." )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${Boost_INCLUDE_DIRS})

add_executable(pepsirf src/main.cpp)

# Library common to all modules
add_library(pepsirf_common STATIC src/modules/parsers/options_parser.cpp src/modules/core/options.cpp 
            src/modules/core/sequence.cpp 
            src/modules/core/file_io.cpp 
			src/modules/parsers/fastq_parser.cpp src/modules/parsers/fasta_parser.cpp
			src/modules/core/module.cpp
			src/modules/core/sequence_indexer.cpp
			src/modules/core/fastq_score.cpp
			src/modules/core/fs_tools.cpp
            src/modules/core/peptide.cpp
            src/modules/core/module_initializer.cpp
            src/modules/core/cli_validator.cpp
            src/modules/core/modules.cpp
			src/modules/core/peptide_scoring.cpp
			src/modules/core/predicate.cpp
           )

add_library(modules_demux STATIC src/modules/demux/options_demux.cpp src/modules/demux/options_parser_demux.cpp src/modules/parsers/options_parser.cpp 
            src/modules/core/options.cpp 
            src/modules/core/file_io.cpp 
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
			 src/modules/core/peptide_scoring.cpp
           )

add_library( modules_subjoin STATIC src/modules/subjoin/options_subjoin.cpp
             src/modules/subjoin/options_parser_subjoin.cpp
             src/modules/core/module.cpp src/modules/subjoin/module_subjoin.cpp
           )

add_library( modules_zscore STATIC src/modules/zscore/module_zscore.cpp
			 src/modules/zscore/options_zscore.cpp
			 src/modules/zscore/options_parser_zscore.cpp
			 src/modules/bin/peptide_bin.cpp
		   )
add_library( modules_bin STATIC src/modules/bin/module_bin.cpp 
			 src/modules/bin/options_bin.cpp
			 src/modules/bin/options_parser_bin.cpp
			 src/modules/bin/peptide_bin.cpp
			 src/modules/bin/probe_rank.cpp
           )

add_library( modules_link STATIC src/modules/link/module_link.cpp
				src/modules/link/options_link.cpp
				src/modules/link/options_parser_link.cpp
		   )
add_library( modules_s_enrich STATIC src/modules/s_enrich/module_s_enrich.cpp
				src/modules/s_enrich/options_s_enrich.cpp
				src/modules/s_enrich/options_parser_s_enrich.cpp
		   )

add_library( modules_p_enrich STATIC src/modules/p_enrich/module_p_enrich.cpp
				src/modules/p_enrich/options_p_enrich.cpp
				src/modules/p_enrich/options_parser_p_enrich.cpp
				src/modules/core/predicate.cpp
		   )

add_library( modules_info STATIC src/modules/info/module_info.cpp
			 src/modules/info/options_info.cpp
			 src/modules/info/options_parser_info.cpp
		   )


option( BUILD_TESTS "Determines whether to build tests." ON )

if(BUILD_TESTS)

  add_executable(pepsirf_test test/pepsirf_test.cpp)
  add_definitions( -DMATRIX_CHECK_BOUNDS )
  include_directories(pepsirf_test lib)
  enable_testing()
  target_link_libraries(pepsirf_test PUBLIC ${PepSIRF_LINK_LIBS} )
  add_test( NAME test COMMAND pepsirf_test ) 

endif()

target_link_libraries(pepsirf PUBLIC ${PepSIRF_LINK_LIBS} )
