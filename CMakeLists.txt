cmake_minimum_required(VERSION 3.9)
project(PepSIRF)

find_package(OpenMP REQUIRED)
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_compile_options( "-Wall" "-std=c++11" "-Wextra" 
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable( pep_sirf ${SOURCES})
target_link_libraries(pep_sirf PRIVATE OpenMP::OpenMP_CXX)
