cmake_minimum_required(VERSION 3.9)
project(PepSIRF)

find_package(OpenMP REQUIRED)
include_directories(include)
add_compile_options( "-Wall" "-std=c++11" "-Wextra" 
                     "$<$<CONFIG:DEBUG>:-O0;-g>" 
					 "$<$<CONFIG:RELEASE>:-O3>")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(pep_sirf src/main.cpp)
add_library(pepsirf_lib SHARED src/Student.cpp)
add_executable(pepsirf_test test/pepsirf_test.cpp)
target_link_libraries(pep_sirf PRIVATE OpenMP::OpenMP_CXX pepsirf_lib)
target_link_libraries(pepsirf_test PRIVATE OpenMP::OpenMP_CXX pepsirf_lib)
