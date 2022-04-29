cmake_minimum_required (VERSION 3.12)
# set (CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_COMPILER C:/Prog/Cpp/MinGW/bin/g++.exe -std=c++20 )

project (seq-cov VERSION 1.00.00 LANGUAGES CXX)
 
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/nana" "${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/cmake-nana-build-${CMAKE_BUILD_TYPE}" )

 
# build app with seqan3
add_executable (seq-cov source/seq-cov.cpp)
#set_property(TARGET seq-cov PROPERTY CXX_STANDARD C++20)
target_link_libraries (seq-cov PUBLIC -static)
target_link_libraries(seq-cov
        PUBLIC
            $<$<CXX_COMPILER_ID:GNU>:-static-libgcc>
            $<$<CXX_COMPILER_ID:GNU>:-static-libstdc++>
            
            $<$<CXX_COMPILER_ID:Clang>:-static-libgcc>
            $<$<CXX_COMPILER_ID:Clang>:-static-libstdc++>
            )
target_link_libraries (seq-cov PUBLIC seqan3::seqan3 nana)
install(TARGETS seq-cov DESTINATION "bin")