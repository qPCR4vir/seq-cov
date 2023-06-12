cmake_minimum_required (VERSION 3.12)
# set (CMAKE_BUILD_TYPE "Release")

# set(CMAKE_CXX_COMPILER C:/MinGW/bin/g++.exe -std=c++20 )
# set(CMAKE_CXX_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.31.31103/bin/Hostx64/x64/cl.exe" /std:c++latest )

project (seq-cov VERSION 1.10.02 LANGUAGES CXX)
 
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/nana" 
                 "${CMAKE_CURRENT_SOURCE_DIR}/../ExtLib/cmake-nana-build-${CMAKE_BUILD_TYPE}" )

# build app with seqan3
add_executable (seq-cov source/seq-cov.cpp)  
set_source_files_properties(source/seq-cov.cpp PROPERTIES COMPILE_OPTIONS "-Wa,-mbig-obj")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set_property(TARGET seq-cov PROPERTY CXX_STANDARD 20)
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