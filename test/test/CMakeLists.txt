
cmake_minimum_required(VERSION 2.8.7)
if (NOT CMAKE_BUILD_TYPE)
                message(STATUS "No build type selected, default to Release")
                set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

project(LearnToCompress CXX C)
set(PROJECT_URL "https://github.com/jxxiaoshaoye/Learn-to-Compress")
set(PROJECT_DESCRIPTION "C++ library to compress with learned methods")

### Next line would be ideal, but it does not work on all platforms:
###set(CMAKE_CXX_STANDARD 11)
####
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")


# library target
include_directories(include)
add_library(LearnToCompress STATIC ./test.cpp
                                                                ${HEADERS}
                                                                )
enable_testing()


add_executable(${PROJECT_NAME} test.cpp )
target_link_libraries( LearnToCompress)
###add_test(  ${PROJECT_SOURCE_DIR}/sampledata.txt )
