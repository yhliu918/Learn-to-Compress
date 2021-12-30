
cmake_minimum_required(VERSION 2.8.7)
if (NOT CMAKE_BUILD_TYPE)
                message(STATUS "No build type selected, default to Release")
                set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

project(example)
set(PROJECT_URL "https://github.com/jxxiaoshaoye/Learn-to-Compress")
set(PROJECT_DESCRIPTION "C++ library to compress with learned methods")
find_package(Snappy)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11 -O2 -g -Wall -ansi -ggdb1 -Iinclude -I /usr/include/Eigen")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    add_definitions("-std=c++11")

    set(CMAKE_CXX_FLAGS_MAIN "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC -march=native -O0 -Wall -Wextra -pedantic -Wshadow -ansi -mfma -DNDEBUG=1  -g")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


set(HEADER_DIR ${PROJECT_SOURCE_DIR}/headers)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# .c/cpp document
aux_source_directory(./experiments SRC_LIST)
add_library(lib ${SRC_LIST})
##set_target_properties(lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
#  .so/.a

# .h document
include_directories(${HEADER_DIR})

# make executable doc
add_executable(${PROJECT_NAME}  experiments/example.cpp)
target_link_libraries( ${PROJECT_NAME} docs Snappy::snappy gmp gmpxx)

include(CTest)
enable_testing()
add_test(NAME ${PROJECT_NAME}
  COMMAND ${PROJECT_NAME}
  )
