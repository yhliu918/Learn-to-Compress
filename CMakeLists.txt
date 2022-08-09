cmake_minimum_required(VERSION 3.16)

project(example
    DESCRIPTION  "C++ library to compress with learned methods"
    HOMEPAGE_URL "https://github.com/yhliu918/Learn-to-Compress"
)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
find_package(Snappy)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
# add_compile_options(-O3 -g )#-fsanitize=address
# link_libraries(-fsanitize=address)

# .h document
set(HEADER_DIR ${PROJECT_SOURCE_DIR}/headers)
include_directories(${HEADER_DIR} /usr/include/eigen3)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/fsst)
# .c/cpp document
aux_source_directory(./experiments SRC_LIST)
add_library(lib ${SRC_LIST})
#  .so/.a
# make executable doc
add_executable(fix_int experiments/fix_int.cpp)
target_link_libraries(fix_int docs Snappy::snappy gmp gmpxx )
add_executable(auto_int experiments/auto_int.cpp)
target_link_libraries(auto_int docs Snappy::snappy gmp gmpxx )
add_executable(leco_int experiments/leco_int_template.cpp)
target_link_libraries(leco_int docs Snappy::snappy gmp gmpxx )
add_executable(leco_int_fix experiments/leco_int_fix_template.cpp)
target_link_libraries(leco_int_fix docs Snappy::snappy gmp gmpxx )
add_executable(leco_int_double experiments/leco_double_int_template.cpp)
target_link_libraries(leco_int_double docs Snappy::snappy gmp gmpxx )
add_executable(delta_int experiments/delta_int_template.cpp)
target_link_libraries(delta_int docs Snappy::snappy gmp gmpxx )
add_executable(leco_string experiments/example_fix_padding_template.cpp)
target_link_libraries(leco_string docs Snappy::snappy gmp gmpxx )
add_executable(leco_string_subset experiments/leco_string_subset.cpp)
target_link_libraries(leco_string_subset docs Snappy::snappy gmp gmpxx )
add_executable(leco_string_subset_shift experiments/leco_string_subset_shift.cpp)
target_link_libraries(leco_string_subset_shift docs Snappy::snappy gmp gmpxx )
add_executable(test_bitmap experiments/test_bitmap.cpp)
target_link_libraries(test_bitmap docs Snappy::snappy gmp gmpxx )
add_executable(fsst_string experiments/example_fsst.cpp)
target_link_libraries(fsst_string docs fsst gmp gmpxx )


add_executable(leco_dp experiments/leco_dp.cpp)
target_link_libraries(leco_dp docs Snappy::snappy gmp gmpxx )

include(CTest)
enable_testing()
add_test(NAME ${PROJECT_NAME}
  COMMAND ${PROJECT_NAME}
)
