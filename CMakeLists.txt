cmake_minimum_required(VERSION 3.16)

project(example
    DESCRIPTION  "C++ library to compress with learned methods"
)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 17)
find_package(Snappy)
set(CMAKE_CXX_FLAGS_NDEBUG "-g -O0  -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
# add_compile_options(-fsanitize=address)
# link_libraries(-fsanitize=address)

# .h document
set(HEADER_DIR ${PROJECT_SOURCE_DIR}/headers)
include_directories(${HEADER_DIR} /usr/include/eigen3)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/fsst)
# .c/cpp document
aux_source_directory(./experiments SRC_LIST)
add_library(leco ${SRC_LIST})
target_link_libraries(leco Snappy::snappy gmp gmpxx fsst)
#  .so/.a
# make executable doc

file(GLOB_RECURSE CPP_FILES ${PROJECT_SOURCE_DIR}/experiments/*.cpp )
foreach(A_CPP_FILE ${CPP_FILES})
    get_filename_component(EXECUTABLE_NAME ${A_CPP_FILE} NAME_WE)
    add_executable(${EXECUTABLE_NAME} ${A_CPP_FILE})
    target_link_libraries(${EXECUTABLE_NAME} docs Snappy::snappy gmp gmpxx fsst)
endforeach(A_CPP_FILE ${CPP_FILES})


# add_executable(RNS_int experiments/RNS_int.cpp)
# target_link_libraries(RNS_int docs Snappy::snappy gmp gmpxx )
# add_executable(huffman_int experiments/huffman_int.cpp)
# target_link_libraries(huffman_int docs Snappy::snappy gmp gmpxx )
# add_executable(snappy_int experiments/snappy_int.cpp)
# target_link_libraries(snappy_int docs Snappy::snappy gmp gmpxx )
# add_executable(fix_int experiments/fix_int.cpp)
# target_link_libraries(fix_int docs Snappy::snappy gmp gmpxx )
# add_executable(auto_int experiments/auto_int.cpp)
# target_link_libraries(auto_int docs Snappy::snappy gmp gmpxx )
# add_executable(leco_int experiments/leco_int_template.cpp)
# target_link_libraries(leco_int docs Snappy::snappy gmp gmpxx )
# add_executable(poly_mix experiments/pol_mix_template.cpp)
# target_link_libraries(poly_mix docs Snappy::snappy gmp gmpxx )
# add_executable(poly_mix_metric experiments/pol_mix_metric_template.cpp)
# target_link_libraries(poly_mix_metric docs Snappy::snappy gmp gmpxx )
# add_executable(poly_int experiments/poly_int.cpp)
# target_link_libraries(poly_int docs Snappy::snappy gmp gmpxx )

# add_executable(leco_int_fix experiments/leco_int_fix_template.cpp)
# target_link_libraries(leco_int_fix docs Snappy::snappy gmp gmpxx )
# add_executable(poly_int_fix experiments/leco_int_pol_fix_template.cpp)
# target_link_libraries(poly_int_fix docs Snappy::snappy gmp gmpxx )
# add_executable(Delta_int_fix experiments/Delta_int_fix_template.cpp)
# target_link_libraries(Delta_int_fix docs Snappy::snappy gmp gmpxx )
# add_executable(FOR_int_fix experiments/FOR_int_fix_template.cpp)
# target_link_libraries(FOR_int_fix docs Snappy::snappy gmp gmpxx )
# add_executable(leco_int_double experiments/leco_double_int_template.cpp)
# target_link_libraries(leco_int_double docs Snappy::snappy gmp gmpxx )
# add_executable(leco_int_double_wo_round experiments/leco_double_int_template_wo_round.cpp)
# target_link_libraries(leco_int_double_wo_round docs Snappy::snappy gmp gmpxx )
# add_executable(leco_lp_cost experiments/leco_lp_cost.cpp)
# target_link_libraries(leco_lp_cost docs Snappy::snappy gmp gmpxx )

# add_executable(delta_int experiments/delta_int_template.cpp)
# target_link_libraries(delta_int docs Snappy::snappy gmp gmpxx )
# add_executable(leco_string_subset experiments/leco_string_subset.cpp)
# target_link_libraries(leco_string_subset docs Snappy::snappy gmp gmpxx )
# add_executable(leco_string_subset_shift experiments/leco_string_subset_shift.cpp)
# target_link_libraries(leco_string_subset_shift docs Snappy::snappy gmp gmpxx )
# add_executable(test_bitmap experiments/test_bitmap.cpp)
# target_link_libraries(test_bitmap docs Snappy::snappy gmp gmpxx )
# add_executable(fsst_string experiments/example_fsst.cpp)
# target_link_libraries(fsst_string docs fsst gmp gmpxx )


# add_executable(leco_dp experiments/leco_dp.cpp)
# target_link_libraries(leco_dp docs Snappy::snappy gmp gmpxx )

# add_executable(bit_pack experiments/Bit_packing_template.cpp)
# target_link_libraries(bit_pack docs Snappy::snappy gmp gmpxx )

# add_executable(leco_dict experiments/leco_test_dictionary.cpp)
# target_link_libraries(leco_dict docs Snappy::snappy gmp gmpxx )
# add_executable(for_dict experiments/for_test_dictionary.cpp)
# target_link_libraries(for_dict docs Snappy::snappy gmp gmpxx )
# add_executable(raw_dict experiments/raw_test_dictionary.cpp)
# target_link_libraries(raw_dict docs Snappy::snappy gmp gmpxx )

include(CTest)
enable_testing()
add_test(NAME ${PROJECT_NAME}
  COMMAND ${PROJECT_NAME}
)
