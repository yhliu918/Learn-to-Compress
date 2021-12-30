add_definitions("-msse3 -msse4.1 -mssse3 -msse2")
add_library(docs SHARED bpacking.cpp varintencode.cpp varintdecode.cpp Utils.cpp LinearRegression.cpp RANSAC.cpp bignumber.cpp  )#MLP.cpp
##set_target_properties(docs PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
##add_library(docs SHARED ${DOC_LIST})
#  .so/.a
## link_directories(.)




