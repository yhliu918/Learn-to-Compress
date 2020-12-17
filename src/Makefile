.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h
# replace the YOURCXX variable with a path to a C++11 compatible compiler.
ifeq ($(INTEL), 1)
# if you wish to use the Intel compiler, please do "make INTEL=1".
    YOURCXX ?= /opt/intel/bin/icpc
ifeq ($(DEBUG),1)
    CXXFLAGS =  -std=c++11 -O3 -Wall -ansi  -DDEBUG=1 -D_GLIBCXX_DEBUG   -ggdb
else
    CXXFLAGS =  -std=c++11 -O3 -Wall -ansi  -DNDEBUG=1  -ggdb
endif # debug
else #intel
    YOURCXX ?= g++
ifeq ($(DEBUG),1)
    CXXFLAGS =   -std=c++11   -ggdb -DDEBUG=1 -D_GLIBCXX_DEBUG -Wall -Wextra  -Wcast-align  
else
    CXXFLAGS =  -std=c++11   -O3 -Wall -Wextra  -Wcast-align  
endif #debug
endif #intel


 
all:  testfix
	@echo "type ./test " 
	@echo "if you have OpenMP support, type make testmp" 


bpacking.o: ../headers/bpacking.h bpacking.cpp
	$(CXX) $(CXXFLAGS) -c bpacking.cpp -Iinclude
    
varintencode.o: ../headers/varintencode.h varintencode.cpp
	$(CCX) $(CXXFLAGS) -c varintencode.cpp -Iinclude  

varintdecode.o: ../headers/varintdecode.h varintencode.cpp
	$(CCX) $(CXXFLAGS) -c varintdecode.cpp -Iinclude  

Utils.o: ../headers/Utils.h Utils.cpp
	$(CCX) $(CXXFLAGS) -c Utils.cpp -Iinclude  

LinearRegression.o: ../headers/LinearRegression.h LinearRegression.cpp
	$(CCX) $(CXXFLAGS) -c LinearRegression.cpp -Iinclude


obj = test.cpp ../headers/common.h ../headers/codecfactory.h \
 ../headers/common.h ../headers/codecs.h ../headers/util.h \
 ../headers/FOR.h  bpacking.o ../headers/forutil.h \
 ../headers/piecewise.h ../headers/time.h ../headers/bit_opt.h \
 ../headers/maskvbyte.h varintdecode.o varintencode.o

fixobj = testpiecewise_fix.cpp ../headers/common.h \
 ../headers/codecfactory.h ../headers/common.h ../headers/codecs.h \
 ../headers/util.h ../headers/FOR.h bpacking.o \
 ../headers/forutil.h ../headers/piecewise.h ../headers/time.h \
 ../headers/bit_opt.h ../headers/piecewise_fix.h \
 LinearRegression.o ../headers/maskvbyte.h \
 varintdecode.o varintencode.o ../headers/time.h Utils.o
 


test: $(obj)
	$(CXX)  $(CXXFLAGS)  -o test $(obj)  -Iinclude
    
testfix:  $(fixobj)
	$(CXX)  $(CXXFLAGS)  -o testfix $(fixobj)  -Iinclude

clean: 
	rm -f *.o testfix





