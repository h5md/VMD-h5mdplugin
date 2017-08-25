############## test related ##############
# the following can be removed from the makefile and lib and plugin is still built (remove also "gtest_main.o libh5md_unittest" from the make target "all")

CXX=h5c++
CXXFLAGS=-Wall -O0 -fPIC -g

# Points to the root of fused googletest, relative to where this file is. make sure to adjust the 3 needed files (which are the only needed files) in this directory (correct #include "/gtest/gtest.h" to #include "gtest.h" since the header file is found in the same directory)
#change this directory
FUSED_GTEST_DIR = $(HOME)/src/gtest-1.7.0/fused-src/gtest

# Paths to the fused gtest files.
FUSED_GTEST_H = $(FUSED_GTEST_DIR)/gtest.h
FUSED_GTEST_ALL_CC = $(FUSED_GTEST_DIR)/gtest-all.cc

# Where to find the tests.
TEST_DIR = .

# Where to find gtest_main.cc.
GTEST_MAIN_CC = $(FUSED_GTEST_DIR)/gtest_main.cc

# Flags passed to the preprocessor.
# We have no idea here whether pthreads is available in the system, so
# disable its use.
CPPFLAGS += -I$(FUSED_GTEST_DIR) -DGTEST_HAS_PTHREAD=0 

PARENT_DIR = $(shell dirname `pwd`) 

all : gtest_main.o libh5md_unittest

# makes everything and runs the built sample test.
check : 
	-$(MAKE)
	-mkdir -p ./samples
	-./libh5md_unittest

gtest-all.o : $(FUSED_GTEST_H) $(FUSED_GTEST_ALL_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(FUSED_GTEST_ALL_CC)

gtest_main.o : $(FUSED_GTEST_H) $(GTEST_MAIN_CC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_MAIN_CC)

libh5md_test.o : $(TEST_DIR)/libh5md_test.cc $(TEST_DIR)/libh5md_test.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_DIR)/libh5md_test.cc

libh5md_unittest.o : $(TEST_DIR)/libh5md_unittest.cc \
                     $(TEST_DIR)/libh5md_test.hpp $(FUSED_GTEST_H)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_DIR)/libh5md_unittest.cc

libh5md_unittest : libh5md_unittest.o libh5md_test.o  gtest-all.o gtest_main.o
	- cd .. && $(MAKE)
	-$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wl,-rpath,$(PARENT_DIR) $^ -o $@ -L .. -lh5md
	
	
clean :
	$(RM) *~
	$(RM) libh5md_unittest
	$(RM) -rf samples
	$(RM) *.o
