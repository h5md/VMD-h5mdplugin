

# Where VMD is installed on your system
VMDDIR=/usr/local/lib/vmd

# Where the HDF5 library is installed
HDF5LDFLAGS=-L/usr/lib64
HDF5CPPFLAGS=-I/usr/include

# Which libraries to use
HDF5LIBS=-lhdf5_hl -lhdf5
HDF5LIBSCplusplus=-lhdf5_cpp -lhdf5 

# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# plugins/include of the VMD installation directory.
VMDCPPFLAGS=-I$(VMDDIR)/plugins/include 

CC=gcc
CFLAGS=-Wall -std=c99 -g -O0 -pedantic -fPIC
CXXFLAGS=$(CFLAGS) -g
CPPFLAGS=$(VMDCPPFLAGS) $(HDF5CPPFLAGS)
LDFLAGS=$(HDF5LDFLAGS) -L. 
LDLIBS=$(HDF5LIBS)
SHLD=$(CC)
SHLDFLAGS=-shared -Wl,--no-undefined $(LDFLAGS)

all: h5mdtest libh5md.so h5mdplugin.so test_libh5md

#nasty
h5mdtest: h5mdtest.c h5mdplugin.so libh5md.so
	$(CC) h5mdtest.c -o h5mdtest h5mdplugin.c $(VMDCPPFLAGS) -std=c99 $(LDFLAGS) $(HDF5LIBS) -lh5md -Wl,-rpath,$(shell pwd)


libh5md.so: libh5md.o
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) $(HDF5LIBS) $< -o libh5md.so

h5mdplugin.so: h5mdplugin.o libh5md.so
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) $(LDLIBS) -lh5md -Wl,-rpath,$(shell pwd) $< -o h5mdplugin.so

#nasty
test_libh5md: test_libh5md.c libh5md.so
	$(CC) test_libh5md.c -o test_libh5md $(VMDCPPFLAGS) -std=c99 $(LDFLAGS) $(HDF5LIBS) -lh5md -Wl,-rpath,$(shell pwd)

clean:
	-rm h5mdplugin.o
	-rm h5mdplugin.so
	-rm h5mdtest
	-rm libh5md.o
	-rm libh5md.so
	-rm test_libh5md
	-rm *~
