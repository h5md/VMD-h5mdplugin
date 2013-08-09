# Where VMD is installed on your system
VMDDIR=/usr/local/lib/vmd

# Where the HDF5 library is installed
HDF5LDFLAGS=-L/usr/lib64
# Where the HDF5 include files are installed
HDF5INCLUDES=/usr/include

# Which libraries to use
HDF5LIBS=-lhdf5_cpp -lhdf5

# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# plugins/include of the VMD installation directory.
VMDINCLUDES=$(VMDDIR)/plugins/include 

CC=gcc
CFLAGS=-Wall -std=c99 -g -O0 -pedantic -fPIC
CXXFLAGS=$(CFLAGS) -g
CPPFLAGS=-I$(VMDINCLUDES) -I$(HDF5INCLUDES)
LDFLAGS=-L$(HDF5LDFLAGS)
LDLIBS=$(HDF5LIBS)
SHLD=$(CC)
SHLDFLAGS=-shared $(LDFLAGS)

all: h5mdtest h5mdplugin.so

h5mdtest: h5mdplugin.o

h5mdplugin.so: h5mdplugin.o
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) $(LDLIBS) $< -o h5mdplugin.so

clean:
	-rm h5mdplugin.o
	-rm h5mdplugin.so
	-rm h5mdtest
