# Where VMD is installed on your system
VMDDIR=/usr/local/lib/vmd

# Where the HDF5 library is installed
HDF5LDFLAGS=-L/usr/lib64/mpi/gcc/openmpi/lib64
HDF5CPPFLAGS=-I/usr/include

# Which libraries to use
HDF5LIBS=-lhdf5_hl -lhdf5
HDF5LIBSCplusplus=-lhdf5_cpp -lhdf5 

# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# plugins/include of the VMD installation directory.
VMDCPPFLAGS=-I$(VMDDIR)/plugins/include 

#CC=clang #for better warnings
CC=gcc
#compiler switches http://gcc.gnu.org/onlinedocs/gcc-3.4.3/gnat_ugn_unw/Switches-for-gcc.html
CFLAGS=-Wall -Wuninitialized -std=c99 -O0 -pedantic -fPIC -g
CPPFLAGS=$(VMDCPPFLAGS) $(HDF5CPPFLAGS)
LDFLAGS=$(HDF5LDFLAGS) -L. 
LDLIBS=$(HDF5LIBS)
SHLD=$(CC)
SHLDFLAGS=-shared -Wl,--no-undefined $(LDFLAGS)

all: h5mdtest libh5md.so h5mdplugin.so

#nasty
h5mdtest: h5mdtest.c h5mdplugin.so libh5md.so
	$(CC) h5mdtest.c -o h5mdtest h5mdplugin.c $(CPPFLAGS) -std=c99 $(LDFLAGS) -Wl,-rpath,$(shell pwd) $(HDF5LIBS) -lh5md


libh5md.so: libh5md.o
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) $< -o libh5md.so $(HDF5LIBS) -lm

h5mdplugin.so: h5mdplugin.o libh5md.so
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) -Wl,-rpath,$(shell pwd) $< -o h5mdplugin.so $(LDLIBS) -lh5md

check:
	cd tests && $(MAKE) check

#remove all files generated by make
clean:
	-rm h5mdplugin.so
	-rm h5mdtest
	-rm libh5md.so
	-rm *.o
	-rm *~
	-cd tests && $(MAKE) clean
	

	
