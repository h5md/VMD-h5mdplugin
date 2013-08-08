# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# plugins/include of the VMD installation directory.

# ICP:
VMDDIR=/usr/local/lib/vmd

VMDINCLUDES=$(VMDDIR)/plugins/include 

# comment this line, if zlib is not available
_USE_ZLIB=1
#DEBUG=1

CC=gcc
CFLAGS=-Wall -std=c99 -g -O0 -pedantic -fPIC
CXXFLAGS=$(CFLAGS) -g
CPPFLAGS=-I$(VMDINCLUDES) 
LDFLAGS=-L/usr/lib64  
LDLIBS=-ltcl8.5 -lhdf5_cpp -lhdf5
SHLD=$(CC)
SHLDFLAGS=-shared $(LDFLAGS)

ifdef _USE_ZLIB
# if you want to enable compressed files, use these
CFLAGS += -D_USE_ZLIB
LDLIBS += -lz
endif

ifdef DEBUG
CFLAGS += -DDEBUG
endif

all: h5mdtest h5mdplugin.so

h5mdtest: h5mdplugin.o

h5mdplugin.so: h5mdplugin.o
	$(SHLD) $(CFLAGS) $(SHLDFLAGS) $(LDLIBS) $< -o h5mdplugin.so

clean:
	-rm h5mdplugin.o
	-rm h5mdplugin.so
	-rm h5mdtest
