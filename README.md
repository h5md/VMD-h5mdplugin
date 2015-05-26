Rationale
---------

This plugin enables [VMD](http://www.ks.uiuc.edu/Research/vmd>) to display
dynamic molecular data stored in [H5MD](http://nongnu.org/h5md) files.

**H5MD is a structured, binary file format** for molecular data, such as
simulation trajectories, molecular structures, or particle subsets with
associated data. In addition, there is scope for custom data.

**H5MD is built on the [HDF5 library](http://www.hdfgroup.org/HDF5)** and
provides, among others, multi-dimensional datasets, data compression,
high-performance parallel I/O, and portability.

**H5MD files are created** by these simulation packages:

  - [HAL’s MD package](http://halmd.org)
  - [ESPResSo](http://espressomd.org)
  - [LAMMPS](http://lammps.sandia.gov)
    ([h5md dump style](https://github.com/pdebuyl/lammps/tree/start_dump_h5md/src/USER-H5MD))

Find out about more [H5MD-related software](http://nongnu.org/h5md/software.html).

**The H5MD file specification** is hosted at http://nongnu.org/h5md and has
been published originally in
P. de Buyl, P. H. Colberg, and F. Höfling,
[Comput. Phys. Commun. **185**, 1546 (2014)](http://dx.doi.org/10.1016/j.cpc.2014.01.018>),
see also [arXiv:1308.6382](http://arxiv.org/abs/1308.6382).


Installation
------------

To compile the plugin, adapt the Makefile and run `make`. This should work on
any Unix that has VMD and HDF5 installed. For example under Ubuntu you need the
package *libhdf5-dev* and the VMD sources.

Alternatively, you may build the plugin with CMake, which requires the
environment variable `VMDDIR` to be set.

To load the plugin in VMD add the following line to your `~/.vmdrc` (replace the
directory to where you have compiled it):

        vmd_plugin_scandirectory /PATH_TO_SO_DIRECTORY/ h5mdplugin.so


Remarks
-------

  - The plugin only works with three-dimensional position data.

  - The files `libh5md.c` and `libh5md.h` form a small library which
    provides common C routines to read and write H5MD files. This library is
    licensed under GPLv3.

