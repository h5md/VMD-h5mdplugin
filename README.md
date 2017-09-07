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

To compile the plugin, run `make`. This should work on
any Unix that has VMD and HDF5 installed. Under Ubuntu, the latter is
provided by the package *libhdf5-dev*.

Alternatively, you may build the plugin via CMake (`cmake && make`).

To load the plugin in VMD add the following line to your `~/.vmdrc` (replace the
directory to where you have compiled it):

        vmd_plugin_scandirectory /PATH_TO_SO_DIRECTORY/ h5mdplugin.so

Alternatively you may also load the plugin by linking h5mdplugin.so and libh5md.so into the vmd molfile plugins folder (e.g. under .../vmd/plugins/LINUXAMD64/molfile/). You can do this automatically via `make install`.

Remarks
-------

  - The plugin only works with three-dimensional position data.

  - The files `libh5md.c` and `libh5md.h` form a small library which
    provides common C routines to read and write H5MD files. This library is
    licensed under GPLv3.
  - In order to make full usage of VMD (for example to make VMDs output more "colourful") you need to specify the dataset /parameters/vmd_structure in your h5md file, the documentation for this dataset can be found in the file "Documentation VMD parameters".


Examples
--------
Example h5md files can be found in the folder "samples".
![Example of the representation of a h5md file in VMD](https://lists.gnu.org/archive/html/h5md-user/2013-08/pngf5euRoAsmj.png)

