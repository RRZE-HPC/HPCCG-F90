The build system is configurable. Intermediate build products are all placed
in a top level directory named according to the configured tool chain. It is possible
to build multiple tool chain configurations in the same tree.

Dependencies are correctly configured only for the source files. If the build configuration is
changed a make clean is required before calling make.

======
Build instructions:
======

* Edit config.mk
* Review and adapt ./make/include_<TOOLCHAIN>.mk
* Call make


Detailed build instructions:

Options in config.mk and their meaning are:

TOOLCHAIN  -  The toolchain to use. According to this option the file
./make/include_<TOOLCHAIN>.mk will be included which contains all paths,
options, includes and library definitions. Possible values: IFORT or CRAY

ENABLE_MPI -  Enable MPI distributed memory parallelism. Possible values: YES or NO.

ENABLE_OMP - Enable usage of threads in OpenMP work sharing constructs. Possible values: YES or NO.

HAVE_COOP - Indicate if the compiler has support for collective intrinsics. Possible values: YES or NO.

TIMER_LEVEL - The granularity of timings in the code. This option will enable
different levels of timing instrumentation. The higher the number the more
instrumentation. Default is 1. At the moment timer routines up to level 2 are supported.

The following build variants are possible:

* CAF only: Set ENABLE_MPI and ENABLE_OMP to NO. HAVE_COOP must be set to YES.
* Hybrid CAF/OpenMP: Set ENABLE_OMP to YES
* Hybrid CAF/MPI parallelization: SET ENABLE_MPI to YES

Additional build targets:

make clean - Clean up all intermediate build products (object and module files)
make distclean - Clean up all build products
make - Build executable

======
Usage
======

The executable has three required options which are the dimension in x, y and z direction.

./HPCGG-IFORT  <nx>  <ny>  <nz>


