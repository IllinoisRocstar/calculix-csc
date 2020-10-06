CalculiX -- IR Extras
---------------------
This is the vanilla CalculiX project (forked from v2.15) with a few improvments and additions as follows:
* IMPACT-based multiphysics coupling module (CSC)
* CMake-based build system
* Homogenization module
* ExodusII output and GPU solver support (GPU is not completely tested, using patches released by Peter A. Gustafson on 8/14/2019).

## Version
Version 2.15.2

The project follows the versioning convention of the original CalculiX. The versioning will be major.minor.patch and increments are as following:
* major: shows major version of the CalculiX project 
* minor: shows minor version of the CalculiX project
* patch: will be incremented as we release more functionality to the project

## Options
The following table contains all available CMake options to configure the project. The necessary third-party libraries are listed in the notes section.

| Option name            | Option description              | Default | Notes                            |
|------------------------|---------------------------------|---------|----------------------------------|
| ENABLE_TESTING         | Enable Testing system           | ON      |                                  |
| ENABLE_PAR             | Enable MPI support              | ON      | Requires MPI compiler            |
| ENABLE_CSC             | Enable coupling client          | ON      | Requires IMPACT                  |
| ENABLE_CFD_SUPPORT     | Enable offline CFD data I/O     | OFF     | Requires flann                   |

## Linux Build Instructions ##

First install project dependencies. For Ubuntu systems:
```
sudo apt install libarpack2-dev libspooles-dev libexodusii-dev libboost-serialization-dev libboost-iostreams-dev libflann-dev libspooles-dev libarpack2-dev
```
For CentOS systems (make sure EPEL repository is installed):

```
sudo yum install arpack-devel exodusii-devel boost-serialization boost-iostreams boost-devel flann-devel
```
You will need to build the IMPACT dependency from source:

* IMPACT (available at [here](https://github.com/IllinoisRocstar/IMPACT))

For instructions regarding building IMPACT please refer to [here](https://github.com/IllinoisRocstar/IMPACT). In case other libraries needed to be compiled from source:

* ARPACK (available [here](https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz))
* SPOOLES (available [here](http://www.netlib.org/linalg/spooles/spooles.2.2.tgz))

Now, you can compile the project:
```
$ export IMPACT_INSTALL_PATH=full_path_to_impact_install_dir
$ mkdir build
$ cd build
$ cmake .. \
-DCMAKE_PREFIX_PATH=${IMPACT_INSTALL_PATH} \
-DBUILD_SHARED_LIBS=ON 
$ make -j$(nproc) 
$ make install 
```
if you have manually compiled ARPACK and SPOOLES, and they are not installed on default system location then you can make specify the location of the projects before running cmake:
```
$ export SPOOLES_DIR=full_path_to_spooles_root_dir
$ export ARPACK_DIR=full_path_to_arpack_root_dir
$ cmake .. 
$ make -j$(nproc) 
$ make install 
```
then, you can test the project by:
```
make test
```
and tests pass, you can install by:
```
make install
```

## Building ARPACK ##

**NOTE:** These instructions are adopted from the standard ARPACK installation guide.

Before you can compile anything, you must first edit and correct the file `ARmake.inc`. Sample ARmake.inc's can be found in the ARMAKES directory. We recommend at least following changes:
* set the definition `home` to the root of the source tree (Top level of ARPACK directory).
* set `PLAT` to `PLAT = x64`
* set `FFLAG` to `FFLAG = -O3`
* make sure `MAKE` is set to full path of the `make` command (find full path by typing `which make` in shell)

Once changes are applied make the library by running
```
$ make lib
```
Depending on your environment, additional adjustments can be performed. Please consult `README` in the ARPACK folder for more details.


## Building SPOOLES ##

**NOTE:** These instructions are adopted from the standard SPOOLES installation guide.

Before you can compile anything, you must first edit and correct the file `Make.inc` and `timings.h`. We recommend at least following changes:

* uncomment and set `CC` to the full path of your `mpicc` compiler
* set `OPTLEVEL` to `OPTLEVEL = -O3`
* set `CFLAGS` to `CFLAGS = $(OPTLEVEL) -fPIC`
* set `MPI_INSTALL_DIR` to the fullpath location of your MPI library (for Ubuntu users with OpenMPI it will be `/usr/lib/x86_64-linux-gnu/openmpi/`)
* set `MPI_LIB_PATH` to empty `MPI_LIB_PATH =`
* set `MPI_LIBS` to `MPI_LIBS = $(MPI_LIB_PATH) -lpthread`

Change `timings.h` to:
```
#ifndef _TIMINGS_
#define _TIMINGS_
#include <sys/time.h>
static struct timeval  TV ;
//static struct timezone TZ ;
#define MARKTIME(t) \
gettimeofday(&TV, NULL) ; \
t = (TV.tv_sec + 0.000001*TV.tv_usec)
#endif
```
Once changes are applied make the library by running
```
$ make lib
```
Note, this option requires the Perl language installed on the system. Note that this option makes the serial version of the library. We recommend Multithread (MT) or MPI version of the library for use with CalculiX. To make the MT version, in the base SPOOLES folder run:
```
$ cd MT/src
$ make -f makeGlobalLib
```
for MPU version, in the base SPOOLES folder run:
```
$ cd MPI/src
$ make -f makeGlobalLib
```
Depending on your environment, additional adjustments can be performed. Please consult SPOOLES documentation (available [here](http://www.netlib.org/linalg/spooles/Install.ps.gz)) for details.
