CalculiX -- IR Extras
---------------------
This is the vanila CalculiX project (forked from v2.15) with a few improvmenets and additions as follows:
* IMPACT-based multiphysics coupling module (CSC)
* CMake-based build system
* Homogeniziation module
* ExodusII output and GPU solver support (GPU is not completely tested, using patches released by Peter A. Gustafson on 8/14/2019).

## Version
Version 2.15.2

The project follows the versioning convention of the original CalculiX. The versioning will be major.minor.patch and increaments are as following:
* major: shows major version of the CalculiX project 
* minor: shows minor version of the CalculiX project
* patch: will be increamented as we release more functionality to the project

## Options
The following table contains all available CMake options to configure the project. The necessary third-party libraries are listed in the notes section.

| Option name            | Option description              | Default | Notes                            |
|------------------------|---------------------------------|---------|----------------------------------|
| ENABLE_TESTING         | Enable Testing system           | ON      |                                  |
| ENABLE_PAR             | Enable MPI support              | ON      | Requires MPI compiler            |
| ENABLE_CSC             | Enable coupling client          | ON      | Requires IMPACT                  |
| ENABLE_CFD_SUPPORT     | Enable offline CFD data I/O     | OFF     | Requires flann                   |

## Linux Build Instructions ##
First install project depencies (on Ubuntu):

* libspooles-dev
* libexodusii-dev
* libboost-serialization-dev
* libboost-iostreams-dev
* IMPACT (available at https://github.com/IllinoisRocstar/IMPACT)
* ARPACK (available at https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz)
* SPOOLES (available at http://www.netlib.org/linalg/spooles/spooles.2.2.tgz)
* libflann-dev (optional, only if `ENABLE_CFD_SUPPORT=ON`)

For instructions regarding building IMPACT please refer to [here](https://github.com/IllinoisRocstar/IMPACT).

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
if you have manually compiled ARPACK and SPOOLES, and they are not installed on defalut system location then you can make specify the location of the projects before running cmake:

```
$ export SPOOLES_DIR=full_path_to_spooles_root_dir
$ export ARPACK_DIR=full_path_to_arpack_root_dir
$ cmake .. 
$ make -j$(nproc) 
$ make install 
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

**NOTE:** On Ubuntu systems, simply install spooles available from system repository.
**NOTE:** These instructions are adopted from the standard SPOOLES installation guide.


Before you can compile anything, you must first edit and correct the file `Make.inc`. We recommend at least following changes:
* uncomment and set `CC` to the full path of your `mpicc` compiler
* set `OPTLEVEL` to `OPTLEVEL = -O3`
* set `CFLAGS` to `CFLAGS = $(OPTLEVEL) -fPIC`
* set `MPI_INSTALL_DIR` to the fullpath location of your MPI library (for Ubuntu users with OpenMPI it will be `/usr/lib/x86_64-linux-gnu/openmpi/`)
* set `MPI_LIB_PATH` to empty `MPI_LIB_PATH =`
* set `MPI_LIBS` to `MPI_LIBS = $(MPI_LIB_PATH) -lpthread`
 
Once changes are applied make the library by running
```
$ make lib
```
Note, this option requires the Perl language installed on the system. Note that this option makes the serial version of the library. If we recomment Multithread (MT) or MPI version of the library for use with CalculiX. To make the MT version, in the base SPOOLES folder run:
```
$ cd MT/src
$ make -f makeGlobalLib
```

for MPU version, in the base SPOOLES folder run:

```
$ cd MPI/src
$ make -f makeGlobalLib
```

Depending on your environment, additional adjustments can be performed. Please consult SPOOLES documentation (here: http://www.netlib.org/linalg/spooles/Install.ps.gz) for details.
