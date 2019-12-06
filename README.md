CalculiX -- IR Extras
---------------------
This is plain, vanila CalculiX project (v2.15) with a few additions as follows:
* IMPACT coupling module (CSC)
* Totally revamped build system (based on CMake)
* CUDA based solvers (currently in-active) 
* ExodusII output format 
The CUDA and ExodusII support is made available using patches released by Peter A. Gustafson (last update 8/14/2019).

## Version
Version 2.15.1

The project follows the versioning convention of the original CalculiX. The versioning will be major.minor.patch and increaments are as following:
* major: shows major version of the CalculiX project 
* minor: shows minor version of the CalculiX project
* patch: will be increamented as we release more functionality to the project

## Options
The following table contains all available CMake options to configure the project. The necessary third-party libraries are listed in the notes section.

| Option name            | Option description              | Default | Notes                            |
|------------------------|---------------------------------|---------|----------------------------------|
| ENABLE_PAR             | Enable MPI support              | ON     | Requires MPI compiler            |
| ENABLE_CSC         | Enable coupling client                  | ON      | Requires IMPACT                            |

## Unix Building Instructions ##
First install project depencies (on Ubuntu):

* libspooles-dev
* libexodusii-dev
* IMPACT

For instructions regarding building IMPACT please refer to [here](https://github.com/IllinoisRocstar/IMPACT).

Now, you can compile the project:

```
export IMPACT_INSTALL_PATH=full_path_to_impact_install_dir
mkdir build
cd build
cmake .. \
    -DCMAKE_PREFIX_PATH=${IMPACT_INSTALL_PATH} \
make -j$(nproc) (or however many threads you'd like to use)
$ make install (sudo if install location requires it)
```

