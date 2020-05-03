#include "com.h"
#include "com_devel.hpp"
#include "Rocout.h"
#include "mpi.h"
#include <iostream>

COM_EXTERN_MODULE(clcxcsc);
COM_EXTERN_MODULE(Rocout);

// static vars
int _wrank = 0;
int _wsize = 0;
int _vrb = 1;
int _hndl;
std::string _jn = "case_name";
bool _step = false;
bool _prep = false;
double _final_time = 0.;
double _step_time = 0.;

// aux declarations
void usage(char **);
void exitme(std::string msg = std::string());
void procArgs(int argc, char *arg[]);

// main driver program
int main(int argc, char *argv[]) {

  // MPI communicator setup
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &_wsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &_wrank);

  // sanity checking
  if (argc < 2)
    usage(argv);
  else
    procArgs(argc, argv);

  // initializing COM
  COM_init(&argc, &argv);

  // loading
  COM_LOAD_MODULE_STATIC_DYNAMIC(clcxcsc, "clcx");

  // set job name
  _hndl = COM_get_function_handle("clcx.set_jobName");
  if (_hndl < 0)
    exitme("Error obtaining set_jobName handle");
  COM_call_function(_hndl, &_jn);

  if (_prep) {
      // preprocess only
      std::string task = "";
      _hndl = COM_get_function_handle("clcx.preprocess");
      if (_hndl < 0)
        exitme("Error obtaining preprocess handle.\n");
      COM_call_function(_hndl, &task);
      // writing window to CGNS
      std::cout << "Writing the window\n";
      COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
      int OUT_set = COM_get_function_handle("OUT.set_option");
      int OUT_write = COM_get_function_handle("OUT.write_dataitem");
      COM_call_function(OUT_set, "format", "HDF");
      int IN_all = COM_get_dataitem_handle("clcx_vol.all");
      std::cout << "IN_all handle = " << IN_all << std::endl;
      char time_level[33] = "";
      COM_call_function(OUT_write, "./test", &IN_all, "clcx",time_level);
      COM_call_function(OUT_set, "mode", "a");
      COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
  } else {
      // initialize
      _hndl = COM_get_function_handle("clcx.initialize");
      if (_hndl < 0)
        exitme("Error obtaining initialize handle.\n");
      COM_call_function(_hndl, _vrb);

      if (!_step) {
        // run
        _hndl = COM_get_function_handle("clcx.run");
        if (_hndl < 0)
          exitme("Error obtaining run handle.\n");
        COM_call_function(_hndl);
      } else {

        // change final simulation time if needed
        if (_final_time > 0) {
          _hndl = COM_get_function_handle("clcx.set_final_time");
          if (_hndl < 0)
            exitme("Error obtaining step handle.\n");
          COM_call_function(_hndl, &_final_time);
        }

        // step in time
        _hndl = COM_get_function_handle("clcx.step");
        if (_hndl < 0)
          exitme("Error obtaining step handle.\n");
        COM_call_function(_hndl, &_step_time);
      }
  }

  // finalize
  _hndl = COM_get_function_handle("clcx.finalize");
  if (_hndl < 0)
    exitme("Error obtaining finalize handle.\n");
  COM_call_function(_hndl);

  // unloading
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(clcxcsc, "clcx");

  // finalizing COM
  COM_finalize();

  // peaceful exit
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

// aux implementations
void usage(char *argv[]) {
  if (_wrank == 0) {
    std::cout << "Calculix Component Side Client (CSC) Driver" << std::endl;
    std::cout << "NOTE: This is an experimental interface to the Calculix and "
                 "subject to change\n";
    if (_wsize > 1)
      std::cout << "Parallel run with " << _wsize << " processes.\n";
    std::cout << std::endl << "Usage: " << std::endl;
    std::cout << "\t" << argv[0] << " [[-flags] [flag_dependent_value]]"
              << "\n\nwhere flags are \n"
              << "\t-v\n"
              << "\tverbose output\n"
              << "\t-e\n"
              << "\texodus output\n"
              << "\t-c\n"
              << "\tcase file name [case_name]\n"
              << "\t-p\n"
              << "\tpre-process and stop\n"
              << "\t-s\n"
              << "\tstarts and steps for a given amount of time passed "
                 "[target_time]\n"
              << "\t-f\n"
              << "\tproceed to the final simulation time passed [final_time]\n"
              << "\tthis switch only works when -s is used\n"
              << std::endl
              << std::endl
              << std::endl;
  }
  exitme("Completed!");
}

void exitme(std::string msg) {
  if (_wrank == 0) {
    std::cerr << "Calculix-Driver : " << msg << std::endl;
  }
  MPI_Finalize();
  exit(-1);
}

void procArgs(int argc, char *argv[]) {
  if (argc > 2) {
    for (int iArg = 1; iArg < argc; iArg++) {
      std::string args(argv[iArg]);
      if (args.compare("-c") == 0)
        _jn = std::string(argv[iArg + 1]);
      if (args.compare("-v") == 0)
        _vrb = 1;
      if (args.compare("-s") == 0) {
        _step = true;
        _step_time = std::stod(std::string(argv[iArg + 1]));
      }
      if (args.compare("-f") == 0)
        _final_time = std::stod(std::string(argv[iArg + 1]));
      if (args.compare("-p") == 0)
        _prep = true;
    }
  }
}
