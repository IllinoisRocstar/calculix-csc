#include "Rocout.h"
#include "com.h"
#include "com_devel.hpp"
#include "mpi.h"
#include <iostream>

COM_EXTERN_MODULE(clcxcsc);
COM_EXTERN_MODULE(Rocout);
COM_EXTERN_MODULE(Rocin);

// static vars
int _wrank = 0;
int _wsize = 0;
int _vrb = 1;
int _hndl;
std::string _jn = "case_name";
bool _step = false;
bool _prep = false;
std::string _prep_task = "rocstar";
double _final_time = 0.;
double _step_time = 0.;
bool _rstrt = false;
std::string _srf_fn = "";
std::string _vol_fn = "";
#ifdef HAVE_CFD      
bool _cfd_data = false;
std::string _cfd_data_fn = "";
#endif

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
  COM_set_default_communicator(MPI_COMM_WORLD);

  // loading
  COM_LOAD_MODULE_STATIC_DYNAMIC(clcxcsc, "clcx");

  // set job name
  _hndl = COM_get_function_handle("clcx.set_jobName");
  if (_hndl < 0)
    exitme("Error obtaining set_jobName handle");
  COM_call_function(_hndl, &_jn);

  if (_prep) {

    // set CFD dataset if requested
    // it is assumed CFD data are steady state and converged
#ifdef HAVE_CFD      
    // setting CFD data file name if needed
    if (_cfd_data) {
      // handle to driver friendly initialize method
      _hndl = COM_get_function_handle("clcx.set_cfd_data");
      if (_hndl < 0)
        exitme("Error obtaining set_cfd_data handle.\n");
      COM_call_function(_hndl, &_cfd_data_fn);
    }
#endif

    // preprocessing
    _hndl = COM_get_function_handle("clcx.preprocess");
    if (_hndl < 0)
      exitme("Error obtaining preprocess handle.\n");
    COM_call_function(_hndl, &_prep_task);

    // writing rocstar files
    if (_prep_task.find("rocstar") != std::string::npos) {
      std::cout << "Writing the window\n";
      COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
      int OUT_set = COM_get_function_handle("OUT.set_option");
      int OUT_write = COM_get_function_handle("OUT.write_dataitem");
      int OUT_write_control = COM_get_function_handle("OUT.write_rocin_control_file");
      int IN_all = COM_get_dataitem_handle("clcx_vol.all");
      // std::cout << "IN_all handle = " << IN_all << std::endl;
      char time_level[33] = "0";
      COM_call_function(OUT_write, "./test_vol_", &IN_all, "clcx", time_level);
      COM_call_function(OUT_write_control, "clcx_vol", "./test_vol_", "./solid_in_00.000000.txt");
      std::cout << "Finished writting volume window" << std::endl;
      IN_all = COM_get_dataitem_handle("clcx_srf.all");
      // std::cout << "IN_all handle = " << IN_all << std::endl;
      COM_call_function(OUT_write, "./test_srf_", &IN_all, "clcx", time_level);
      COM_call_function(OUT_write_control, "clcx_srf", "./test_srf_", "./isolid_in_00.000000.txt");
      std::cout << "Finished writting surface window" << std::endl;
      COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
      std::cout << "Unloaded SIMOUT" << std::endl;
    }

  } else if (_rstrt) {

    // mimicing agent's initialization process
    // load SimIN and load snapshots to input windows
    std::string srfWinName = "srf_in";
    std::string volWinName = "vol_in";
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    int _read_hndl = COM_get_function_handle("IN.read_window");
    if (_read_hndl < 0)
      exitme("Error obtaining read_window handle");

    // loading surface and volume windows
    COM_call_function(_read_hndl, (_srf_fn + " *_srf_*").c_str(),
                      srfWinName.c_str());
    COM_call_function(_read_hndl, (_vol_fn + " *_vol_*").c_str(),
                      volWinName.c_str());

    // handle to driver friendly initialize method
    _hndl = COM_get_function_handle("clcx.initialize");
    if (_hndl < 0)
      exitme("Error obtaining initialize handle.\n");
    double initT = 0.0;
    MPI_Comm comm = MPI_COMM_WORLD;
    int initHndl = -1;
    int obtHndl = -1;
    COM_call_function(_hndl, &initT, &comm, &initHndl, srfWinName.c_str(),
                      volWinName.c_str(), &obtHndl);

    if (_step) {
      // update_solution
      _hndl = COM_get_function_handle("clcx.update_solution");
      if (_hndl < 0)
        exitme("Error obtaining run handle.\n");
      double currTime = 0.0;
      int updHndl = -1;
      COM_call_function(_hndl, &currTime, &_step_time, &updHndl);
    }

    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");

  } else {

    // old initializer
    // handle to driver friendly initialize method
    _hndl = COM_get_function_handle("clcx.initialize_drv");
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
  //std::cout << "Process " << _wrank << " passed the barrier." << std::endl;
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
    std::cout
        << "\t" << argv[0] << " [[-flags] [flag_dependent_input(s)]]"
        << "\n\nwhere flags are \n"
        << "\t-v\n"
        << "\tverbose output\n"
        << "\t-e\n"
        << "\texodus output\n"
        << "\t-i\n"
        << "\tinput file name without .inp [input_name]\n"
        << "\t-p\n"
        << "\tpreprocess the task [task]\n"
        << "\ttask can be \"rocstar\" for generating Rocstar input "
           "data\n"
        << "\ttask can be \"fsi\" for the addition of Surface keyword to input "
           "deck\n"
        << "\t-s\n"
        << "\tstarts and steps for a given amount of time passed "
           "[target_time]\n"
        << "\t-f\n"
        << "\tproceed to the final simulation time passed [final_time]\n"
        << "\tthis switch only works when -s is used\n"
        << "\t-r\n"
        << "\tsurface file name [surf_file] and volume file name [vol_file]\n"
        << "\trestarts the from snapshot files and continues the simulation\n"
#ifdef HAVE_CFD      
        << "\t-cfd-data\n"
        << "\tuse cfd-data passed in file name [cfd_file]\n"
#endif
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
      if (args.compare("-i") == 0)
        _jn = std::string(argv[iArg + 1]);
      if (args.compare("-v") == 0)
        _vrb = 1;
      if (args.compare("-s") == 0) {
        _step = true;
        _step_time = std::stod(std::string(argv[iArg + 1]));
      }
      if (args.compare("-f") == 0)
        _final_time = std::stod(std::string(argv[iArg + 1]));
      if (args.compare("-p") == 0) {
        _prep = true;
        _prep_task = std::string(argv[iArg + 1]);
      }
      if (args.compare("-r") == 0) {
        _rstrt = true;
        _srf_fn = std::string(argv[iArg + 1]);
        _vol_fn = std::string(argv[iArg + 2]);
      }
#ifdef HAVE_CFD      
      if (args.compare("-cfd-data") == 0) {
        _cfd_data = true;
        _cfd_data_fn = std::string(argv[iArg + 1]);
      }
#endif
    }
  }
}
