#ifndef CLCXCSC_H
#define CLCXCSC_H

#include <string>
#include <memory>

// IMPACT
#include "com.h"

// Calculix inteface
#include "clcxInterface.H"

class clcx_module : public COM_Object {
public:
  clcx_module();
  ~clcx_module(){};

public:
  static void Load(const std::string &name);
  static void Unload(const std::string &name);

  // top level methods used to control the simulation
  // workflow. These methods are registered on the COM
  // window object
public:
  // sets the name of the job
  void set_jobName(std::string);

  // initializes the calculix from scratch by
  // reading the input file and loading native data structure
  // accordingly
  void initialize(int vrb);

  // runs the simulation end-to-end without applying
  // any control. This method has full access to all
  // possible combination of simulation configurations.
  void run();

  // finalizes the simulation
  void finalize();

  // sets the simulation final time. Simulation time can not exceed
  // this value. This method will overide the final time prescribed
  // by the input file.
  // NOTE: currently this method should be called only after initialize
  //   , otherwise the default value (controlled from input file)
  //   will be used
  void set_final_time(const double &t_f);

  // NOTE: this method is only tested with certain simulation configurations
  // steps the simulation for the requested amount of time (\delta time).
  // If user has set the input file properly (with direct integration)
  // the simulation steps in the time accurately for the given amount of time,
  // otherwise the value is used as target and simulation proceed to its
  // neighborhood. If target time goes beyond what prescribed at input file
  // user will be informed. Input file will take the precedence. In order to
  // go beyond the input file value use set_final_time command.
  // in case of perturbation simulations, only one step is taken and value
  // passed is not used
  void step(const double &step_tt);

  // performs a variety of pre-processing actions
  // currently only generating GENX data files
  void preprocess(std::string actStr);

  // helper methods
public:
  // writes out message on processs rank 0
  int message(std::string msg, bool enforce_vrb = false);

  // set output type
  void set_outputType(std::string = "exo");

  // write splash to stdout
  void splash();

  // perform post-initialization steps (if multiple steps needed should be
  // called)
  void postInit();

  // check if step is completed
  bool step_completed();

  // internal management
private:
  void warm_up(int vrb = 0);
  void init_clcx();
  void step_initialize();
  void step_check();
  void step_run(const double step_tt);
  bool step_ready_to_finalize();
  void step_finalize();
  void register_global_data();
  void register_volume_data();
  void register_surface_data();

  // management level variables
private:
  int _vrb;
  int _rank;
  int _nproc;
  std::string _wname;
  std::string _wname_vol;
  std::string _wname_srf;
  MPI_Comm _comm;

// calculix interface
private:
  std::shared_ptr<clcx_interface> _ci;

};

#endif // CLCXCSC_H
