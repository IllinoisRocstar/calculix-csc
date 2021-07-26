#ifndef CLCXCSC_H
#define CLCXCSC_H

#include "mpi.h"
#include <memory>
#include <string>

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
  // accordingly. This initialize method is for driver user only
  // not for restart operations.
  void initialize_drv(int vrb);

  // initializes the calculix from restart sanpshot windows
  // This initialize method is for rocstar/agent use. SimIN
  // first have to load in the surface and volume snapshot files
  // and pass the IN window names to this method. For full reference to
  // argument list please refer to Rocman_design.pdf available in IMPACT
  // project.
  void initialize(double &initTime, MPI_Comm &inComm, int &manInitHndl,
                  char *srfWinIn, char *volWinIn, int &obtHndl);

  // updates the calculix solution for the time interval timeStep and
  // subcycles if necessary. It will call a updateHndl (update_inbuff_X)
  // provided by Rocman before each subcycle.
  // At the end, it will update the outgoing interface buffers.
  void update_solution(double &currTime, double &timeStep, int &updHndl);

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

  // sets the file name with FSI CFD data
  // the input file should have pressures at first column, and
  // x,y, and z coordinates at the columns 2-4. The files should
  // be CSV formated.
  // this operation is currently only supported for driver/testing purposes
  void set_cfd_data(std::string cfd_data_fn);

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

private:
  std::vector<std::string> tokenize(const char *lineIn, const char &delim);

  // internal management
private:
  bool main_process() { return (_rank == 0 ? true : false); }
  void warm_up(int vrb = 0);
  void set_rank(int irnk = 0) { _rank = irnk; };
  void init_clcx();
  void step_initialize();
  void step_check();
  void step_run(const double step_tt);
  bool step_ready_to_finalize();
  void step_finalize();
  // scan for FSI node sets and registers them
  void scan_actions();
  void scan_fsi();
  void scan_fsi_2();
  void scan_elements();
  void register_global_data();
  void register_volume_data();
  void register_surface_data();
  void fsi_augment_inp();
  void fsi_compute_norms();
  void fsi_compute_face_centers();
  void fsi_sync_loads();
  void fsi_sync_other_quantities();
  void fsi_sync_loads_from_file();
  void debug_print(const std::string &str, int pane, int pe, MPI_Comm comm,
                   const char *memo);

public:
  // for external access to the update handler
  void update_bc(double dt);

  // management level variables
private:
  int _vrb;
  int _rank;
  int _nproc;
  std::string _wname;
  std::string _wname_vol;
  std::string _wname_srf;
  MPI_Comm _comm;
  bool _cfd_data;
  std::string _cfd_data_fn;

  // calculix interface
private:
  std::shared_ptr<clcx_interface> _ci;
  std::map<std::string, std::vector<size_t>> _elmTypIdx;
  std::vector<size_t> _tetConn;
  std::vector<size_t> _hexConn;
  // fsi surfaces
  size_t _nFsiNde;
  size_t _nFsiTri;
  size_t _nFsiQuad;
  size_t _nFsiFct;
  size_t _totElmNum;
  size_t _lowTetElmNum;  // to address element index offset 
  size_t _lowHexElmNum;  // to address element index offset
  std::map<std::string, std::vector<size_t>> _fsiNSetIdx;
  std::map<std::string, std::vector<size_t>> _elmSetIdx;
  std::map<std::string, std::vector<size_t>> _srfSideIdx;
  std::multimap<size_t, size_t> _fsiTetFaceIdx;
  std::multimap<size_t, size_t> _fsiHexFaceIdx;
  std::vector<size_t> _fsiNdeIdx;
  std::vector<double> _fsiNdeCrd;
  std::vector<int> _fsiTriConn;
  std::vector<int> _fsiQuadConn;
  std::vector<int> _fsiTriConnInterface;
  std::vector<int> _fsiQuadConnInterface;
  std::vector<double> _fsiTriNorm;
  std::vector<double> _fsiTriCenCrd;
  std::vector<double> _fsiU;
  std::vector<double> _fsiUHat;
  std::vector<double> _fsiV;
  std::vector<double> _fsiCs;
  std::vector<double> _fsiTs;
  std::vector<double> _fsiRhos;
  std::vector<double> _fsiTs_alp;
  std::vector<double> _fsiVbar_alp;
  std::vector<double> _fsiQs_alp;
  std::map<size_t, size_t> _fsiClcxLoadMap;
  std::vector<int> _fsiLoadArray;
  std::vector<std::vector<double>> _fsiClcxCenters;
  std::vector<std::vector<double>> _fsiSurfaceCenters;
  int _fsibcflag;
  int *_updHndl;
  double _upd_start_time;
  double _upd_time_step;
};

#endif // CLCXCSC_H
