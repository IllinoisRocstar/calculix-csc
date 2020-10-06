#include "clcxCSC.H"
#include "stepNonlinGeo.H"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#ifdef HAVE_CFD
#include <flann/flann.hpp>
#endif

void clcx_module_upd_bc_wrapper(void *context, double dt) {
  static_cast<clcx_module *>(context)->update_bc(dt);
}

// defaults constructor
clcx_module::clcx_module()
    : _vrb(0), _rank(0), _nproc(1), _cfd_data(false), _nFsiNde(0), _nFsiTri(0),
      _nFsiQuad(0), _nFsiFct(0), _totElmNum(0), _lowTetElmNum(0),
      _lowHexElmNum(0), _updHndl(nullptr) {
  _ci = std::make_shared<clcx_interface>();
  _ci->register_update_bc(&clcx_module_upd_bc_wrapper, this);
}

// COM module loader
void clcx_module::Load(const std::string &name) {

  // anouncing default communicator
  MPI_Comm inComm;
  inComm = COM_get_default_communicator();

  // TODO: at some point we need to split
  // communicator here.
  int rank, nproc;
  MPI_Comm_rank(inComm, &rank);
  MPI_Comm_size(inComm, &nproc);

  // register module with COM
  clcx_module *module_pointer = new clcx_module();
  COM_new_window(name, MPI_COMM_NULL);
  module_pointer->_wname = name;
  module_pointer->set_rank(rank);
  MPI_Comm_dup(inComm, &(module_pointer->_comm));
  MPI_Comm_rank(module_pointer->_comm, &(module_pointer->_rank));
  MPI_Comm_size(module_pointer->_comm, &(module_pointer->_nproc));
  std::string global_name(name + ".global");
  COM_new_dataitem(global_name.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object(global_name.c_str(), 0, module_pointer);

  // sanity check
  module_pointer->message("Loading on " + std::to_string(nproc) +
                          " processes.");
  module_pointer->message("Loading Calculix on the window " + name);

  // registering functions
  COM_Type types[7];

  types[0] = COM_RAWDATA;
  types[1] = COM_INT;
  COM_set_member_function((name + ".initialize_drv").c_str(),
                          (Member_func_ptr)(&clcx_module::initialize_drv),
                          global_name.c_str(), "bi", types);

  COM_set_member_function((name + ".run").c_str(),
                          (Member_func_ptr)(&clcx_module::run),
                          global_name.c_str(), "b", types);

  COM_set_member_function((name + ".finalize").c_str(),
                          (Member_func_ptr)(&clcx_module::finalize),
                          global_name.c_str(), "b", types);

  types[1] = COM_STRING;
  COM_set_member_function((name + ".set_jobName").c_str(),
                          (Member_func_ptr)(&clcx_module::set_jobName),
                          global_name.c_str(), "bi", types);

  types[1] = COM_DOUBLE;
  COM_set_member_function((name + ".set_final_time").c_str(),
                          (Member_func_ptr)(&clcx_module::set_final_time),
                          global_name.c_str(), "bi", types);

  types[1] = COM_DOUBLE;
  COM_set_member_function((name + ".step").c_str(),
                          (Member_func_ptr)(&clcx_module::step),
                          global_name.c_str(), "bi", types);

  types[1] = COM_STRING;
  COM_set_member_function((name + ".preprocess").c_str(),
                          (Member_func_ptr)(&clcx_module::preprocess),
                          global_name.c_str(), "bi", types);

  types[1] = COM_STRING;
  COM_set_member_function((name + ".set_cfd_data").c_str(),
                          (Member_func_ptr)(&clcx_module::set_cfd_data),
                          global_name.c_str(), "bi", types);

  // initialize(double initTime, MPI_Comm inComm, int manInitHndl, char
  // *srfWinIn, char *volWinIn, int obtHndl)
  types[0] = COM_RAWDATA;
  types[1] = COM_DOUBLE;
  types[2] = COM_MPI_COMM;
  types[3] = COM_INT;
  types[4] = COM_STRING;
  types[5] = COM_STRING;
  types[6] = COM_INT;
  COM_set_member_function((name + ".initialize").c_str(),
                          (Member_func_ptr)(&clcx_module::initialize),
                          global_name.c_str(), "biiiiii", types);

  // update_solution(double& currTime, double& timeStep, int& updateHndl);
  types[0] = COM_RAWDATA;
  types[1] = COM_DOUBLE;
  types[2] = COM_DOUBLE;
  types[3] = COM_INT;
  COM_set_member_function((name + ".update_solution").c_str(),
                          (Member_func_ptr)(&clcx_module::update_solution),
                          global_name.c_str(), "biii", types);

  COM_window_init_done(name);
}

// COM module unloader
void clcx_module::Unload(const std::string &name) {
  clcx_module *module_pointer = NULL;
  std::string global_name(name + ".global");
  COM_get_object(global_name.c_str(), 0, &module_pointer);
  COM_assertion_msg(module_pointer->validate_object() == 0, "Invalid object");

  module_pointer->message("Unloading Calculix from " + name);
  delete module_pointer;
  COM_delete_window(std::string(name));
}

// C/C++ bindings to load Calculix
extern "C" void clcxcsc_load_module(const char *name) {
  clcx_module::Load(name);
}

// C/C++ bindings to unload Calculix
extern "C" void clcxcsc_unload_module(const char *name) {
  clcx_module::Unload(name);
}

int clcx_module::message(std::string msg, bool enfv) {

  if (_rank == 0 && (!enfv || (_vrb > 0 && enfv)))
    std::cout << "Calculix : " << msg << std::endl;

  return (1);
}

void clcx_module::set_jobName(std::string jn) {
  message("Job name was set to " + jn);
  strcpy(_ci->jobnamec, jn.c_str());
  jn = jn + " ";
  strcpy(_ci->jobnamef, jn.c_str());
}

void clcx_module::set_outputType(std::string oti) {
  std::string ot;
  std::transform(oti.begin(), oti.end(), std::back_inserter(ot), ::tolower);
  if (ot.compare("exo") == 0) {
    message("Output file was set to " + ot);
    strcpy(_ci->output, "exo");
  } else {
    strcpy(_ci->output, "asc");
  }
}

void clcx_module::set_final_time(const double &t_f) { _ci->timepar[1] = t_f; }

void clcx_module::warm_up(int vrbin) {
  _vrb = vrbin;
  _ci->set_verb(vrbin);
  _ci->myid = _rank;
  _ci->nproc = _nproc;

  // window names for internal use
  _wname_vol = _wname + "_vol";
  _wname_srf = _wname + "_srf";
}

void clcx_module::initialize_drv(int vrb) {
  // TODO: There are several hard coded variables with funky values.
  //       Some work needed to tidy up the implementation.
  warm_up(vrb);
  message("Initializing ...");

  // initializing calculix from native input file
  // this has to be refactored. We only route to this if
  // for some reason GENX interface files are not available
  _ci->init();
}

void clcx_module::run() {
  message("Running ...");
  message("Method = " + std::to_string(_ci->nmethod));
  message("Done running");
}

void clcx_module::finalize() {
  message("Finilizing ...");
  if (main_process())
    _ci->finalize();
  message("Done finalizing");
}

void clcx_module::step(const double &stp_tt) { _ci->step(stp_tt); }

void clcx_module::preprocess(std::string actStr) {
  warm_up();
  message("Action is " + actStr);
  message("Pre-processing....");

  if (main_process()) {
    // read input file and generate respective data structure
    _ci->init();

    // preparation operations
    scan_actions();
  }

  if (actStr.find("rocstar") != std::string::npos) {

    if (main_process()) {
      // read CFD data if requested
      if (_cfd_data)
        fsi_sync_loads_from_file();
      else
        // setting FSI loads to zero
        fsi_sync_loads();
    }

    // register global data
    register_global_data();

    // register volume data
    register_volume_data();

    // register global data
    register_surface_data();

  } else if (actStr.find("fsi") != std::string::npos) {
    if (main_process()) {
      scan_fsi_2();
      fsi_augment_inp();
    }
  }
}

void clcx_module::register_global_data() {

  // registering solver state
  COM_new_window(_wname_vol, MPI_COMM_NULL);

  // state
  auto ss = std::make_shared<std::stringstream>();
  // auto ss = std::ofstream("test_of.dat");
  _ci->save(*ss);
  const std::string &ss_str = ss->str();
  // std::cout << "Size is " << ss_str.size() << std::endl;
  COM_new_dataitem(_wname_vol + ".state", 'w', COM_CHAR, 1, "");
  COM_set_size(_wname_vol + ".state", 0, ss_str.size());
  COM_resize_array(_wname_vol + ".state", 0);
  char *buf;
  COM_get_array((_wname_vol + ".state").c_str(), 0, &buf);
  memcpy(buf, ss_str.c_str(), ss_str.size());

  // for future debugging
  // std::ofstream fo;
  // fo.open("pre.dat");
  // fo << ss_str.c_str();
  // fo.close();
  // std::ifstream fi;
  // fi.open("pre.dat");
  // std::stringstream ss2;
  // ss2 << fi.rdbuf();
  // fi.close();
  //_ci->load(ss2);

  COM_window_init_done(_wname_vol);
}

void clcx_module::register_volume_data() {
  // nodal coordinates
  if (main_process()) {
    message("Number of nodes " + std::to_string(_ci->nk));
    COM_set_size(_wname_vol + ".nc", _rank + 1, _ci->nk);
    COM_set_array(_wname_vol + ".nc", _rank + 1, _ci->co, 3);
  }
  //} else
  //  COM_set_size(_wname_vol + ".nc", _rank + 1, 0);

  // connectivity

  // supporting tets for now
  if (main_process()) {
    if (_elmTypIdx.find("C3D4") == _elmTypIdx.end() &&
        _elmTypIdx.find("C3D8") == _elmTypIdx.end()) {
      message("Only tetrahedral and hexahedral elements are supported!");
      message("Element types found");
      for (auto const &eti : _elmTypIdx)
        message(eti.first);
      throw;
    }
  }

  // registering all tet or hex elements under pane 1
  // note: we do not allow for mixed meshes yet
  size_t nTet = 0;
  size_t nHex = 0;
  if (main_process()) {
    nTet = _elmTypIdx["C3D4"].size();
    nHex = _elmTypIdx["C3D8"].size();
  }

  if (nTet > 0) {
    COM_set_size(_wname_vol + ".:T4", _rank + 1, nTet);
    COM_set_array(_wname_vol + ".:T4", _rank + 1, &_tetConn[0], 4);
  } else if (nHex > 0) {
    COM_set_size(_wname_vol + ".:H8", _rank + 1, nHex);
    COM_set_array(_wname_vol + ".:H8", _rank + 1, &_hexConn[0], 8);
  }

  // output variables
  // registering displacements (u)_{SN}^{n+1}
  COM_new_dataitem(_wname_vol + ".u", 'n', COM_DOUBLE, 3, "");
  if (main_process()) {
    COM_set_array(_wname_vol + ".u", _rank + 1, _ci->vold, 3);
  }
  // bound values for debugging
  double lb = 0.;
  double ub = 0.;
  for (size_t idx = 0; idx < (_ci->nk); idx++) {
    ub = std::max(ub, _ci->vold[idx]);
    lb = std::min(lb, _ci->vold[idx]);
  }
  message("Min vol.u value is " + std::to_string(lb));
  message("Max vol.u value is " + std::to_string(ub));

  // registering velocity (v_{s})_{SN}^{n+1}
  COM_new_dataitem(_wname_vol + ".vs", 'n', COM_DOUBLE, 3, "");
  if (main_process()) {
    COM_set_array(_wname_vol + ".vs", _rank + 1, _ci->veold, 3);
  }
  // bound values for debugging
  ub = 0.;
  lb = 0.;
  for (size_t idx = 0; idx < (_ci->nk); idx++) {
    ub = std::max(ub, _ci->veold[idx]);
    lb = std::min(lb, _ci->veold[idx]);
  }
  message("Min vol.vs value is " + std::to_string(lb));
  message("Max vol.vs value is " + std::to_string(ub));

  COM_window_init_done(_wname_vol);
}

// interface data registration
void clcx_module::register_surface_data() {

  // surface data registration
  COM_new_window(_wname_srf, MPI_COMM_NULL);

  if (_nFsiNde > 0) {
    // nodal coordinates
    COM_set_size(_wname_srf + ".nc", _rank + 1, _nFsiNde);
    COM_set_array(_wname_srf + ".nc", _rank + 1, &_fsiNdeCrd[0], 3);
  }

  // connectivity

  // registering all fsi facets under pane 1
  if (_nFsiTri > 0) {
    COM_set_size(_wname_srf + ".:t3", _rank + 1, _nFsiTri);
    COM_set_array(_wname_srf + ".:t3", _rank + 1, &_fsiTriConnInterface[0], 3);
  } else if (_nFsiQuad > 0) {
    COM_set_size(_wname_srf + ".:q4", _rank + 1, _nFsiQuad);
    COM_set_array(_wname_srf + ".:q4", _rank + 1, &_fsiQuadConnInterface[0], 4);
  }

  // physics quantities
  _fsiU.clear();
  _fsiUHat.resize(_nFsiNde * 3, 0.0);
  _fsiV.clear();
  _fsiCs.resize(_nFsiFct, 0.0);
  _fsiTs.resize(_nFsiNde, 0.0);
  _fsiRhos.resize(_nFsiFct, 0.0);
  _fsiTs_alp.resize(_nFsiFct * 3, 0.0);
  _fsiVbar_alp.resize(_nFsiNde * 3, 0.0);
  _fsiQs_alp.resize(_nFsiFct, 0.0);
  if (main_process())
    _fsibcflag = 0;
  else
    _fsibcflag = 2;

  // updating tractions from previous run
  for (int iFct = 0; iFct < _nFsiFct; iFct++) {
    // updating interface quantities
    double pressure = _ci->xload[iFct * 2];
    _fsiTs_alp[iFct * 3] = _fsiTriNorm[iFct * 3] * pressure;
    _fsiTs_alp[iFct * 3 + 1] = _fsiTriNorm[iFct * 3 + 1] * pressure;
    _fsiTs_alp[iFct * 3 + 2] = _fsiTriNorm[iFct * 3 + 2] * pressure;
  }

  // supporting only displacement and velocity
  // for (auto &ndeIdx : _fsiNSetIdx["C3D4"])
  for (auto &ndeIdx : _fsiNdeIdx)
    for (auto iCmp = 0; iCmp < 3; iCmp++) {
      _fsiU.push_back(_ci->vold[(ndeIdx - 1) * 3 + iCmp]);
      _fsiV.push_back(_ci->veold[(ndeIdx - 1) * 3 + iCmp]);
    }

  // output variables
  // registering displacements (u)_{SN}^{n+1}
  COM_new_dataitem(_wname_srf + ".u", 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem(_wname_srf + ".uhat", 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem(_wname_srf + ".vs", 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem(_wname_srf + ".Cs", 'e', COM_DOUBLE, 1, "");
  COM_new_dataitem(_wname_srf + ".Ts", 'n', COM_DOUBLE, 1, "");
  COM_new_dataitem(_wname_srf + ".rhos", 'e', COM_DOUBLE, 1, "");
  COM_new_dataitem(_wname_srf + ".bcflag", 'p', COM_INT, 1, "");
  if (_nFsiNde > 0) {
    // COM_set_size(_wname_srf + ".u", 1, _nFsiNde, 3);
    COM_set_array(_wname_srf + ".u", _rank + 1, &_fsiU[0], 3);

    // registering displacements (u_{hat})_{SN}^{n+1}
    // According to Rocman Design, it is not completely clear
    // what is the difference between u and u_{hat}. It seems
    // u_hat is reserved for some sort of particle which we are
    // not supporting yet
    // COM_set_size(_wname_srf + ".uhat", 1, _nFsiNde);
    COM_set_array(_wname_srf + ".uhat", _rank + 1, &_fsiUHat[0], 3);

    // registering velocity (v_{s})_{SN}^{n+1}
    // COM_set_size(_wname_srf + ".vs", 1, _nFsiNde);
    COM_set_array(_wname_srf + ".vs", _rank + 1, &_fsiV[0], 3);

    // registering temperature (T_s)_{SN}^{n+1}
    // COM_set_size(_wname_srf + ".Ts", 1, _nFsiNde);
    COM_set_array(_wname_srf + ".Ts", _rank + 1, &_fsiTs[0], 1);
  }

  if (_nFsiFct > 0) {
    // registering specific heat (C_s)_{SF}^{n+1}
    // COM_set_size(_wname_srf + ".Cs", 1, _nFsiTri);
    COM_set_array(_wname_srf + ".Cs", _rank + 1, &_fsiCs[0], 1);

    // registering density (rho_s)_{SF}^{n+1}
    // COM_set_size(_wname_srf + ".rhos", 1, _nFsiTri);
    COM_set_array(_wname_srf + ".rhos", _rank + 1, &_fsiRhos[0], 1);

    // registering bcflag
    COM_set_size(_wname_srf + ".bcflag", _rank + 1, 1);
    COM_set_array(_wname_srf + ".bcflag", _rank + 1, &_fsibcflag, 1);
  }

  //// input variables
  COM_new_dataitem(_wname_srf + ".vbar_alp", 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem(_wname_srf + ".ts_alp", 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem(_wname_srf + ".qs_alp", 'e', COM_DOUBLE, 1, "");

  // registering tractions (t_s)_{SN}^{n+alpha}
  // seems like tractions are passed for nodes according to the manual
  // rather than face centers this may need to be double checked

  if (_nFsiNde > 0) {
    // registering mesh velocity (vbar_alp)_{SN}^{n+alpha}
    // COM_set_size(_wname_srf + ".vbar_alp", 1, _nFsiTri);
    COM_set_array(_wname_srf + ".vbar_alp", _rank + 1, &_fsiVbar_alp[0], 3);
  }

  if (_nFsiFct > 0) {
    // COM_set_size(_wname_srf + ".ts_alp", 1, _nFsiTri);
    COM_set_array(_wname_srf + ".ts_alp", _rank + 1, &_fsiTs_alp[0], 3);

    // registering heat flux (q_s)_{SF}^{n}
    // annother inconsistency in the Rocman Design, the variable is called
    // qs_alpha but symbol show ^{n} instead of ^{n+alpha}
    // COM_set_size(_wname_srf + ".qs_alp", 1, _nFsiTri);
    COM_set_array(_wname_srf + ".qs_alp", _rank + 1, &_fsiQs_alp[0], 1);
  }

  COM_window_init_done(_wname_srf);
}

// TODO: this method is only supporting tetrahedral meshes
// searchs for fsi boundary conditions within nodesets. Any nodeset with fsi
// in the name will be recorded. This implementation is poor-man's solution
// to cases without surface sidesets. The implementation only works for
// tetrahedral elements. Known issue:
// * if a face share all fsi nodes, it will be picked up as a fsi surface
//   regardgless.
void clcx_module::scan_fsi_2() {
  message("Collecting FSI set information");
  std::vector<std::string> setNames(tokenize(_ci->set, ' '));

  // only node sets
  int iSet = 0;
  _fsiNSetIdx.clear();
  for (auto const &sn : setNames) {
    if (sn[sn.size() - 1] == 'N')
      if (sn.find("FSI") != std::string::npos) {
        std::vector<size_t> idx;
        for (int i = _ci->istartset[iSet] - 1; i < _ci->iendset[iSet]; i++)
          idx.push_back(_ci->ialset[i]);
        _fsiNSetIdx[sn] = idx;
      }
    iSet++;
  }
  message("Found " + std::to_string(_fsiNSetIdx.size()) + " FSI nodeset(s).");

  // total number of nodes in FSI patches
  _nFsiNde = 0;
  for (auto const &patch : _fsiNSetIdx)
    _nFsiNde += patch.second.size();
  message("Number of FSI nodes " + std::to_string(_nFsiNde));

  // coodinates of all fsi nodes
  for (auto const &patch : _fsiNSetIdx)
    for (auto const &ndeIdx : patch.second)
      _fsiNdeIdx.push_back(ndeIdx);

  // loop through all elements and find tri faces that share fsi nodes
  _nFsiTri = 0;
  std::set<size_t> fsiNdeIdxSet(_fsiNdeIdx.begin(), _fsiNdeIdx.end());
  std::vector<size_t> tetFaceConnPtrn = {0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0};
  _fsiTriConn.clear();
  _fsiTetFaceIdx.clear();
  for (size_t tetIdx = 0; tetIdx < _elmTypIdx["C3D4"].size(); tetIdx++) {
    for (size_t facIdx = 0; facIdx < 4; facIdx++) {
      bool faceIsOn = true;
      std::vector<size_t> tetFaceConnTmp;
      for (size_t ndeLoc = 0; ndeLoc < 3; ndeLoc++) {
        size_t faceNdeIdx =
            _tetConn[tetIdx * 4 + tetFaceConnPtrn[facIdx * 3 + ndeLoc]];
        tetFaceConnTmp.push_back(faceNdeIdx);
        if (fsiNdeIdxSet.find(faceNdeIdx) == fsiNdeIdxSet.end()) {
          faceIsOn = false;
          break;
        }
      }

      if (!faceIsOn)
        continue;

      // std::cout << "Tet Idx " << _elmTypIdx["C3D4"][tetIdx] << " face "
      //          << facIdx << " is on\n";
      //_fsiTetFaceIdx[_elmTypIdx["C3D4"][tetIdx]] = facIdx + 1;
      _fsiTetFaceIdx.insert(
          std::pair<size_t, size_t>(_elmTypIdx["C3D4"][tetIdx], facIdx + 1));
      _fsiTriConn.insert(_fsiTriConn.end(), tetFaceConnTmp.begin(),
                         tetFaceConnTmp.end());
    }
  }
  _nFsiTri = _fsiTriConn.size() / 3;
  message("Founds " + std::to_string(_nFsiTri) +
          " Tri facets on FSI surfaces.");

  // coodinates of all fsi nodes
  // renumbering face connectivities
  _fsiNdeCrd.clear();
  std::map<size_t, size_t> vol2srfConn;
  size_t newIndx = 1;
  for (auto const &indx : _fsiTriConn) {
    auto ret = vol2srfConn.insert(std::pair<size_t, size_t>(indx, newIndx));
    if (ret.second == true) {
      newIndx++;
      for (size_t iCrd = 0; iCrd < 3; iCrd++)
        _fsiNdeCrd.push_back(_ci->co[(indx - 1) * 3 + iCrd]);
    }
  }
  // std::cout << "FSI node coordinates " << std::endl;
  // for (int iNde = 0; iNde < _fsiNdeCrd.size() / 3; iNde++)
  //  std::cout << "Node " << iNde << " " << _fsiNdeCrd[iNde * 3] << " "
  //            << _fsiNdeCrd[iNde * 3 + 1] << " " << _fsiNdeCrd[iNde * 3 + 2]
  //            << std::endl;
  _fsiTriConnInterface.clear();
  for (auto const &indx : _fsiTriConn)
    _fsiTriConnInterface.push_back(vol2srfConn[indx]);

  // std::cout << "New connectivities : \n";
  // for (auto const &idx : _fsiTriConnInterface)
  //  std::cout << idx << std::endl;
  message("FSI surface connectivity number of elements is " +
          std::to_string(_fsiTriConnInterface.size() / 3));
}

// searchs for fsi boundary conditions within nodesets. Any nodeset with fsi
// in the name will be recorded
void clcx_module::scan_fsi() {
  message("Collecting FSI set information");
  std::vector<std::string> setNames(tokenize(_ci->set, ' '));

  // for (auto const &name : setNames)
  //  std::cout << name << std::endl;

  // only fsi element and sidesets considered
  int iSet = 0;
  _elmSetIdx.clear();
  _srfSideIdx.clear();
  for (auto const &sn : setNames) {
    if (sn[sn.size() - 1] == 'E') {
      std::vector<size_t> idx;
      for (int i = _ci->istartset[iSet] - 1; i < _ci->iendset[iSet]; i++)
        idx.push_back(_ci->ialset[i]);
      _elmSetIdx[sn] = idx;
    } else if (sn[sn.size() - 1] == 'T') {
      std::string elmSetName = sn;
      elmSetName.pop_back();
      std::vector<size_t> idx;
      for (int i = _ci->istartset[iSet] - 1; i < _ci->iendset[iSet]; i++)
        idx.push_back(_ci->ialset[i]);
      _srfSideIdx[elmSetName] = idx;
    }
    iSet++;
  }
  message("Found " + std::to_string(_elmSetIdx.size()) + " element set(s).");
  message("Found " + std::to_string(_srfSideIdx.size()) + " surface set(s).");
  // std::cout << "Face numbers :\n";
  // for (auto const & patch: _srfSideIdx)
  //{
  //    std::cout << patch.first << std::endl;
  //    for (auto const & elmSideNum: patch.second)
  //        std::cout << elmSideNum << std::endl;
  //}

  // loop through all FSI element surface faces
  std::vector<size_t> tetFaceConnPtrn = {0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0};
  std::vector<size_t> hexFaceConnPtrn = {0, 1, 2, 3, 4, 7, 6, 5, 0, 4, 5, 1,
                                         1, 5, 6, 2, 2, 6, 7, 3, 3, 7, 4, 0};
  _fsiTriConn.clear();
  _fsiQuadConn.clear();

  for (auto const &elmFIdx : _srfSideIdx["FSI"]) {
    size_t elmIdx = elmFIdx / 10;
    size_t facIdx = elmFIdx - 10 * elmIdx;
    // std::cout << "Elmenet " << elmIdx << " Face " << facIdx << std::endl;

    std::vector<size_t> elmFaceConnTmp;
    // if tetrahedral, this number should be above zero
    if (_lowTetElmNum > 0) {
      for (size_t ndeLoc = 0; ndeLoc < 3; ndeLoc++) {
        size_t faceNdeIdx =
            _tetConn[(elmIdx - _lowTetElmNum) * 4 +
                     tetFaceConnPtrn[(facIdx - 1) * 3 + ndeLoc]];
        elmFaceConnTmp.push_back(faceNdeIdx);
      }

      _fsiTetFaceIdx.insert(std::pair<size_t, size_t>(elmIdx, facIdx));
      _fsiTriConn.insert(_fsiTriConn.end(), elmFaceConnTmp.begin(),
                         elmFaceConnTmp.end());
    }

    // if hexahedral, this number should be above zero
    if (_lowHexElmNum > 0) {
      for (size_t ndeLoc = 0; ndeLoc < 4; ndeLoc++) {
        auto tidx = (elmIdx - _lowHexElmNum) * 8 +
                    hexFaceConnPtrn[(facIdx - 1) * 4 + ndeLoc];
        size_t faceNdeIdx = _hexConn[tidx];
        elmFaceConnTmp.push_back(faceNdeIdx);
      }

      _fsiHexFaceIdx.insert(std::pair<size_t, size_t>(elmIdx, facIdx));
      _fsiQuadConn.insert(_fsiQuadConn.end(), elmFaceConnTmp.begin(),
                          elmFaceConnTmp.end());
    }
  }
  _nFsiTri = _fsiTriConn.size() / 3;
  message("Found " + std::to_string(_nFsiTri) + " Tri facets on FSI surfaces.");
  _nFsiQuad = _fsiQuadConn.size() / 4;
  message("Found " + std::to_string(_nFsiQuad) +
          " Quad facets on FSI surfaces.");
  _nFsiFct = _nFsiTri + _nFsiQuad;

  // coodinates of all fsi nodes
  // renumbering face connectivities
  std::vector<int> _fsiElmConn(_fsiTriConn);
  _fsiElmConn.insert(_fsiElmConn.end(), _fsiQuadConn.begin(),
                     _fsiQuadConn.end());
  _fsiNdeCrd.clear();
  std::map<size_t, size_t> vol2srfConn;
  _nFsiNde = 0;
  for (auto const &indx : _fsiElmConn) {
    // std::cout << " FSI node indx " << indx << std::endl;
    // CGNS is 1 indexed
    auto ret =
        vol2srfConn.insert(std::pair<size_t, size_t>(indx, _nFsiNde + 1));
    if (ret.second == true) {
      _fsiNdeIdx.push_back(indx);
      _nFsiNde++;
      for (size_t iCrd = 0; iCrd < 3; iCrd++) {
        _fsiNdeCrd.push_back(_ci->co[(indx - 1) * 3 + iCrd]);
      }
    }
  }
  message("Found " + std::to_string(_nFsiNde) + " FSI nodes.");
  // std::cout << "FSI node coordinates " << std::endl;
  // for (int iNde = 0; iNde < _fsiNdeCrd.size() / 3; iNde++)
  //  std::cout << "Node " << iNde << " " << _fsiNdeCrd[iNde * 3] << " "
  //            << _fsiNdeCrd[iNde * 3 + 1] << " " << _fsiNdeCrd[iNde * 3 + 2]
  //            << std::endl;
  _fsiTriConnInterface.clear();
  for (auto const &indx : _fsiTriConn)
    _fsiTriConnInterface.push_back(vol2srfConn[indx]);
  ////std::cout << "New connectivities : \n";
  ////for (auto const &idx : _fsiTriConnInterface)
  ////  std::cout << idx << std::endl;

  _fsiQuadConnInterface.clear();
  for (auto const &indx : _fsiQuadConn)
    _fsiQuadConnInterface.push_back(vol2srfConn[indx]);

  message("FSI surface connectivity number of elements is " +
          std::to_string(_fsiTriConnInterface.size() / 3 +
                         _fsiQuadConnInterface.size() / 4));
}

// decomposes a string into tokens using single delimiter
// the implementation allows for arbitrary number of blanks for a blank
// delimitor
std::vector<std::string> clcx_module::tokenize(const char *lineIn,
                                               const char &delim) {
  std::istringstream buf(lineIn);
  if (delim == ' ') {
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end);
    return tokens;
  } else {
    std::string tk;
    std::vector<std::string> tokens;
    while (std::getline(buf, tk, delim))
      tokens.push_back(tk);
    return tokens;
  }
}

void clcx_module::initialize(double &initTime, MPI_Comm &inComm,
                             int &manInitHndl, char *srfWinInChar,
                             char *volWinInChar, int &obtHndl) {
  // since we are multithread only, we limit to process zero
  // communicator here.
  MPI_Comm_rank(inComm, &_rank);
  MPI_Comm_size(inComm, &_nproc);
  MPI_Comm_dup(inComm, &_comm);

  // obtain restart info from volume window and update interface module
  message("Initializing ...");
  message("Time " + std::to_string(initTime));
  message("Surface window " + std::string(srfWinInChar));
  message("Volume window " + std::string(volWinInChar));
  message(std::string("ManInitHandle is ") +
          std::string((manInitHndl < 0) ? ("not set") : ("set")));
  message(std::string("ObtainHandle is ") +
          std::string((obtHndl < 0) ? ("not set") : ("set")));

  // preparations
  warm_up(1);

  std::string volWinIn(volWinInChar);
  std::string srfWinIn(srfWinInChar);

  // sanity checking
  char *names;
  int nItems = 0;
  COM_get_dataitems(volWinIn.c_str(), &nItems, &names);
  std::string dataItemNames(names);
  message("Number of dataitems passed " + std::to_string(nItems));
  if (nItems > 0)
    message("Dataitem names : " + dataItemNames);

  // volume data
  // restoring state from volume window

  if (main_process()) {
    int sz;
    COM_get_size(volWinIn + ".state", 0, &sz);
    // std::cout << "size is " << sz << std::endl;
    char *buf;
    COM_get_array((volWinIn + ".state").c_str(), 0, &buf);
    auto ss = std::make_shared<std::ofstream>("restart_data.bin");
    for (int i = 0; i < sz; i++)
      (*ss) << buf[i];
    ss->close();
    std::ifstream ifs("restart_data.bin");
    _ci->load(ifs);
    _ci->init_files();
  }
  message("Successfully restored the state.");

  // for future debugging
  // std::string data(buf);
  // std::cout << "Data : " << buf << std::endl;
  // std::ofstream fo;
  // fo.open("post.dat");
  // fo << data;
  // fo.close();
  // std::ifstream fi;
  // fi.open("post.dat");
  // ss << fi.rdbuf();
  // fi.close();

  if (main_process()) {
    // preparation operations
    scan_actions();
    fsi_compute_norms();
  }

  // register volume data
  COM_new_window(_wname_vol, _comm);
  register_volume_data();
  message("Successfully restored volume window.");

  // register global data
  register_surface_data();
  message("Successfully restored surface window.");

  // callback to manInitHandle and pass surface names
  if (manInitHndl > 0)
    COM_call_function(manInitHndl, _wname_srf.c_str(), _wname_vol.c_str());
}

void clcx_module::update_solution(double &currTime, double &timeStep,
                                  int &updHndl) {

  if (main_process()) {
    _upd_start_time = currTime;
    _upd_time_step = timeStep;
    _updHndl = &updHndl;
    std::ostringstream out;
    out.precision(4);
    message("Updating solution.");
    out << std::scientific << _ci->get_current_time();
    message("Current time " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << timeStep;
    message("Requested time step " + out.str());
    message("Update_inbuff_handle is " +
            std::string((updHndl < 0) ? ("not set") : ("set")));

    // time step setup
    double clcx_ttime = _ci->get_current_time();
    double &clcx_tinc = _ci->timepar[0];
    double clcx_tper = _ci->timepar[1];
    double clcx_tper_min;
    double clcx_tper_max = _ci->timepar[3];

    // at the begining of the nonlinear steps, the min time step is
    // not scaled by the step size but it will immediately after
    if (clcx_ttime < 100 * std::numeric_limits<double>::epsilon())
      clcx_tper_min = _ci->timepar[2];
    else
      clcx_tper_min = (_ci->timepar[2]) * clcx_tper;

    // testing
    // message(std::to_string(_ci->timepar[0]));
    // message(std::to_string(_ci->timepar[1]));
    // message(std::to_string(_ci->timepar[2]));
    // message(std::to_string(_ci->timepar[3]));

    out.str("");
    out.clear();
    out << std::scientific << clcx_ttime;
    message("Solver current time is " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << clcx_tinc;
    message("User specified time step is " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << clcx_tper_min;
    message("User specified min time step is " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << clcx_tper_max;
    message("User specified max time step is " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << clcx_tper;
    message("User specified final step time " + out.str());

    // sanity checks
    if (timeStep <= 0) {
      message("Negative and zero timestep is not possible.");
      throw;
    }
    if (timeStep < clcx_tper_min) {
      message(
          "Cannot step less than minimum value specified in the input file. "
          "Adjust input file.");
      throw;
    }

    // check if timestep is divisable otherwise change step size
    double check =
        abs(floor(timeStep / clcx_tinc) * 10 - timeStep / clcx_tinc * 10);
    bool divisable = check < 100 * std::numeric_limits<double>::epsilon();
    // message("Check is " + std::to_string(check));
    message("Time step is divisable " + std::to_string(divisable));
    if (!divisable) {
      message("Time step is not divisable, re-computing it.");
      // set to timestep requested by the caller if feasible
      if (timeStep <= clcx_tper_max && timeStep >= clcx_tper_min)
        clcx_tinc = timeStep;
      else if (timeStep > clcx_tper_max) {
        // adjust timestep to exactly meet requested time
        double ratio = timeStep / clcx_tper_max;
        clcx_tinc = timeStep / (1.2 * ratio);
      }
      out.str("");
      out.clear();
      out << std::scientific << clcx_tinc;
      message("Time step was adjusted to " + out.str());
    }

    if (clcx_tinc < clcx_tper_min) {
      message("Cannot step less than the minimum value specified in the input "
              "file. Adjust input file.");
      throw;
    }
  }

  if (main_process()) {
    // synch loads
    fsi_sync_loads();

    // proceed in time
    _ci->step(timeStep);

    // synch other quantities
    fsi_sync_other_quantities();
    // debug_print(_wname_srf + ".u", _rank + 1, 0, _comm, "");
  }
}

void clcx_module::scan_elements() {
  // element labels start from C, D, or E
  for (int iElm = 0; iElm < (_ci->ne); iElm++) {
    if (_ci->lakon[iElm * 8] != 'C' && _ci->lakon[iElm * 8] != 'D' &&
        _ci->lakon[iElm * 8] != 'E')
      continue;
    std::string elmTyp;
    for (int iC = 0; iC < 4; iC++)
      elmTyp.insert(elmTyp.end(), _ci->lakon[iElm * 8 + iC]);
    _elmTypIdx[elmTyp].push_back(iElm);
  }
  message("Highest element number " + std::to_string(_ci->ne));
  message("Element types and numbers ");
  _totElmNum = 0;
  for (auto it = _elmTypIdx.begin(); it != _elmTypIdx.end(); it++) {
    message(it->first + " -> " + std::to_string((it->second).size()));
    _totElmNum += (it->second).size();
  }

  // connectivity
  const std::vector<size_t> &tetIdx = _elmTypIdx["C3D4"];
  _tetConn.clear();
  for (auto it = tetIdx.begin(); it != tetIdx.end(); it++)
    for (int iCon = _ci->ipkon[*it]; iCon <= (_ci->ipkon[*it] + 3); iCon++)
      _tetConn.push_back(_ci->kon[iCon]);

  const std::vector<size_t> &hexIdx = _elmTypIdx["C3D8"];
  _hexConn.clear();
  for (auto it = hexIdx.begin(); it != hexIdx.end(); it++)
    for (int iCon = _ci->ipkon[*it]; iCon <= (_ci->ipkon[*it] + 7); iCon++)
      _hexConn.push_back(_ci->kon[iCon]);

  //// trying to factor out non-tetrahedral elements
  //// two possibilities for this conditional to come true:
  //// 1- there are non-tetrahedral elements in the simulation
  //// 2- the lowest element label is not 1
  //// we ignore other elements if the first conditions holds, but
  //// we need to find lowest tet element label, we assume tets are
  //// defined as the last set of elements
  //_lowTetElmNum = 1;
  // if (tetIdx.size() != (_ci->ne)) {
  //  message("There are non-tetrahedral elements in the input file.");
  //  message("Make sure tetrahedral elements are last set of elements.");
  //  message("Other elements will be ignored.");
  //  _lowTetElmNum = (_ci->ne) - tetIdx.size() + 1;
  //  message("Lowest tetrahedral element label is " +
  //          std::to_string(_lowTetElmNum));
  //}

  // changing the behavior to allow for only tetrahedral and hexahedral
  // elements.
  if (tetIdx.size() != 0 && hexIdx.size() != 0) {
    message("Only tetrahedral and hexahedral elements are supported.");
    throw;
  }

  // removing element index offsets, if any
  // if offset exist it will be hard to deferentiate mixed meshes
  // from offset-only condition. We warn the user about this.
  if (tetIdx.size() > 0) {
    if (tetIdx.size() != (_ci->ne)) {
      message("Mixed meshes are not supported.");
      _lowTetElmNum = (_ci->ne) - tetIdx.size() + 1;
    } else
      _lowTetElmNum = 1;
  }

  if (hexIdx.size() > 0) {
    if (hexIdx.size() != (_ci->ne)) {
      message("Mixed meshes are not supported.");
      _lowHexElmNum = (_ci->ne) - hexIdx.size() + 1;
    } else
      _lowHexElmNum = 1;
  }
}

void clcx_module::fsi_augment_inp() {
  // augment the input file with Surface keyword for
  // fsi surfaces
  std::string jnRoot = std::string(_ci->jobnamec);
  std::ifstream inpf;
  std::ofstream inpfAug;
  inpf.open(jnRoot + ".inp");
  inpfAug.open(jnRoot + "_fsi.inp");
  if (!(inpf.good() || inpfAug.good())) {
    message("Problem openning input file " + jnRoot + ".inp");
    throw;
  }

  // it can be either defining a surface and then defining a dload
  // for the surface. Another possibility, define dload for individual
  // elemenet faces. In the first solution, the face number (x in Px)
  // does not matter (anything btw 1-8). In the second case face numbers
  // should be assigned individually. Another limit, calculix only accepts
  // normal pressure as distributed load.
  // we go with solution 1 and comment out solution 2
  std::string line;
  while (std::getline(inpf, line)) {
    if (line.find("*STEP") != std::string::npos) {
      if (_fsiTetFaceIdx.size() > 0) {
        inpfAug << "*SURFACE, NAME=FSI, TYPE=ELEMENT\n";
        for (auto const &elm : _fsiTetFaceIdx)
          inpfAug << std::to_string(elm.first + 1) << ", S"
                  << std::to_string(elm.second) << std::endl;
      }
    }
    // solution 1, 0.0 arbitrarily selected
    if (line.find("*END STEP") != std::string::npos)
      if (_fsiTetFaceIdx.size() > 0)
        inpfAug << "*DLOAD\nFSI, P1, 0.0\n";
    // solution 2
    // if (line.find("*END STEP") != std::string::npos) {
    //  if (_fsiTetFaceIdx.size() > 0) {
    //    inpfAug << "*DLOAD\n";
    //    for (auto const &elm : _fsiTetFaceIdx)
    //      inpfAug << std::to_string(elm.first + 1) << ", P"
    //              << std::to_string(elm.second) << ", 1.0" << std::endl;
    //  }
    //}
    inpfAug << line << std::endl;
  }

  inpf.close();
  inpfAug.close();
}

void clcx_module::fsi_sync_loads() {

  // current loads
  // for (int iload = 0; iload < (_ci->nload); iload++)
  //  std::cout << "FSI load " << iload << " Element "
  //            << _ci->nelemload[2 * iload] << " load " << _ci->xload[2 *
  //            iload]
  //            << " sideload " << _ci->sideload[20 * iload]
  //            << _ci->sideload[20 * iload + 1] << "\n";

  // sanity checking
  if ((_ci->nload) != (_nFsiFct)) {
    message("Number of FSI facets does not mach with input file!");
    message("Input file count " + std::to_string(_ci->nload) +
            " preprocessor finds " + std::to_string(_nFsiFct));
    throw;
  }

  if (_fsiTs_alp.size() != 3 * (_nFsiFct)) {
    message("The FSI traction vector seems to be not initialized properly. "
            "Presetting with zero load.");
    _fsiTs_alp.resize(3 * (_nFsiFct), 0.0);
  }

  // NOTE: this is just for testing comment out in production
  // for (auto &val : _fsiTs_alp)
  // val = 1.0;

  // computing surface normals and applying pressure loads
  fsi_compute_norms();
  for (int iFct = 0; iFct < _nFsiFct; iFct++) {
    double pressure = _fsiTs_alp[iFct * 3] * _fsiTriNorm[iFct * 3] +
                      _fsiTs_alp[iFct * 3 + 1] * _fsiTriNorm[iFct * 3 + 1] +
                      _fsiTs_alp[iFct * 3 + 2] * _fsiTriNorm[iFct * 3 + 2];
    // std::cout << "Tri " << iTri << " pressure " << pressure << std::endl;
    _ci->xload[iFct * 2] = pressure;
  }

  // current loads
  // for (int iload = 0; iload < (_ci->nload); iload++)
  //  std::cout << "FSI load " << iload << " Element "
  //            << _ci->nelemload[2 * iload] << " load " << _ci->xload[2 *
  //            iload]
  //            << " sideload " << _ci->sideload[20 * iload]
  //            << _ci->sideload[20 * iload + 1] << "\n";
}

void clcx_module::scan_actions() {
  // element type check
  scan_elements();

  // fsi treatments
  scan_fsi();
}

void clcx_module::fsi_sync_other_quantities() {
  // assuming no change in mesh topology
  // so only node coordinates and physical quantities
  // may change

  // coodinates of all fsi nodes
  //_fsiNdeIdx.clear();
  //_fsiNdeCrd.clear();
  // for (auto const &patch : _fsiNSetIdx)
  //  for (auto const &ndeIdx : patch.second) {
  //    _fsiNdeIdx.push_back(ndeIdx);
  //    for (size_t iCrd = 0; iCrd < 3; iCrd++)
  //      _fsiNdeCrd.push_back(_ci->co[(ndeIdx - 1) * 3 + iCrd]);
  //  }
  // std::vector<size_t> const &nIdx = _fsiNSetIdx["FSIN"];
  // std::cout << "FSI coordinates " << std::endl;
  // for (int iNde = 0; iNde < _fsiNdeCrd.size() / 3; iNde++)
  //  std::cout << "Node " << nIdx[iNde] << " " << _fsiNdeCrd[iNde * 3] << " "
  //            << _fsiNdeCrd[iNde * 3 + 1] << " " << _fsiNdeCrd[iNde * 3 + 2]
  //            << std::endl;

  message("Synching quanties at the end of step.");
  // physics quantities
  _fsiU.clear();
  _fsiUHat.resize(_nFsiNde * 3, 0.0);
  _fsiV.clear();
  _fsiCs.resize(_nFsiFct, 0.0);
  _fsiTs.resize(_nFsiNde, 0.0);
  _fsiRhos.resize(_nFsiFct, 0.0);
  _fsiTs_alp.resize(_nFsiFct * 3, 0.0);
  _fsiVbar_alp.resize(_nFsiNde * 3, 0.0);
  _fsiQs_alp.resize(_nFsiFct, 0.0);

  // supporting only displacement and velocity
  // for (auto &ndeIdx : _fsiNSetIdx["C3D4"])
  size_t mt = _ci->mt;
  for (auto &ndeIdx : _fsiNdeIdx)
    for (auto iCmp = 0; iCmp < 3; iCmp++) {
      _fsiU.push_back(_ci->vold[(ndeIdx - 1) * mt + iCmp + 1]);
      _fsiV.push_back(_ci->veold[(ndeIdx - 1) * mt + iCmp + 1]);
    }

  // DBG
  // double min_ci, max_ci;
  // min_ci = _ci->vold[0];
  // max_ci = min_ci;
  // for (size_t idx=0; idx<_ci->nk*(_ci->mt); idx++)
  //{
  //    min_ci = std::min(min_ci, _ci->vold[idx]);
  //    max_ci = std::max(max_ci, _ci->vold[idx]);
  //}
  // std::cout << "mt = " << _ci->mt << std::endl;
  // std::cout << "Max ci = " << max_ci << std::endl;
  // std::cout << "Min ci = " << min_ci << std::endl;

  // updating array pointers
  COM_set_array(_wname_srf + ".u", _rank + 1, &_fsiU[0], 3);
  COM_set_array(_wname_srf + ".vs", _rank + 1, &_fsiV[0], 3);

  // values in _fsiU
  if (_vrb > 0) {
    message("Updated surface displacements");
    double minv, maxv;
    minv = _fsiU[0];
    maxv = minv;
    std::ostringstream out;
    out.precision(4);
    for (auto const &v : _fsiU) {
      out.str("");
      out.clear();
      out << std::scientific << v;
      minv = std::min(minv, v);
      maxv = std::max(maxv, v);
      // message(out.str());
    }
    out.str("");
    out.clear();
    out << std::scientific << maxv;
    message("Max vol.u value is " + out.str());
    out.str("");
    out.clear();
    out << std::scientific << minv;
    message("Min vol.u value is " + out.str());
  }
}

void clcx_module::set_cfd_data(std::string fn) {
  _cfd_data = true;
  _cfd_data_fn = fn;
  message("Trying to use CFD data provided in " + fn);
}

// TODO: supports only tetrahedral meshes
void clcx_module::fsi_sync_loads_from_file() {

  // sanity checking
  if ((_ci->nload) != _nFsiTri) {
    message("Number of FSI triangles does not mach with input file!");
    message("Input file count " + std::to_string(_ci->nload) +
            " preprocessor finds " + std::to_string(_nFsiTri));
    throw;
  }

  if (_fsiTs_alp.size() != 3 * _nFsiTri) {
    message("The FSI traction vector seems to be not initialized properly. "
            "Presetting with zero load.");
    _fsiTs_alp.resize(3 * _nFsiTri, 0.0);
  }

  // if filename is dummy do this dummy loading case
  if (_cfd_data_fn.compare("dummy") == 0) {
    fsi_compute_norms();
    fsi_compute_face_centers();

    for (int iTri = 0; iTri < _nFsiTri; iTri++) {

      double pressure = 0;
      if (_fsiTriCenCrd[iTri * 3] < 0)
        pressure = 1.e5;

      // updating interface quantities
      _fsiTs_alp[iTri * 3] = _fsiTriNorm[iTri * 3] * pressure;
      _fsiTs_alp[iTri * 3 + 1] = _fsiTriNorm[iTri * 3 + 1] * pressure;
      _fsiTs_alp[iTri * 3 + 2] = _fsiTriNorm[iTri * 3 + 2] * pressure;

      // updating calculix data structrue
      _ci->xload[iTri * 2] = pressure;

      // message
      message("Triangle " + std::to_string(iTri) + " pressure " +
              std::to_string(pressure));
    }
    return;
  }

#ifdef HAVE_CFD

  message("Loading FSI pressure data from the CFD file. ");

  // parse and obtain data from the file
  std::ifstream cfd;
  cfd.open(_cfd_data_fn);
  if (!cfd.good()) {
    message("CFD data file is corroupt or not in path.");
    throw;
  }

  std::string line;
  size_t ln = 0;
  std::vector<double> pressure, crds;
  while (std::getline(cfd, line)) {
    ln++;
    // throw away first few lines
    if (ln <= 1)
      continue;

    // get pressure and coords data
    std::vector<std::string> data = tokenize(line.c_str(), ',');
    if (data.size() != 4) {
      message("Corroupt CFD file. Missing data in line " + std::to_string(ln));
      continue;
    }
    pressure.push_back(std::stod(data[0]));
    crds.push_back(std::stod(data[1]));
    crds.push_back(std::stod(data[2]));
    crds.push_back(std::stod(data[3]));
  }
  message("Obtained " + std::to_string(pressure.size()) + " pressure data.");

  // create index
  using namespace flann;
  Matrix<double> dataset(&crds[0], crds.size() / 3, 3);

  // construct an randomized kd-tree index using 4 kd-trees
  Index<L2<double>> index(dataset, KDTreeIndexParams(3));
  index.buildIndex();

  // search index for FSI trianlge facet centers and apply pressures
  fsi_compute_norms();
  fsi_compute_face_centers();
  Matrix<double> query(&_fsiTriCenCrd[0], _nFsiTri, 3);

  // do a knn search, using 128 checks
  Matrix<size_t> indices(new size_t[query.rows], query.rows, 1);
  Matrix<double> dists(new double[query.rows], query.rows, 1);
  index.knnSearch(query, indices, dists, 1, SearchParams(128));
  // flann::save_to_file(indices,"result.hdf5","result");

  for (int iTri = 0; iTri < _nFsiTri; iTri++) {
    message("Triangle " + std::to_string(iTri) + " distance " +
            std::to_string(*dists[iTri]) + " pressure " +
            std::to_string(pressure[*indices[iTri]]));

    // updating interface quantities
    _fsiTs_alp[iTri * 3] = _fsiTriNorm[iTri * 3] * pressure[*indices[iTri]];
    _fsiTs_alp[iTri * 3 + 1] =
        _fsiTriNorm[iTri * 3 + 1] * pressure[*indices[iTri]];
    _fsiTs_alp[iTri * 3 + 2] =
        _fsiTriNorm[iTri * 3 + 2] * pressure[*indices[iTri]];

    // updating calculix data structrue
    _ci->xload[iTri * 2] = pressure[*indices[iTri]];
  }

  // delete[] dataset.ptr();
  // delete[] query.ptr();
  delete[] indices.ptr();
  delete[] dists.ptr();

  cfd.close();
#endif
}

void clcx_module::fsi_compute_norms() {
  // computing surface normals and applying pressure loads
  // for simplicity, we will be using _fsiTriNorm for both
  // tri and quad facets
  _fsiTriNorm.clear();
  for (int iFct = 0; iFct < _nFsiFct; iFct++) {

    std::vector<double> p0, p1, p2, v1, v2, n;

    if (_lowTetElmNum > 0) {
      size_t bId = _fsiTriConn[iFct * 3] - 1;
      p0.push_back(_ci->co[bId * 3]);
      p0.push_back(_ci->co[bId * 3 + 1]);
      p0.push_back(_ci->co[bId * 3 + 2]);

      bId = _fsiTriConn[iFct * 3 + 1] - 1;
      p1.push_back(_ci->co[bId * 3]);
      p1.push_back(_ci->co[bId * 3 + 1]);
      p1.push_back(_ci->co[bId * 3 + 2]);

      bId = _fsiTriConn[iFct * 3 + 2] - 1;
      p2.push_back(_ci->co[bId * 3]);
      p2.push_back(_ci->co[bId * 3 + 1]);
      p2.push_back(_ci->co[bId * 3 + 2]);
    } else if (_lowHexElmNum > 0) {
      size_t bId = _fsiQuadConn[iFct * 4] - 1;
      p0.push_back(_ci->co[bId * 3]);
      p0.push_back(_ci->co[bId * 3 + 1]);
      p0.push_back(_ci->co[bId * 3 + 2]);

      bId = _fsiQuadConn[iFct * 4 + 1] - 1;
      p1.push_back(_ci->co[bId * 3]);
      p1.push_back(_ci->co[bId * 3 + 1]);
      p1.push_back(_ci->co[bId * 3 + 2]);

      bId = _fsiQuadConn[iFct * 4 + 2] - 1;
      p2.push_back(_ci->co[bId * 3]);
      p2.push_back(_ci->co[bId * 3 + 1]);
      p2.push_back(_ci->co[bId * 3 + 2]);

    } else {
      message("Only tetrahedral and hexahedral elements are supported.");
      throw;
    }

    // in calculix, face connectivity is based on normals pointing inwards
    // element switching v1, v2 to compensate for this
    v1.push_back(p1[0] - p0[0]);
    v1.push_back(p1[1] - p0[1]);
    v1.push_back(p1[2] - p0[2]);

    v2.push_back(p2[0] - p0[0]);
    v2.push_back(p2[1] - p0[1]);
    v2.push_back(p2[2] - p0[2]);

    n.push_back(v1[1] * v2[2] - v1[2] * v2[1]);
    n.push_back(v1[0] * v2[2] - v1[2] * v2[0]);
    n.push_back(v1[0] * v2[1] - v1[1] * v2[0]);

    double l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

    if (l < 100 * std::numeric_limits<double>::epsilon()) {
      message("Element face area is close to zero or one of the elements is "
              "degenrate.");
      throw;
    }

    n[0] /= l;
    n[1] /= l;
    n[2] /= l;

    // std::cout << "Face connectivity -> " << _fsiTriConn[iFct * 3] << ", "
    //          << _fsiTriConn[iFct * 3 + 1] << ", " << _fsiTriConn[iFct * 3 +
    //          2]
    //          << std::endl;
    // std::cout << "P0 = " << p0[0] << " " << p0[1] << " " << p0[2] <<
    // std::endl; std::cout << "P1 = " << p1[0] << " " << p1[1] << " " <<
    // p1[2]
    // << std::endl; std::cout << "P2 = " << p2[0] << " " << p2[1] << " " <<
    // p2[2] << std::endl;
    // std::cout << "Tri " << iFct << " Normal " << n[0] << " " << n[1] << " "
    //         << n[2] << std::endl;
    _fsiTriNorm.insert(_fsiTriNorm.end(), n.begin(), n.end());
  }
}

// print double data of an dataitem on a pane
void clcx_module::debug_print(const std::string &str, int pane, int pe,
                              MPI_Comm comm, const char *memo) {
  if (COMMPI_Comm_rank(comm) == pe) {
    double *vm = nullptr;
    int strid, cap;
    printf("%s %s: before %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    COM_get_array(str.c_str(), pane, &vm, &strid, &cap);
    printf("%s %s: after  %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    for (int i = 0; i < strid * cap; i++)
      printf("%.17e ", vm[i]);
    printf("\n");
  }
}

void clcx_module::update_bc(double dt) {
  double alpha = (dt - _upd_start_time) / _upd_time_step;
  message("Requesting FSI BC update for time " + std::to_string(dt));
  message("Update for alpha " + std::to_string(alpha));
  if (_updHndl != nullptr && (*_updHndl > 0)) {
    COM_call_function(*_updHndl, &alpha);
    fsi_sync_loads();
  }
};

void clcx_module::fsi_compute_face_centers() {
  _fsiTriCenCrd.clear();
  for (int iTri = 0; iTri < _nFsiTri; iTri++) {

    std::vector<double> p0, p1, p2;

    size_t bId = _fsiTriConn[iTri * 3] - 1;
    p0.push_back(_ci->co[bId * 3]);
    p0.push_back(_ci->co[bId * 3 + 1]);
    p0.push_back(_ci->co[bId * 3 + 2]);

    bId = _fsiTriConn[iTri * 3 + 1] - 1;
    p1.push_back(_ci->co[bId * 3]);
    p1.push_back(_ci->co[bId * 3 + 1]);
    p1.push_back(_ci->co[bId * 3 + 2]);

    bId = _fsiTriConn[iTri * 3 + 2] - 1;
    p2.push_back(_ci->co[bId * 3]);
    p2.push_back(_ci->co[bId * 3 + 1]);
    p2.push_back(_ci->co[bId * 3 + 2]);

    _fsiTriCenCrd.push_back(1. / 3. * (p0[0] + p1[0] + p2[0]));
    _fsiTriCenCrd.push_back(1. / 3. * (p0[1] + p1[1] + p2[1]));
    _fsiTriCenCrd.push_back(1. / 3. * (p0[2] + p1[2] + p2[2]));
  }
}
