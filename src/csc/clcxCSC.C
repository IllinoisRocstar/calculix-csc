#include "clcxCSC.H"
#include "stepNonlinGeo.H"
#include <algorithm>
#include <iostream>
#include <string.h>
#include <sstream>

// defaults constructor
clcx_module::clcx_module()
    : _vrb(0), _rank(0), _nproc(1)
{
    _ci = std::make_shared<clcx_interface>();      
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
  COM_Type types[5];

  types[0] = COM_RAWDATA;
  types[1] = COM_INT;
  COM_set_member_function((name + ".initialize").c_str(),
                          (Member_func_ptr)(&clcx_module::initialize),
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

  if (_rank == 0 && (~enfv || (_vrb > 0 && enfv)))
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

void clcx_module::warm_up(int vrbin)
{
  _vrb = vrbin;
  _ci->set_verb(vrbin);
  _ci->myid = _rank;
  _ci->nproc = _nproc;
  
  // window names for internal use
  _wname_vol = _wname + "_vol";
  _wname_srf = _wname + "_srf";
}

void clcx_module::initialize(int vrb) {
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
  _ci->finalize();
  message("Done finalizing");
}

void clcx_module::step(const double &stp_tt) {
  _ci->step(stp_tt);
}

void clcx_module::preprocess(std::string actStr)
{
    warm_up();
    message("Pre-processing....");

    // read input file and generate respective data structure
    _ci->init();

    // register global data
    register_global_data();
}


void clcx_module::register_global_data()
{
    // volume diata
    // pane 1
    COM_new_window(_wname_vol, MPI_COMM_NULL);
    COM_window_init_done(_wname_vol);
    COM_set_size(_wname_vol+".nc", 1, _ci->nk);
    COM_set_array(_wname_vol+".nc", 1, _ci->co, 3);
    // assuming all hex
    std::cout << "ne = " << _ci->ne << std::endl;
    COM_set_size(_wname_vol+".:H8", 1, _ci->ne);
    COM_set_array(_wname_vol+".:H8", 1, _ci->kon, _ci->ne*8);
}
