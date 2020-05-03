#ifndef CLCXINTERFACE_H
#define CLCXINTERFACE_H

#include <string>
#include <sstream>
#include <map>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>


// Calculix
#include "CalculiX.h"


class stepNonlinGeo;

class clcx_interface {


public:
  clcx_interface();
  ~clcx_interface(){};

  // helper methods
public:
  // initializes the Calculix interface by reading an input file
  // and populating native data structure
  void init();

  // perform post-initialization steps (if multiple steps needed should be
  // called)
  void postInit();

  // runs the simulation from A-Z without any external control
  void run();

  // finalizes the interface and deletes all memory references
  void finalize();

  // steps the simulation in time to reach step_tt (target time)
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

  // helper methods
public:
  // writes out message on processs rank 0
  int message(std::string msg);
  void splash();
  void step_check();
  void step_initialize();
  void step_run(const double step_tt);
  bool step_completed();
  void step_finalize();
  void set_verb(const int vrb) { _vrb = vrb; };

  // internal management
private:
  int _vrb;
  bool _step_initialized;
  bool _step_finalized;
  bool _is_dynamic;
  bool _is_static;
  bool _is_thermal;
  stepNonlinGeo *nlgeo;

  // storing and restore
private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version);
public:
  void save(std::stringstream& ss);
  void load(std::stringstream& ss);

  // sharable data used by calculix
  // TODO: most of these should be turned to private
public:
  ITG myid = 0;
  ITG nproc = 0;
  FILE *f1;

  // raw data used by calculix
public:
  char *sideload = NULL;
  char *set = NULL;
  char *matname = NULL;
  char *orname = NULL;
  char *amname = NULL;
  char *filab = NULL;
  char *lakon = NULL;
  char *labmpc = NULL;
  char *prlab = NULL;
  char *prset = NULL;
  char jobnamec[660] = "";
  char jobnamef[132] = "";
  char output[4] = "exo";
  char *typeboun = NULL;
  char *inpc = NULL;
  char *tieset = NULL;
  char *cbody = NULL;
  char fneig[132] = "";
  char *sideloadtemp = NULL;
  char kind1[2] = "T";
  char kind2[2] = "T";
  char *heading = NULL;
  char *objectset = NULL;

  // raw data used by calculix
public:
  ITG *kon = NULL;
  ITG *nodeboun = NULL;
  ITG *ndirboun = NULL;
  ITG *ipompc = NULL;
  ITG *nodempc = NULL;
  ITG *nodeforc = NULL;
  ITG *ndirforc = NULL;
  ITG *nelemload = NULL;
  ITG im;
  ITG *inodesd = NULL;
  ITG nload1;
  ITG *idefforc = NULL;
  ITG *nactdof = NULL;
  ITG *icol = NULL;
  ITG *ics = NULL;
  ITG *jq = NULL;
  ITG *mast1 = NULL;
  ITG *irow = NULL;
  ITG *rig = NULL;
  ITG *idefbody = NULL;
  ITG *ikmpc = NULL;
  ITG *ilmpc = NULL;
  ITG *ikboun = NULL;
  ITG *ilboun = NULL;
  ITG *nreorder = NULL;
  ITG *ipointer = NULL;
  ITG *idefload = NULL;
  ITG *istartset = NULL;
  ITG *iendset = NULL;
  ITG *ialset = NULL;
  ITG *ielmat = NULL;
  ITG *ielorien = NULL;
  ITG *nrhcon = NULL;
  ITG *nodebounold = NULL;
  ITG *ndirbounold = NULL;
  ITG *nelcon = NULL;
  ITG *nalcon = NULL;
  ITG *iamforc = NULL;
  ITG *iamload = NULL;
  ITG *iamt1 = NULL;
  ITG *namta = NULL;
  ITG *ipkon = NULL;
  ITG *iamboun = NULL;
  ITG *nplicon = NULL;
  ITG *nplkcon = NULL;
  ITG *inotr = NULL;
  ITG *iponor = NULL;
  ITG *knor = NULL;
  ITG *ikforc = NULL;
  ITG *ilforc = NULL;
  ITG *iponoel = NULL;
  ITG *inoel = NULL;
  ITG *nshcon = NULL;
  ITG *ncocon = NULL;
  ITG *ibody = NULL;
  ITG *ielprop = NULL;
  ITG *islavsurf = NULL;
  ITG *ipoinpc = NULL;
  ITG mt;
  ITG nxstate;
  ITG nload0;
  ITG iload;
  ITG *iuel = NULL;

  // raw data used by calculix
public:
  ITG nk;
  ITG ne;
  ITG nboun;
  ITG nmpc;
  ITG nforc;
  ITG nload;
  ITG nprint = 0;
  ITG nset;
  ITG nalset;
  ITG nentries = 17;
  ITG nmethod;
  ITG neq[3] = {0, 0, 0};
  ITG i;
  ITG mpcfree = 1;
  ITG mei[4];
  ITG j;
  ITG nzl;
  ITG nam;
  ITG nbounold = 0;
  ITG nforcold = 0;
  ITG nloadold = 0;
  ITG nbody;
  ITG nbody_ = 0;
  ITG nbodyold = 0;
  ITG network = 0;
  ITG nheading_ = 0;
  ITG k;
  ITG nzs[3];
  ITG nmpc_ = 0;
  ITG nload_ = 0;
  ITG nforc_ = 0;
  ITG istep;
  ITG istat;
  ITG nboun_ = 0;
  ITG nintpoint = 0;
  ITG iperturb[2] = {0, 0};
  ITG nmat;
  ITG ntmat_ = 0;
  ITG norien;
  ITG ithermal[2] = {0, 0};
  ITG nmpcold;
  ITG iprestr;
  ITG kode;
  ITG isolver = 0;
  ITG nslavs = 0;
  ITG nkon_ = 0;
  ITG ne0;
  ITG nkon0;
  ITG mortar = 0;
  ITG jout[2] = {1, 1};
  ITG nlabel;
  ITG nkon = 0;
  ITG idrct;
  ITG jmax[2];
  ITG iexpl;
  ITG nevtot = 0;
  ITG ifacecount = 0;
  ITG iplas = 0;
  ITG npmat_ = 0;
  ITG mi[3] = {0, 3, 1};
  ITG ntrans;
  ITG mpcend = -1;
  ITG namtot_ = 0;
  ITG iumat = 0;
  ITG icascade = 0;
  ITG maxlenmpc;
  ITG mpcinfo[4];
  ITG ne1d = 0;
  ITG ne2d = 0;
  ITG infree[4] = {0, 0, 0, 0};
  ITG callfrommain;
  ITG nflow = 0;
  ITG jin = 0;
  ITG irstrt[2] = {0, 0};
  ITG nener = 0;
  ITG jrstrt = 0;
  ITG nenerold;
  ITG nline;
  ITG *ipoinp = NULL;
  ITG *inp = NULL;
  ITG ntie;
  ITG ntie_ = 0;
  ITG mcs = 0;
  ITG nprop_ = 0;
  ITG nprop = 0;
  ITG itpamp = 0;
  ITG iviewfile;
  ITG nkold;
  ITG nevdamp_ = 0;
  ITG npt_ = 0;
  ITG cyclicsymmetry;
  ITG nmethodl;
  ITG iaxial = 1;
  ITG inext = 0;
  ITG icontact = 0;
  ITG nobject = 0;
  ITG nobject_ = 0;
  ITG iit = -1;
  ITG nzsprevstep[3];
  ITG memmpcref_;
  ITG mpcfreeref = -1;
  ITG maxlenmpcref;
  ITG *nodempcref = NULL;
  ITG *ikmpcref = NULL;
  ITG isens = 0;
  ITG namtot = 0;
  ITG nstam = 0;
  ITG ndamp = 0;
  ITG nef = 0;

  // raw data used by calculix
public:
  ITG *meminset = NULL;
  ITG *rmeminset = NULL;

  // raw data used by calculix
public:
  ITG nzs_;
  ITG nk_ = 0;
  ITG ne_ = 0;
  ITG nset_ = 0;
  ITG nalset_ = 0;
  ITG nmat_ = 0;
  ITG norien_ = 0;
  ITG nam_ = 0;
  ITG ntrans_ = 0;
  ITG ncs_ = 0;
  ITG nstate_ = 0;
  ITG ncmat_ = 0;
  ITG memmpc_ = 0;
  ITG nprint_ = 0;
  ITG nuel_ = 0;

  // raw data used by calculix
public:
  double *co = NULL;
  double *xboun = NULL;
  double *coefmpc = NULL;
  double *xforc = NULL;
  double *clearini = NULL;
  double *xload = NULL;
  double *xbounold = NULL;
  double *xforcold = NULL;
  double *vold = NULL;
  double *sti = NULL;
  double *xloadold = NULL;
  double *xnor = NULL;
  double *reorder = NULL;
  double *dcs = NULL;
  double *thickn = NULL;
  double *thicke = NULL;
  double *offset = NULL;
  double *elcon = NULL;
  double *rhcon = NULL;
  double *alcon = NULL;
  double *alzero = NULL;
  double *t0 = NULL;
  double *t1 = NULL;
  double *prestr = NULL;
  double *orab = NULL;
  double *amta = NULL;
  double *veold = NULL;
  double *accold = NULL;
  double *t1old = NULL;
  double *eme = NULL;
  double *plicon = NULL;
  double *pslavsurf = NULL;
  double *plkcon = NULL;
  double *xstate = NULL;
  double *trab = NULL;
  double *ener = NULL;
  double *shcon = NULL;
  double *cocon = NULL;
  double *cs = NULL;
  double *tietol = NULL;
  double *fmpc = NULL;
  double *prop = NULL;
  double *t0g = NULL;
  double *t1g = NULL;
  double *xbody = NULL;
  double *xbodyold = NULL;
  double *coefmpcref = NULL;
  double *dacon = NULL;
  double *vel = NULL;
  double *velo = NULL;
  double *veloo = NULL;

  // raw data used by calculix
public:
  double ctrl[56] = {4.5,   8.5,   9.5,   16.5,  10.5,  4.5,   0.,    5.5,
                     0.,    0.,    0.25,  0.5,   0.75,  0.85,  0.,    0.,
                     1.5,   0.,    0.005, 0.01,  0.,    0.,    0.02,  1.e-5,
                     1.e-3, 1.e-8, 1.e30, 1.5,   0.25,  1.01,  1.,    1.,
                     5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, -1.,
                     1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.5,
                     0.5,   20.5,  1.5,   1.5,   0.001, 0.1,   100.5, 60.5};

  // raw data used by calculix
public:
  double fei[3];
  double *xmodal = NULL;
  double timepar[5];
  double alpha;
  double ttime = 0.;
  double qaold[2] = {0., 0.};
  double physcon[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

#ifdef BAM
public:
  ITG lop = 0;
  ITG lrestart = 0;
  ITG kstep = 1;
  ITG kinc = 1;
  double time[2];
  double dtime;
#endif
};

#endif // CLCXINTERFACE_H
