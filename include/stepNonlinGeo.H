#ifndef STEPNONLINGEO_H
#define STEPNONLINGEO_H

#include "CalculiX.h"

class stepNonlinGeo {

  typedef void (*fptr_update_bc)(void *, double);

public:
  stepNonlinGeo() = default;
  ~stepNonlinGeo();

  // public access major methods
public:
  void initialize(
      double **cop, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp, ITG *ne,
      ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, ITG **ipompcp,
      ITG **nodempcp, double **coefmpcp, char **labmpcp, ITG *nmpc,
      ITG *nodeforc, ITG *ndirforc, double *xforc, ITG *nforc, ITG **nelemloadp,
      char **sideloadp, double *xload, ITG *nload, ITG *nactdof, ITG **icolp,
      ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, ITG *nmethod, ITG **ikmpcp,
      ITG **ilmpcp, ITG *ikboun, ITG *ilboun, double *elcon, ITG *nelcon,
      double *rhcon, ITG *nrhcon, double *alcon, ITG *nalcon, double *alzero,
      ITG **ielmatp, ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
      double *t0, double *t1, double *t1old, ITG *ithermal, double *prestr,
      ITG *iprestr, double **voldp, ITG *iperturb, double *sti, ITG *nzs,
      ITG *kode, char *filab, ITG *idrct, ITG *jmax, ITG *jout, double *timepar,
      double *eme, double *xbounold, double *xforcold, double *xloadold,
      double *veold, double *accold, char *amname, double *amta, ITG *namta,
      ITG *nam, ITG *iamforc, ITG **iamloadp, ITG *iamt1, double *alpha,
      ITG *iexpl, ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
      ITG *nplkcon, double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
      char *matname, double *qaold, ITG *mi, ITG *isolver, ITG *ncmat_,
      ITG *nstate_, ITG *iumat, double *cs, ITG *mcs, ITG *nkon, double **enerp,
      ITG *mpcinfo, char *output, double *shcon, ITG *nshcon, double *cocon,
      ITG *ncocon, double *physcon, ITG *nflow, double *ctrl, char *set,
      ITG *nset, ITG *istartset, ITG *iendset, ITG *ialset, ITG *nprint,
      char *prlab, char *prset, ITG *nener, ITG *ikforc, ITG *ilforc,
      double *trab, ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
      ITG *ibody, double *xbody, ITG *nbody, double *xbodyold, ITG *ielprop,
      double *prop, ITG *ntie, char *tieset, ITG *itpamp, ITG *iviewfile,
      char *jobnamec, double *tietol, ITG *nslavs, double *thicke, ITG *ics,
      ITG *nintpoint, ITG *mortar, ITG *ifacecount, char *typeboun,
      ITG **islavsurfp, double **pslavsurfp, double **clearinip, ITG *nmat,
      double *xmodal, ITG *iaxial, ITG *inext, ITG *nprop, ITG *network,
      char *orname, double *vel, ITG *nef, double *velo, double *veloo);

  void increment(
      double **cop, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp, ITG *ne,
      ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, ITG **ipompcp,
      ITG **nodempcp, double **coefmpcp, char **labmpcp, ITG *nmpc,
      ITG *nodeforc, ITG *ndirforc, double *xforc, ITG *nforc, ITG **nelemloadp,
      char **sideloadp, double *xload, ITG *nload, ITG *nactdof, ITG **icolp,
      ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, ITG *nmethod, ITG **ikmpcp,
      ITG **ilmpcp, ITG *ikboun, ITG *ilboun, double *elcon, ITG *nelcon,
      double *rhcon, ITG *nrhcon, double *alcon, ITG *nalcon, double *alzero,
      ITG **ielmatp, ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
      double *t0, double *t1, double *t1old, ITG *ithermal, double *prestr,
      ITG *iprestr, double **voldp, ITG *iperturb, double *sti, ITG *nzs,
      ITG *kode, char *filab, ITG *idrct, ITG *jmax, ITG *jout, double *timepar,
      double *eme, double *xbounold, double *xforcold, double *xloadold,
      double *veold, double *accold, char *amname, double *amta, ITG *namta,
      ITG *nam, ITG *iamforc, ITG **iamloadp, ITG *iamt1, double *alpha,
      ITG *iexpl, ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
      ITG *nplkcon, double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
      char *matname, double *qaold, ITG *mi, ITG *isolver, ITG *ncmat_,
      ITG *nstate_, ITG *iumat, double *cs, ITG *mcs, ITG *nkon, double **enerp,
      ITG *mpcinfo, char *output, double *shcon, ITG *nshcon, double *cocon,
      ITG *ncocon, double *physcon, ITG *nflow, double *ctrl, char *set,
      ITG *nset, ITG *istartset, ITG *iendset, ITG *ialset, ITG *nprint,
      char *prlab, char *prset, ITG *nener, ITG *ikforc, ITG *ilforc,
      double *trab, ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
      ITG *ibody, double *xbody, ITG *nbody, double *xbodyold, ITG *ielprop,
      double *prop, ITG *ntie, char *tieset, ITG *itpamp, ITG *iviewfile,
      char *jobnamec, double *tietol, ITG *nslavs, double *thicke, ITG *ics,
      ITG *nintpoint, ITG *mortar, ITG *ifacecount, char *typeboun,
      ITG **islavsurfp, double **pslavsurfp, double **clearinip, ITG *nmat,
      double *xmodal, ITG *iaxial, ITG *inext, ITG *nprop, ITG *network,
      char *orname, double *vel, ITG *nef, double *velo, double *veloo);

  void finalize(
      double **cop, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp, ITG *ne,
      ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, ITG **ipompcp,
      ITG **nodempcp, double **coefmpcp, char **labmpcp, ITG *nmpc,
      ITG *nodeforc, ITG *ndirforc, double *xforc, ITG *nforc, ITG **nelemloadp,
      char **sideloadp, double *xload, ITG *nload, ITG *nactdof, ITG **icolp,
      ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, ITG *nmethod, ITG **ikmpcp,
      ITG **ilmpcp, ITG *ikboun, ITG *ilboun, double *elcon, ITG *nelcon,
      double *rhcon, ITG *nrhcon, double *alcon, ITG *nalcon, double *alzero,
      ITG **ielmatp, ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
      double *t0, double *t1, double *t1old, ITG *ithermal, double *prestr,
      ITG *iprestr, double **voldp, ITG *iperturb, double *sti, ITG *nzs,
      ITG *kode, char *filab, ITG *idrct, ITG *jmax, ITG *jout, double *timepar,
      double *eme, double *xbounold, double *xforcold, double *xloadold,
      double *veold, double *accold, char *amname, double *amta, ITG *namta,
      ITG *nam, ITG *iamforc, ITG **iamloadp, ITG *iamt1, double *alpha,
      ITG *iexpl, ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
      ITG *nplkcon, double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
      char *matname, double *qaold, ITG *mi, ITG *isolver, ITG *ncmat_,
      ITG *nstate_, ITG *iumat, double *cs, ITG *mcs, ITG *nkon, double **enerp,
      ITG *mpcinfo, char *output, double *shcon, ITG *nshcon, double *cocon,
      ITG *ncocon, double *physcon, ITG *nflow, double *ctrl, char *set,
      ITG *nset, ITG *istartset, ITG *iendset, ITG *ialset, ITG *nprint,
      char *prlab, char *prset, ITG *nener, ITG *ikforc, ITG *ilforc,
      double *trab, ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
      ITG *ibody, double *xbody, ITG *nbody, double *xbodyold, ITG *ielprop,
      double *prop, ITG *ntie, char *tieset, ITG *itpamp, ITG *iviewfile,
      char *jobnamec, double *tietol, ITG *nslavs, double *thicke, ITG *ics,
      ITG *nintpoint, ITG *mortar, ITG *ifacecount, char *typeboun,
      ITG **islavsurfp, double **pslavsurfp, double **clearinip, ITG *nmat,
      double *xmodal, ITG *iaxial, ITG *inext, ITG *nprop, ITG *network,
      char *orname, double *vel, ITG *nef, double *velo, double *veloo);

  int set_step_target_time(const double step_tt);
  bool is_completed();
  double get_total_increment_time() { return (theta * (*tper)); };

public:
  void register_update_bc(fptr_update_bc update_bc, void *callee_context) {
    _update_bc = update_bc;
    _callee_context = callee_context;
  };
  void update_bc(double t) { _update_bc(_callee_context, t); }

  // local management
private:
  ITG clcx_nmethod;
  double theta_t = 0;

  // state variables
private:
  char description[13] = "            ";
  char *lakon = NULL;
  char jobnamef[396] = "";
  char *sideface = NULL;
  char *labmpc = NULL;
  char *lakonf = NULL;
  char *sideloadref = NULL;
  char *sideload = NULL;
  char stiffmatrix[132] = "";

  // state variables
private:
  ITG *inum = NULL;
  ITG k;
  ITG iout = 0;
  ITG icntrl;
  ITG iinc = 0;
  ITG jprint = 0;
  ITG iit = -1;
  ITG jnz = 0;
  ITG icutb = 0;
  ITG istab = 0;
  ITG ifreebody;
  ITG uncoupled;
  ITG n1;
  ITG n2;
  ITG itruecontact;
  ITG iperturb_sav[2];
  ITG ilin;
  ITG *icol = NULL;
  ITG *irow = NULL;
  ITG ielas = 0;
  ITG icmd = 0;
  ITG memmpc_;
  ITG mpcfree;
  ITG icascade;
  ITG maxlenmpc;
  ITG *nodempc = NULL;
  ITG *iaux = NULL;
  ITG *nodempcref = NULL;
  ITG memmpcref_;
  ITG mpcfreeref;
  ITG *itg = NULL;
  ITG *ineighe = NULL;
  ITG *ieg = NULL;
  ITG ntg = 0;
  ITG ntr;
  ITG *kontri = NULL;
  ITG *nloadtr = NULL;
  ITG idamping = 0;
  ITG *ipiv = NULL;
  ITG ntri;
  ITG newstep;
  ITG mode = -1;
  ITG noddiam = -1;
  ITG nasym = 0;
  ITG im;
  ITG ntrit;
  ITG *inocs = NULL;
  ITG inewton = 0;
  ITG *ipobody = NULL;
  ITG *nacteq = NULL;
  ITG *nactdog = NULL;
  ITG nteq;
  ITG *itietri = NULL;
  ITG *koncont = NULL;
  ITG istrainfree = 0;
  ITG ncont;
  ITG ne0;
  ITG nkon0;
  ITG *ipkon = NULL;
  ITG *kon = NULL;
  ITG *ielorien = NULL;
  ITG *ielmat = NULL;
  ITG itp = 0;
  ITG symmetryflag = 0;
  ITG inputformat = 0;
  ITG kscale = 1;
  ITG *iruc = NULL;
  ITG iitterm = 0;
  ITG iturbulent;
  ITG ngraph = 1;
  ITG ismallsliding = 0;
  ITG *ipompc = NULL;
  ITG *ikmpc = NULL;
  ITG *ilmpc = NULL;
  ITG i0ref;
  ITG irref;
  ITG icref;
  ITG *itiefac = NULL;
  ITG *islavsurf = NULL;
  ITG *islavnode = NULL;
  ITG *imastnode = NULL;
  ITG *nslavnode = NULL;
  ITG *nmastnode = NULL;
  ITG *imastop = NULL;
  ITG imat;
  ITG *iponoels = NULL;
  ITG *inoels = NULL;
  ITG *islavsurfold = NULL;
  ITG maxlenmpcref;
  ITG *islavact = NULL;
  ITG mt;
  ITG *nactdofinv = NULL;
  ITG *ipe = NULL;
  ITG *ime = NULL;
  ITG *ikactmech = NULL;
  ITG nactmech;
  ITG inode;
  ITG idir;
  ITG neold;
  ITG neini;
  ITG iemchange = 0;
  ITG nzsrad;
  ITG *mast1rad = NULL;
  ITG *irowrad = NULL;
  ITG *icolrad = NULL;
  ITG *jqrad = NULL;
  ITG *ipointerrad = NULL;
  ITG *integerglob = NULL;
  ITG negpres = 0;
  ITG mass[2] = {0, 0};
  ITG stiffness = 1;
  ITG buckling = 0;
  ITG rhsi = 1;
  ITG intscheme = 0;
  ITG idiscon = 0;
  ITG coriolis = 0;
  ITG *ipneigh = NULL;
  ITG *neigh = NULL;
  ITG maxprevcontel;
  ITG nslavs_prev_step;
  ITG *nelemface = NULL;
  ITG *ipoface = NULL;
  ITG *nodface = NULL;
  ITG *ifreestream = NULL;
  ITG iex;
  ITG *isolidsurf = NULL;
  ITG *neighsolidsurf = NULL;
  ITG *iponoel = NULL;
  ITG *inoel = NULL;
  ITG nface;
  ITG nfreestream;
  ITG nsolidsurf;
  ITG i;
  ITG icfd = 0;
  ITG id;
  ITG *neij = NULL;
  ITG node;
  ITG networknode;
  ITG iflagact = 0;
  ITG *nodorig = NULL;
  ITG *ipivr = NULL;
  ITG iglob = 0;
  ITG *inomat = NULL;
  ITG *ipnei = NULL;
  ITG ntrimax;
  ITG *nx = NULL;
  ITG *ny = NULL;
  ITG *nz = NULL;
  ITG *neifa = NULL;
  ITG *neiel = NULL;
  ITG *ielfa = NULL;
  ITG *ifaext = NULL;
  ITG nflnei;
  ITG nfaext;
  ITG idampingwithoutcontact = 0;
  ITG *nactdoh = NULL;
  ITG *nactdohinv = NULL;
  ITG *ipkonf = NULL;
  ITG *ielmatf = NULL;
  ITG *ielorienf = NULL;
  ITG ialeatoric = 0;
  ITG nloadref;
  ITG isym;
  ITG *nelemloadref = NULL;
  ITG *iamloadref = NULL;
  ITG *idefload = NULL;
  ITG nload_;
  ITG *nelemload = NULL;
  ITG *iamload = NULL;
  ITG ncontacts = 0;
  ITG inccontact = 0;
  ITG nrhs = 1;
  ITG j = 0;
  ITG *ifatie = NULL;
  ITG n;
  ITG inoelsize = 0;
  ITG isensitivity = 0;
  ITG *istartblk = NULL;
  ITG *iendblk = NULL;
  ITG *nblket = NULL;
  ITG *nblkze = NULL;
  ITG nblk;
  ITG *konf = NULL;
  ITG *ielblk = NULL;
  ITG *iwork = NULL;
  ITG nelt;
  ITG lrgw;
  ITG *igwk = NULL;
  ITG itol;
  ITG itmax;
  ITG iter;
  ITG ierr;
  ITG iunit;
  ITG ligw;
  ITG mei[4] = {0, 0, 0, 0};
  ITG initial;

  // state variables
private:
  double *stn = NULL;
  double *v = NULL;
  double *een = NULL;
  double cam[5];
  double *epn = NULL;
  double *cg = NULL;
  double *cdn = NULL;
  double *vfa = NULL;
  double *pslavsurfold = NULL;
  double *f = NULL;
  double *fn = NULL;
  double qa[4] = {0., 0., -1., 0.};
  double qam[2] = {0., 0.};
  double dtheta;
  double theta;
  double err;
  double ram[6] = {0., 0., 0., 0., 0., 0.};
  double *areaslav = NULL;
  double *springarea = NULL;
  double ram1[6] = {0., 0., 0., 0., 0., 0.};
  double ram2[6] = {0., 0., 0., 0., 0., 0.};
  double deltmx;
  double ptime;
  double smaxls;
  double sminls;
  double uam[2] = {0., 0.};
  double *vini = NULL;
  double *ac = NULL;
  double qa0;
  double qau;
  double ea;
  double *straight = NULL;
  double *t1act = NULL;
  double qamold[2];
  double *xbounact = NULL;
  double *bc = NULL;
  double *xforcact = NULL;
  double *xloadact = NULL;
  double *fext = NULL;
  double *clearini = NULL;
  double reltime;
  double time;
  double bet = 0.;
  double gam = 0.;
  double *aux2 = NULL;
  double dtime;
  double *fini = NULL;
  double *fextini = NULL;
  double *veini = NULL;
  double *accini = NULL;
  double *xstateini = NULL;
  double *ampli = NULL;
  double scal1;
  double *eei = NULL;
  double *t1ini = NULL;
  double pressureratio;
  double *xbounini = NULL;
  double dev;
  double *xstiff = NULL;
  double *stx = NULL;
  double *stiini = NULL;
  double *enern = NULL;
  double *coefmpc = NULL;
  double *aux = NULL;
  double *xstaten = NULL;
  double *coefmpcref = NULL;
  double *enerini = NULL;
  double *emn = NULL;
  double alpham;
  double betam;
  double *tarea = NULL;
  double *tenv = NULL;
  double *erad = NULL;
  double *fnr = NULL;
  double *fni = NULL;
  double *adview = NULL;
  double *auview = NULL;
  double *qfx = NULL;
  double *cvini = NULL;
  double *cv = NULL;
  double *qfn = NULL;
  double *co = NULL;
  double *vold = NULL;
  double *fenv = NULL;
  double sigma = 0.;
  double *xbodyact = NULL;
  double *cgr = NULL;
  double dthetaref;
  double *vr = NULL;
  double *vi = NULL;
  double *stnr = NULL;
  double *stni = NULL;
  double *vmax = NULL;
  double *stnmax = NULL;
  double *fmpc = NULL;
  double *ener = NULL;
  double *f_cm = NULL;
  double *f_cs = NULL;
  double *adc = NULL;
  double *auc = NULL;
  double *res = NULL;
  double *xstate = NULL;
  double *eenmax = NULL;
  double *adrad = NULL;
  double *aurad = NULL;
  double *bcr = NULL;
  double *xmastnor = NULL;
  double *emeini = NULL;
  double *tinc;
  double *tper;
  double *tmin;
  double *tmax;
  double *tincf;
  double *doubleglob = NULL;
  double *xnoels = NULL;
  double *au = NULL;
  double *resold = NULL;
  double *ad = NULL;
  double *b = NULL;
  double *aub = NULL;
  double *adb = NULL;
  double *pslavsurf = NULL;
  double *pmastsurf = NULL;
  double *x = NULL;
  double *y = NULL;
  double *z = NULL;
  double *xo = NULL;
  double sum1;
  double sum2;
  double flinesearch;
  double *yo = NULL;
  double *zo = NULL;
  double *cdnr = NULL;
  double *cdni = NULL;
  double *fnext = NULL;
  double *fnextini = NULL;
  double allwk = 0.;
  double allwkini;
  double energy[4] = {0., 0., 0., 0.};
  double energyini[4];
  double energyref;
  double denergymax;
  double dtcont;
  double dtvol;
  double *wavespeed;
  double emax;
  double r_abs;
  double enetoll;
  double dampwk = 0.;
  double dampwkini = 0.;
  double temax;
  double *tmp = NULL;
  double energystartstep[4];
  double sizemaxinc;
  double *adblump = NULL;
  double *adcpy = NULL;
  double *aucpy = NULL;
  double *rwork = NULL;
  double *sol = NULL;
  double *rgwk = NULL;
  double tol;
  double *sb = NULL;
  double *sx = NULL;
  double delcon;
  double alea;

  // state variable
private:
  FILE *f1;
  fptr_update_bc _update_bc;
  void *_callee_context;
};

#endif
