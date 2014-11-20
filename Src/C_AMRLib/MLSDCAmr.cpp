/*
 * Multilevel SDC + AMR controller.
 *
 * When RNS is compiled with USE_SDCLIB=TRUE, the time-stepping is
 * done using the multi-level SDC (MLSDC) algorithm, with IMEX
 * sweepers on each level (hydrodynamics are explicit, chemistry is
 * implicit).  The MLSDC algorithm is implemented in C in SDCLib (in
 * SDCLib multi-level SDC is called multi-grid SDC).
 *
 * The interface between SDCLib and RNS is (mostly) contained in
 * SDCAmr (derived from Amr) and SDCAmrEncap.
 *
 * Note that in MLSDC, there is no concept of "sub-cycling" fine
 * levels.  As such, the "timeStep" method defined in SDCAmr is only
 * called on the coarsest level, and it essentially passes control to
 * SDCLib to advance the solution.
 *
 * SDCLib handles interpolation and restriction in time.  SDCLib also
 * takes care of allocating multifabs at each SDC node.  As such the
 * traditional Amr concepts of "old data", "new data", and defining an
 * "advance" routine no longer apply.  In order to reuse as much of
 * the existing Amr code as possible, SDCAmr puts the solution
 * obtained from SDCLib into the "new data" state (this means all the
 * logic to write plotfiles etc still works), but never defines "old
 * data" or any "mid data".
 *
 * Some notes:
 *
 * 1. Since we never use 'old data', any calls to fill patch iterator
 *    will always use 'new data', regardless of any 'time' information
 *    set on the state data.
 *
 * 2. The 'n_factor' and 'n_cycle' logic in computeNewDt dictates
 *    that: if n_cycles is set to 1, 2, 2, ..., then dt_level stores
 *    the cfl timestep that each level wants to take.
 *
 * Known issues:
 *
 * 1. The SDC hierarchy is currently not depth limited.
 *
 * 2. We're using Gauss-Lobatto nodes, but Gauss-Radau would probably
 *    be better for chemistry.
 *
 * 3. mlsdc_amr_interpolate won't work for correcton at wall boundary.
 */

#include <MLSDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>
#include <Interpolater.H>
#include <FabArray.H>
#include <iomanip>

#ifdef USE_COLOROUTPUT
// only works on some systems
#define RESETCOLOR       "\033[0m"
#define BOLDFONT         "\033[1m"
#define REDCOLOR         "\033[31m"      /* Red */
#define GREENCOLOR       "\033[32m"      /* Green */
#define YELLOWCOLOR      "\033[33m"      /* Yellow */
#define BLUECOLOR        "\033[34m"      /* Blue */
#define MAGENTACOLOR     "\033[35m"      /* Magenta */
#define CYANCOLOR        "\033[36m"      /* Cyan */
#else
#define RESETCOLOR       ""
#define BOLDFONT         ""
#define REDCOLOR         ""
#define GREENCOLOR       ""
#define YELLOWCOLOR      ""
#define BLUECOLOR        ""
#define MAGENTACOLOR     ""
#define CYANCOLOR        ""
#endif

using namespace std;

BEGIN_EXTERN_C

/*
 * Spatial interpolation between MultiFabs.  Called by SDCLib.
 */
void mlsdc_amr_interpolate(void *Fp, void *Gp, sdc_state *state, void *ctxF, void *ctxG)
{
  BL_PROFILE("MLSDC_AMR_INTERPOLATE()");

  bool correction = state->flags & SDC_CORRECTION;
  bool evaluation = state->flags & SDC_FEVAL;

  MLSDCAmrEncap& F = *((MLSDCAmrEncap*) Fp);
  MLSDCAmrEncap& G = *((MLSDCAmrEncap*) Gp);
  MLSDCAmrLevel& levelF = *((MLSDCAmrLevel*) ctxF);
  MLSDCAmrLevel& levelG = *((MLSDCAmrLevel*) ctxG);

  levelF.interpolate(F, G, state->t, correction, evaluation, levelG);
}

/*
 * Spatial restriction between solutions and integrals of function
 * evals.  Called by SDCLib.
 */
void mlsdc_amr_restrict(void *Fp, void *Gp, sdc_state *state, void *ctxF, void *ctxG)
{
  BL_PROFILE("mlsdc_amr_restrict()");

  MLSDCAmrEncap& F = *((MLSDCAmrEncap*) Fp);
  MLSDCAmrEncap& G = *((MLSDCAmrEncap*) Gp);
  MLSDCAmrLevel& levelF = *((MLSDCAmrLevel*) ctxF);
  MLSDCAmrLevel& levelG = *((MLSDCAmrLevel*) ctxG);

  BL_ASSERT(G.type==SDC_SOLUTION || G.type==SDC_TAU);

  levelF.restrict(F, G, state->t, levelG);
}

END_EXTERN_C


void MLSDCAmr::FinalIntegrate(double t, double dt, int niter)
{
  int nlevs = mg.nlevels;
  for (int l=0; l<nlevs; l++) {
    sdc_sweeper* swp    = mg.sweepers[l];
    swp->sweep(swp, t, dt, niter, SDC_SWEEP_NOEVAL);
  }
}

/*
 * Rebuild MLSDC hierarchy.
 */
void MLSDCAmr::RebuildMLSDC()
{
  BL_PROFILE("MLSDCAmr::RebuildMLSDC()");

  // reset previous and clear sweepers etc
  sdc_mg_reset(&mg);
  DestroyMLSDC();

  // rebuild
  for (int lev=0; lev<=finest_level; lev++) {
    encaps[lev] = BuildEncap(lev);
    sweepers[lev] = BuildLevel(lev);
    sweepers[lev]->nset->ctx = &getLevel(lev);
    sweepers[lev]->nset->encap = encaps[lev];
    sdc_mg_add_level(&mg, sweepers[lev], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  }

  sdc_mg_setup(&mg, 0);
  sdc_mg_allocate(&mg);

  n_cycle[0] = 1;
  for (int lev=1; lev<=finest_level; lev++)
    n_cycle[lev] = 2; // see notes

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    cout << "Rebuilt MLSDC: " << mg.nlevels << ", nnodes: ";
    for (int lev=0; lev<=finest_level; lev++)
      cout << sweepers[lev]->nset->nnodes << " ";
    cout << endl;
    if (verbose > 2)
      sdc_mg_print(&mg, 2);
  }
}

/*
 * Initialize SDC multilevel sweeper, set parameters.
 */
MLSDCAmr::MLSDCAmr ()
{
  InitializeMLSDC();
}

void MLSDCAmr::InitializeMLSDC ()
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 4;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 2;
  if (!ppsdc.query("nnodes0",   nnodes0))   nnodes0 = 3;
  if (!ppsdc.query("trat",      trat))      trat = 2;

  if (verbose > 2)
    sdc_log_set_stdout(SDC_LOG_DEBUG);
  else if (verbose > 1)
    sdc_log_set_stdout(SDC_LOG_INFO);

  sdc_mg_build(&mg, max_level+1);

  sweepers.resize(max_level+1);
  encaps.resize(max_level+1);
  nsweeps.resize(max_trefs+1);

  if (!ppsdc.queryarr("nsweeps", nsweeps)) {
    nsweeps[0] = (max_level > 0) ? 2 : 1;
    for (int l=1; l<max_trefs+1; l++)
      nsweeps[l] = 1;
  }

  for (int i=0; i<=max_level; i++)
    sweepers[i] = NULL;
}

/*
 * Destroy SDC multigrid sweeper, set parameters.
 */

MLSDCAmr::~MLSDCAmr()
{
  DestroyMLSDC();
  sdc_mg_destroy(&mg);
}

void MLSDCAmr::DestroyMLSDC()
{
  for (unsigned int lev=0; lev<=max_level; lev++) {
    if (sweepers[lev] != NULL) {
      sweepers[lev]->destroy(sweepers[lev]);
      sweepers[lev] = NULL;
      delete (MLSDCAmrEncapCtx*) encaps[lev]->ctx;
      delete encaps[lev];
    }
  }
}
