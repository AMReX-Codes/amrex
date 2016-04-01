
#ifndef _MLSDCAMR_H_
#define _MLSDCAMR_H_

#include <MultiFab.H>
#include <FluxRegister.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <sdc.h>
#include <BL_CXX11.H>

using namespace std;

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C }


/*
 * Multilevel SDC+AMR level class.
 */

class MLSDCAmr;
class MLSDCAmrLevel;

/*
 * MultiFab encapsulation context for SDCLib.
 */

typedef struct {
  const BoxArray *ba;
  IntVect crse_ratio;
  int ncomp, ngrow, level, finest;
  MLSDCAmrLevel *amrlevel;
} MLSDCAmrEncapCtx;

typedef struct {
  MultiFab* U;			// solution or function eval
  FluxRegister* fine_flux;	// current levels flux register
  FluxRegister* crse_flux;	// coarse levels flux register
  MLSDCAmrLevel* amrlevel;
  int type;
} MLSDCAmrEncap;

/*
 * Multilevel SDC+AMR level class.
 */

class MLSDCAmrLevel : public AmrLevel
{
public:

  MLSDCAmrLevel () : AmrLevel() { };

  MLSDCAmrLevel (Amr&            papa,
                 int             lev,
                 const Geometry& level_geom,
                 const BoxArray& bl,
                 Real            time) : AmrLevel(papa, lev, level_geom, bl, time) { };

  MLSDCAmr* getMLSDCAmr()
  {
    return (MLSDCAmr*) parent;
  }

  virtual void restrict(MLSDCAmrEncap& F, MLSDCAmrEncap& G, Real t, MLSDCAmrLevel& levelG) = 0;
  virtual void interpolate(MLSDCAmrEncap& F, MLSDCAmrEncap& G, Real t, bool correction, bool evaluation, MLSDCAmrLevel& levelG) = 0;
};

/*
 * Multilevel SDC+AMR class.
 */

class MLSDCAmr : public Amr
{
public:

  int max_iters;                // maximum number of sdc iterations
  int max_trefs;                // maximum number of time refinements
  int nnodes0;			// number of nodes on the coarsest level
  int trat;			// time refinement ratio
  int first_refinement_level;   // first amr level that has time refinement

  sdc_mg mg;
  vector<sdc_sweeper*> sweepers;
  vector<sdc_encap*>   encaps;
  vector<int>          nsweeps;

  MLSDCAmr();
  ~MLSDCAmr();

  virtual sdc_encap*   BuildEncap(int lev);
  virtual sdc_sweeper* BuildLevel(int lev) = 0;
  virtual void         InitializeMLSDC();
  virtual void         RebuildMLSDC();
  virtual void         DestroyMLSDC();
  virtual void         FinalIntegrate(double t, double dt, int niter);

  virtual void coarseTimeStep (Real stop_time) = 0;

  virtual void timeStep(int level,
			Real time,
			int iteration,
			int niter ,
			Real stop_time) BL_OVERRIDE 
	{ BL_ASSERT(false); }

};

#endif
