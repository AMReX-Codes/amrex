
#include <MLSDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>

#include <cassert>

#include <stdio.h>

using namespace std;


void MLSDCAmr::timeStep (int  level,
                         Real time,
                         int  iteration,
                         int  niter,
                         Real stop_time)
{
  assert(level == 0);

  // XXX

  // for (int k=0; k<max_iters; k++) {
  //   sdc_mg_sweep(&mg, time, dt_level[0], 0);
  // }
  // Amr::timeStep(level, time, iteration, niter, stop_time);
}


MLSDCAmr::MLSDCAmr()
{
  Amr::Amr();
  printf("DOUBLE YAR!\n");

  //
  // get parameters
  //
  ParmParse ppsdc("sdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 22;
  
  //
  // build MultiFab encapsulation for SDCLib
  //
  encap_ctx.ncomp = 1;
  encap_ctx.ngrow = 1;
  build_encap();

  //
  // build multigrid sdc sweeper
  //
  // sdc_mg_build(&mg, max_level);
  
  // XXX: add levels to mg
  
  // sdc_mg_setup(&mg);
  // for (int l=0; l<mg.nlevels; l++)
  //   mg.sweepers[l]->nset->flags |= SDC_ALLOC_FAKE;


  //
  // set subcycling_mode to "None"
  //
  sub_cycle = false;
  for (int i=0; i<=max_level; i++)
    n_cycle[i] = 1;
}


MLSDCAmr::~MLSDCAmr()
{
  // sdc_mg_destroy(&mg);
}
