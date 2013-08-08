
#ifndef _MLSDCAMR_H_
#define _MLSDCAMR_H_

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C   }

#include <MultiFab.H>
#include <Amr.H>
#include <sdc.h>


/*
 * MultiFab encapsulation context for SDCLib.
 */

typedef struct {
  BoxArray& ba;
  int       ncomp;
  int       ngrow;
} mf_encap_t;


/*
 * Multigrid SDC+AMR hybrid AMR class.
 */

class MLSDCAmr : public Amr
{
public:

  int max_iters;
  int current_level;

  sdc_mg_t     mg;
  sdc_encap_t  encap;
  mf_encap_t   encap_ctx;


  MLSDCAmr();
  ~MLSDCAmr();

  void build_encap();

  virtual void timeStep(int  level,
                        Real time,
                        int  iteration,
                        int  niter,
                        Real stop_time);


};

#endif
