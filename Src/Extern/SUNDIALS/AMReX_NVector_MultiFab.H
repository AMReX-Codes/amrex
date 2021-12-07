/*--------------------------------------------------------------------
  Time Integration and Nonlinear Solvers
  Hands-on Lessons with SUNDIALS + AMReX
  2019 Argonne Training Program in Extreme-Scale Computing

  Authors (alphabetical):
    David Gardner (gardner48@llnl.gov)
    John Loffeld (loffeld1@llnl.gov)
    Daniel Reynolds (reynolds@smu.edu)
    Donald Willcox (dewillcox@lbl.gov)

  --------------------------------------------------------------------
  Header file for N_Vector wrap of AMReX 'MultiFab' structure.
  --------------------------------------------------------------------*/

#ifndef AMREX_NVECTOR_MULTIFAB_H
#define AMREX_NVECTOR_MULTIFAB_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include <sundials/sundials_nvector.h>
#include <cstdio>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

namespace amrex {
namespace sundials {

/* -----------------------------------------------------------------
 * Vector content structure
 * -----------------------------------------------------------------*/

struct N_VectorContent_MultiFab_notptr {
  sunindextype length;   /* vector length           */
  booleantype own_mf;    /* MultiFab ownership flag */
  amrex::MultiFab *mf;   /* wrapped MultiFab        */
};

typedef struct N_VectorContent_MultiFab_notptr *N_VectorContent_MultiFab;

/* -----------------------------------------------------------------
 * Accessor Macros
 * -----------------------------------------------------------------*/

#define AMREX_NV_CONTENT_M(v) ( (amrex::sundials::N_VectorContent_MultiFab)(v->content) )

#define AMREX_NV_LENGTH_M(v)  ( AMREX_NV_CONTENT_M(v)->length )

#define AMREX_NV_OWN_MF_M(v)  ( AMREX_NV_CONTENT_M(v)->own_mf )

/* -----------------------------------------------------------------
 * Exported functions
 * -----------------------------------------------------------------*/

N_Vector N_VNewEmpty_MultiFab(sunindextype vec_length);
N_Vector N_VNew_MultiFab(sunindextype vec_length,
                                         const amrex::BoxArray &ba,
                                         const amrex::DistributionMapping &dm,
                                         sunindextype nComp,
                                         sunindextype nGhost);
N_Vector N_VMake_MultiFab(sunindextype vec_length,
                                          amrex::MultiFab *mf);
sunindextype N_VGetLength_MultiFab(N_Vector v);
N_Vector N_VCloneEmpty_MultiFab(N_Vector w);
N_Vector N_VClone_MultiFab(N_Vector w);
void N_VDestroy_MultiFab(N_Vector v);
void N_VSpace_MultiFab(N_Vector v, sunindextype *lrw,
                                       sunindextype *liw);
N_VectorContent_MultiFab N_VGetContent_MultiFab(N_Vector v);
MultiFab*& getMFptr(N_Vector v);
MultiFab* N_VGetVectorPointer_MultiFab(N_Vector v);
MultiFab N_VGetVectorAlias_MultiFab(N_Vector v);
realtype NormHelper_MultiFab(N_Vector x, N_Vector w, N_Vector id, int use_id, bool rms);
/* standard vector operations */
void N_VLinearSum_MultiFab(realtype a, N_Vector x,
                                           realtype b, N_Vector y, N_Vector z);
void N_VConst_MultiFab(realtype c, N_Vector z);
void N_VProd_MultiFab(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_MultiFab(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_MultiFab(realtype c, N_Vector x, N_Vector z);
void N_VAbs_MultiFab(N_Vector x, N_Vector z);
void N_VInv_MultiFab(N_Vector x, N_Vector z);
void N_VAddConst_MultiFab(N_Vector x, realtype b, N_Vector z);
realtype N_VDotProd_MultiFab(N_Vector x, N_Vector y);
realtype N_VMaxNorm_MultiFab(N_Vector x);
realtype N_VWrmsNorm_MultiFab(N_Vector x, N_Vector w);
realtype N_VWrmsNormMask_MultiFab(N_Vector x, N_Vector w,
                                                  N_Vector id);
realtype N_VMin_MultiFab(N_Vector x);
realtype N_VWL2Norm_MultiFab(N_Vector x, N_Vector w);
realtype N_VL1Norm_MultiFab(N_Vector x);
void N_VCompare_MultiFab(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_MultiFab(N_Vector x, N_Vector z);
booleantype N_VConstrMask_MultiFab(N_Vector c, N_Vector x,
                                                   N_Vector m);
realtype N_VMinQuotient_MultiFab(N_Vector num, N_Vector denom);

}

}
#ifdef __cplusplus
}
#endif

#endif
