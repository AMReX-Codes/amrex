/*------------------------------------------------------------------------------
  Header file for N_Vector wrap of an AMReX 'MultiFab'. Based on example codes
  for the 2019 Argonne Training Program in Extreme-Scale Computing with SUNDIALS
  and AMReX.

  Authors (alphabetical):
    David Gardner (gardner48@llnl.gov)
    John Loffeld (loffeld1@llnl.gov)
    Daniel Reynolds (reynolds@smu.edu)
    Donald Willcox (dewillcox@lbl.gov)
  ----------------------------------------------------------------------------*/
#ifndef AMREX_NVECTOR_MULTIFAB_H_
#define AMREX_NVECTOR_MULTIFAB_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Sundials_Core.H>

#include <sundials/sundials_nvector.h>
#include <cstdio>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
namespace amrex::sundials {
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Vector content structure
 * -----------------------------------------------------------------*/

struct N_VectorContent_MultiFab_notptr {
  sunindextype length;   /* vector length           */
  int own_mf;            /* MultiFab ownership flag */
  amrex::MultiFab *mf;   /* wrapped MultiFab        */
};

typedef struct N_VectorContent_MultiFab_notptr *N_VectorContent_MultiFab;

/* -----------------------------------------------------------------
 * Exported functions
 * -----------------------------------------------------------------*/

N_Vector N_VNewEmpty_MultiFab(sunindextype vec_length,
                              ::sundials::Context* sunctx = The_Sundials_Context());
N_Vector N_VNew_MultiFab(sunindextype vec_length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost,
                         ::sundials::Context* sunctx = The_Sundials_Context());
N_Vector N_VMake_MultiFab(sunindextype vec_length,
                          amrex::MultiFab *mf,
                          ::sundials::Context* sunctx = The_Sundials_Context());
sunindextype N_VGetLength_MultiFab(N_Vector v);
int N_VGetOwnMF_MultiFab(N_Vector v);
void N_VSetOwnMF_MultiFab(N_Vector v, int own_mf);
N_Vector N_VCloneEmpty_MultiFab(N_Vector w);
N_Vector N_VClone_MultiFab(N_Vector w);
void N_VDestroy_MultiFab(N_Vector v);
void N_VSpace_MultiFab(N_Vector v, sunindextype *lrw,
                       sunindextype *liw);
N_VectorContent_MultiFab N_VGetContent_MultiFab(N_Vector v);
MultiFab* N_VGetVectorPointer_MultiFab(N_Vector v);
amrex::Real NormHelper_NVector_MultiFab(N_Vector x, N_Vector w, N_Vector id, int use_id, bool rms);
/* standard vector operations */
void N_VLinearSum_MultiFab(amrex::Real a, N_Vector x,
                           amrex::Real b, N_Vector y, N_Vector z);
void N_VConst_MultiFab(amrex::Real c, N_Vector z);
void N_VProd_MultiFab(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_MultiFab(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_MultiFab(amrex::Real c, N_Vector x, N_Vector z);
void N_VAbs_MultiFab(N_Vector x, N_Vector z);
void N_VInv_MultiFab(N_Vector x, N_Vector z);
void N_VAddConst_MultiFab(N_Vector x, amrex::Real b, N_Vector z);
amrex::Real N_VDotProd_MultiFab(N_Vector x, N_Vector y);
amrex::Real N_VMaxNorm_MultiFab(N_Vector x);
amrex::Real N_VWrmsNorm_MultiFab(N_Vector x, N_Vector w);
amrex::Real N_VWrmsNormMask_MultiFab(N_Vector x, N_Vector w,
                                     N_Vector id);
amrex::Real N_VMin_MultiFab(N_Vector x);
amrex::Real N_VWL2Norm_MultiFab(N_Vector x, N_Vector w);
amrex::Real N_VL1Norm_MultiFab(N_Vector x);
void N_VCompare_MultiFab(amrex::Real a, N_Vector x, N_Vector z);
int N_VInvTest_MultiFab(N_Vector x, N_Vector z);
int N_VConstrMask_MultiFab(N_Vector a_a, N_Vector x,
                           N_Vector m);
amrex::Real N_VMinQuotient_MultiFab(N_Vector num, N_Vector denom);

#ifdef __cplusplus
} // extern "C"

MultiFab*& getMFptr(N_Vector v);
MultiFab N_VGetVectorAlias_MultiFab(N_Vector v);

}//namespace amrex::sundials
#endif

#endif
