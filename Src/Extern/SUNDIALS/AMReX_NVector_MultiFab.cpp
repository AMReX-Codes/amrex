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
  Implementation file for N_Vector wrap of AMReX 'MultiFab' structure.
  --------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include <sundials/sundials_math.h>

#include "AMReX_NVector_MultiFab.H"

namespace amrex {
namespace sundials {

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Function to create a new empty multifab vector
 */

N_Vector N_VNewEmpty_MultiFab(sunindextype length)
{
   N_Vector v;
   N_Vector_Ops ops;
   N_VectorContent_MultiFab content;

   /* Create vector */
   v = NULL;
   v = (N_Vector) malloc(sizeof *v);
   if (v == NULL) return(NULL);

   /* Create vector operation structure */
   ops = NULL;
   ops = (N_Vector_Ops) malloc(sizeof *ops);
   if (ops == NULL) { free(v); return(NULL); }

   ops->nvgetvectorid     = NULL;
   ops->nvclone           = N_VClone_MultiFab;
   ops->nvcloneempty      = N_VCloneEmpty_MultiFab;
   ops->nvdestroy         = N_VDestroy_MultiFab;
   ops->nvspace           = N_VSpace_MultiFab;
   ops->nvgetarraypointer = NULL;
   ops->nvsetarraypointer = NULL;
   ops->nvgetlength       = N_VGetLength_MultiFab;

   /* standard vector operations */
   ops->nvlinearsum    = N_VLinearSum_MultiFab;
   ops->nvconst        = N_VConst_MultiFab;
   ops->nvprod         = N_VProd_MultiFab;
   ops->nvdiv          = N_VDiv_MultiFab;
   ops->nvscale        = N_VScale_MultiFab;
   ops->nvabs          = N_VAbs_MultiFab;
   ops->nvinv          = N_VInv_MultiFab;
   ops->nvaddconst     = N_VAddConst_MultiFab;
   ops->nvdotprod      = N_VDotProd_MultiFab;
   ops->nvmaxnorm      = N_VMaxNorm_MultiFab;
   ops->nvwrmsnormmask = N_VWrmsNormMask_MultiFab;
   ops->nvwrmsnorm     = N_VWrmsNorm_MultiFab;
   ops->nvmin          = N_VMin_MultiFab;
   ops->nvwl2norm      = N_VWL2Norm_MultiFab;
   ops->nvl1norm       = N_VL1Norm_MultiFab;
   ops->nvcompare      = N_VCompare_MultiFab;
   ops->nvinvtest      = N_VInvTest_MultiFab;
   ops->nvconstrmask   = N_VConstrMask_MultiFab;
   ops->nvminquotient  = N_VMinQuotient_MultiFab;

   /* fused vector operations (optional, NULL means disabled by default) */
   ops->nvlinearcombination = NULL;
   ops->nvscaleaddmulti     = NULL;
   ops->nvdotprodmulti      = NULL;

   /* vector array operations (optional, NULL means disabled by default) */
   ops->nvlinearsumvectorarray         = NULL;
   ops->nvscalevectorarray             = NULL;
   ops->nvconstvectorarray             = NULL;
   ops->nvwrmsnormvectorarray          = NULL;
   ops->nvwrmsnormmaskvectorarray      = NULL;
   ops->nvscaleaddmultivectorarray     = NULL;
   ops->nvlinearcombinationvectorarray = NULL;

   /* Create content */
   content = NULL;
   content = (N_VectorContent_MultiFab) malloc(sizeof *content);
   if (content == NULL) { free(ops); free(v); return(NULL); }

   content->length = length;
   content->own_mf = SUNFALSE;
   content->mf     = NULL;

   /* Attach content and ops */
   v->content = content;
   v->ops     = ops;

   return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new MultiFab vector
 */

N_Vector N_VNew_MultiFab(sunindextype length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost)
{
   N_Vector v;

   v = NULL;
   v = N_VNewEmpty_MultiFab(length);
   if (v == NULL) return(NULL);

   // Create and attach new MultiFab
   if (length > 0)
   {
      amrex::MultiFab *mf_v = new amrex::MultiFab(ba, dm, nComp, nGhost);
      AMREX_NV_OWN_MF_M(v) = SUNTRUE;
      amrex::sundials::getMFptr(v)     = mf_v;
   }

   return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a MultiFab N_Vector with user-specific MultiFab
 */

N_Vector N_VMake_MultiFab(sunindextype length, amrex::MultiFab *v_mf)
{
   N_Vector v;

   v = NULL;
   v = N_VNewEmpty_MultiFab(length);
   if (v == NULL) return(NULL);

   if (length > 0)
   {
      // Attach MultiFab
      AMREX_NV_OWN_MF_M(v) = SUNFALSE;
      amrex::sundials::getMFptr(v)     = v_mf;
   }

   return(v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_MultiFab(N_Vector v)
{
   return AMREX_NV_LENGTH_M(v);
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_MultiFab(N_Vector w)
{
   N_Vector v;
   N_Vector_Ops ops;
   N_VectorContent_MultiFab content;

   if (w == NULL) return(NULL);

   /* Create vector */
   v = NULL;
   v = (N_Vector) malloc(sizeof *v);
   if (v == NULL) return(NULL);

   /* Create vector operation structure */
   ops = NULL;
   ops = (N_Vector_Ops) malloc(sizeof *ops);
   if (ops == NULL) { free(v); return(NULL); }

   ops->nvgetvectorid     = w->ops->nvgetvectorid;
   ops->nvclone           = w->ops->nvclone;
   ops->nvcloneempty      = w->ops->nvcloneempty;
   ops->nvdestroy         = w->ops->nvdestroy;
   ops->nvspace           = w->ops->nvspace;
   ops->nvgetarraypointer = w->ops->nvgetarraypointer;
   ops->nvsetarraypointer = w->ops->nvsetarraypointer;
   ops->nvgetlength       = w->ops->nvgetlength;

   /* standard vector operations */
   ops->nvlinearsum    = w->ops->nvlinearsum;
   ops->nvconst        = w->ops->nvconst;
   ops->nvprod         = w->ops->nvprod;
   ops->nvdiv          = w->ops->nvdiv;
   ops->nvscale        = w->ops->nvscale;
   ops->nvabs          = w->ops->nvabs;
   ops->nvinv          = w->ops->nvinv;
   ops->nvaddconst     = w->ops->nvaddconst;
   ops->nvdotprod      = w->ops->nvdotprod;
   ops->nvmaxnorm      = w->ops->nvmaxnorm;
   ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
   ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
   ops->nvmin          = w->ops->nvmin;
   ops->nvwl2norm      = w->ops->nvwl2norm;
   ops->nvl1norm       = w->ops->nvl1norm;
   ops->nvcompare      = w->ops->nvcompare;
   ops->nvinvtest      = w->ops->nvinvtest;
   ops->nvconstrmask   = w->ops->nvconstrmask;
   ops->nvminquotient  = w->ops->nvminquotient;

   /* fused vector operations */
   ops->nvlinearcombination = w->ops->nvlinearcombination;
   ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
   ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

   /* vector array operations */
   ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
   ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
   ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
   ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
   ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
   ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
   ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

   /* Create content */
   content = NULL;
   content = (N_VectorContent_MultiFab) malloc(sizeof *content);
   if (content == NULL) { free(ops); free(v); return(NULL); }

   content->length = AMREX_NV_LENGTH_M(w);
   content->own_mf = SUNFALSE;
   content->mf     = NULL;

   /* Attach content and ops */
   v->content = content;
   v->ops     = ops;

   return(v);
}

N_Vector N_VClone_MultiFab(N_Vector w)
{
   N_Vector v;
   sunindextype length;

   v = NULL;
   v = N_VCloneEmpty_MultiFab(w);
   if (v == NULL) return(NULL);

   length = AMREX_NV_LENGTH_M(w);

   if (length > 0)
   {
      // Copy the multifab
      amrex::MultiFab *mf_w = amrex::sundials::getMFptr(w);
      const amrex::BoxArray &ba = mf_w->boxArray();
      const amrex::DistributionMapping &dm = mf_w->DistributionMap();
      int nComp = mf_w->nComp();
      int nGhost = mf_w->nGrow();  // same number of ghost cells in the clone
      amrex::MultiFab *mf_v = new amrex::MultiFab(ba, dm, nComp, nGhost);

      // Attach multifab
      AMREX_NV_OWN_MF_M(v) = SUNTRUE;
      amrex::sundials::getMFptr(v)     = mf_v;
   }

   return(v);
}

void N_VDestroy_MultiFab(N_Vector v)
{
   if (AMREX_NV_OWN_MF_M(v) == SUNTRUE)
   {
      delete amrex::sundials::getMFptr(v);
      amrex::sundials::getMFptr(v) = NULL;
   }
   free(v->content); v->content = NULL;
   free(v->ops); v->ops = NULL;
   free(v); v = NULL;

   return;
}

void N_VSpace_MultiFab(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
   *lrw = AMREX_NV_LENGTH_M(v);
   *liw = 1;

   return;
}

N_VectorContent_MultiFab N_VGetContent_MultiFab(N_Vector v)
{
  return (N_VectorContent_MultiFab)(v->content);
}

/* ----------------------------------------------------------------
 * Extract MultiFab*
 */
amrex::MultiFab*& getMFptr(N_Vector v)
{
  return ((N_VectorContent_MultiFab)(v->content) )->mf;
}

amrex::MultiFab* N_VGetVectorPointer_MultiFab(N_Vector v)
{
  return ((N_VectorContent_MultiFab)(v->content) )->mf;
}

/* ----------------------------------------------------------------
 * Extract alias MultiFab
 */

amrex::MultiFab N_VGetVectorAlias_MultiFab(N_Vector v)
{
    return amrex::MultiFab(*((N_VectorContent_MultiFab)(v->content) )->mf,amrex::make_alias,0,(((N_VectorContent_MultiFab)(v->content) )->mf)->nComp());
}

void N_VLinearSum_MultiFab(realtype a, N_Vector x, realtype b, N_Vector y,
                           N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::LinComb(*mf_z, a, *mf_x, 0, b, *mf_y, 0, 0, ncomp, nghost);
}

void N_VConst_MultiFab(realtype c, N_Vector z)
{
   sunindextype i, N;
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);
   *mf_z = c;
}

void N_VProd_MultiFab(N_Vector x, N_Vector y, N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   amrex::MultiFab::Multiply(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VDiv_MultiFab(N_Vector x, N_Vector y, N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   amrex::MultiFab::Divide(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VScale_MultiFab(realtype c, N_Vector x, N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->mult(c, 0, ncomp, nghost);
}

void N_VAbs_MultiFab(N_Vector x, N_Vector z)
{
   using namespace amrex;

   MultiFab *mf_x = amrex::sundials::getMFptr(x);
   MultiFab *mf_z = amrex::sundials::getMFptr(z);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& z_fab = mf_z->array(mfi);

      amrex::ParallelFor(bx, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
      {
          z_fab(i,j,k,c) = SUNRabs(x_fab(i,j,k,c));
      });
   }
}

void N_VInv_MultiFab(N_Vector x, N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->invert(1.0, 0, ncomp, nghost);
}

void N_VAddConst_MultiFab(N_Vector x, realtype b, N_Vector z)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->plus(b, nghost);
}

realtype N_VDotProd_MultiFab(N_Vector x, N_Vector y)
{
   using namespace amrex;

   MultiFab *mf_x = amrex::sundials::getMFptr(x);
   MultiFab *mf_y = amrex::sundials::getMFptr(y);
   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells in dot product

   realtype dotproduct = amrex::MultiFab::Dot(*mf_x, 0, *mf_y, 0, ncomp, nghost);

   return dotproduct;
}

realtype N_VMaxNorm_MultiFab(N_Vector x)
{
   using namespace amrex;

   MultiFab *mf_x = amrex::sundials::getMFptr(x);
   sunindextype ncomp = mf_x->nComp();
   sunindextype startComp = 0;
   sunindextype nghost = 0;  // do not include ghost cells in the norm

   realtype max = mf_x->max(startComp, nghost);

   // continue with rest of comps
   for (int c = 1; c < ncomp; ++c)
   {
      realtype comp_max = mf_x->max(c, nghost); // comp c, no ghost zones
      if (comp_max > max)
      {
         max = comp_max;
      }
   }

   // no reduction needed, done in multifab
   return max;
}

realtype N_VWrmsNorm_MultiFab(N_Vector x, N_Vector w)
{

   using namespace amrex;
   sunindextype N = AMREX_NV_LENGTH_M(x);
   return N_VWL2Norm_MultiFab(x, w)*std::sqrt(1.0_rt/N);
}

realtype N_VWrmsNormMask_MultiFab(N_Vector x, N_Vector w, N_Vector id)
{
   using namespace amrex;

   return NormHelper_MultiFab(x, w, id, true, true);
}

realtype N_VMin_MultiFab(N_Vector x)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   sunindextype ncomp = mf_x->nComp();

   sunindextype startComp = 0;
   sunindextype nghost = 0;  // ghost zones not included in min

   realtype min = mf_x->min(startComp, nghost);

   // continue with rest of comps
   for (int c = 1; c < ncomp; ++c)
   {
      realtype comp_min = mf_x->min(c, nghost); // comp c, no ghost zones
      if (comp_min < min)
      {
         min = comp_min;
      }
   }

   // no reduction needed, done in multifab
   return min;
}

realtype NormHelper_MultiFab(N_Vector x, N_Vector w, N_Vector id, int use_id, bool rms)
{
   using namespace amrex;

   MultiFab *mf_x = amrex::sundials::getMFptr(x);
   MultiFab *mf_w = amrex::sundials::getMFptr(w);
   MultiFab *mf_id = use_id ? amrex::sundials::getMFptr(id) : NULL;
   sunindextype numcomp = mf_x->nComp();
   sunindextype N = AMREX_NV_LENGTH_M(x);
   realtype prodi;
   bool local = true;
   int nghost = 0;

   // ghost cells not included
    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = mf_x->const_arrays();
        auto const& wma = mf_w->const_arrays();
        auto const& idma = use_id ? mf_id->const_arrays() : mf_x->const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, *mf_x, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            auto const& x_fab = xma[box_no];
            auto const& w_fab = wma[box_no];
            auto const& id_fab = use_id ? idma[box_no] : xma[box_no];
            for (int n = 0; n < numcomp; ++n) {
	      if(use_id)
		t += (id_fab(i,j,k,n) > ZERO) ? std::sqrt(x_fab(i,j,k,n) * w_fab(i,j,k,n)) : 0.0;
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(*mf_x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& x_fab = mf_x->const_array(mfi);
            Array4<Real const> const& w_fab = mf_w->const_array(mfi);
	    Array4<Real const> const& id_fab = use_id ? mf_id->const_array(mfi) : mf_x->const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
      	      if(use_id)
		sm += id_fab(i,j,k,n) > ZERO ? std::sqrt(x_fab(i,j,k,n) * w_fab(i,j,k,n)) : 0.0_rt;
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    Real sum=sm;
   ParallelDescriptor::ReduceRealSum(sum);

   return rms ? SUNRsqrt(sum/N) : SUNRsqrt(sum);
}

realtype N_VWL2Norm_MultiFab(N_Vector x, N_Vector w)
{
   using namespace amrex;

   return NormHelper_MultiFab(x, w, N_VCloneEmpty_MultiFab(x), false, false);
   /*
   MultiFab *mf_x = amrex::sundials::getMFptr(x);
   MultiFab *mf_w = amrex::sundials::getMFptr(w);
   sunindextype numcomp = mf_x->nComp();
   sunindextype N = AMREX_NV_LENGTH_M(x);
   realtype prodi;
   bool local = true;
   int nghost = 0;

   // ghost cells not included
    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = mf_x->const_arrays();
        auto const& wma = mf_w->const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, x, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            auto const& xfab = xma[box_no];
            auto const& wfab = wma[box_no];
            for (int n = 0; n < numcomp; ++n) {
	      sm += std::sqrt(x_fab(i,j,k,n) * w_fab(i,j,k,n));
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(*mf_x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& x_fab = mf_x->const_array(mfi);
            Array4<Real const> const& w_fab = mf_w->const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
	        sm += std::sqrt(x_fab(i,j,k,n) * w_fab(i,j,k,n));
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    Real sum=sm;
   ParallelDescriptor::ReduceRealSum(sum);

   return SUNRsqrt(sum);
   */
}

realtype N_VL1Norm_MultiFab(N_Vector x)
{
   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   sunindextype ncomp = mf_x->nComp();

   sunindextype startComp = 0;
   sunindextype nghost = 0;  // ghost zones not included in norm

   realtype sum = mf_x->norm1(startComp, nghost);

   // continue with rest of comps
   for (int c = 1; c < ncomp; ++c)
   {
      realtype comp_sum = mf_x->norm1(c, nghost);
      sum += comp_sum;
   }

   // no reduction needed, it was done in multifab
   return sum;
}

void N_VCompare_MultiFab(realtype a, N_Vector x, N_Vector z)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& z_fab = mf_z->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      amrex::ParallelFor(bx, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
      {
        z_fab(i,j,k,c) = (SUNRabs(x_fab(i,j,k,c)) >= a) ? ONE : ZERO;
      });
   }

   return;
}

booleantype N_VInvTest_MultiFab(N_Vector x, N_Vector z)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);
   sunindextype ncomp = mf_x->nComp();
   int nghost = 0;

   //Lazily assume we're ok dividing by 0
   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->invert(1.0, 0, ncomp, nghost);

   // ghost cells not included
   bool val = amrex::ReduceLogicalAnd(*mf_x, *mf_z, nghost,
	       [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& x_fab, Array4<Real const> const& z_fab) -> bool
    {
    bool val_loc = true;
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  if (x_fab(i,j,k,c) == ZERO)
                  {
		    val_loc &= false;
		  }
		  /*
		  else
		  {
                     z_fab(i,j,k,c) = ONE / x_fab(i,j,k,c);
                  }*/
               }
            }
         }
      }
      return val_loc;
    });

   amrex::ParallelDescriptor::ReduceBoolAnd(val);

   if (val == false)
   {
      return SUNFALSE;
   }
   else
   {
      return SUNTRUE;
   }
}

booleantype N_VConstrMask_MultiFab(N_Vector a, N_Vector x, N_Vector m)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
   amrex::MultiFab *mf_a = amrex::sundials::getMFptr(a);
   amrex::MultiFab *mf_m = amrex::sundials::getMFptr(m);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   realtype temp = ZERO;
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& a_fab = mf_a->array(mfi);
      Array4<Real> const& m_fab = mf_a->array(mfi);
      //Changing continue to if check, temp calculation should be changed to reduction
      amrex::ParallelFor(bx, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
      {
                  m_fab(i,j,k,c) = ZERO;

                  /* Continue if no constraints were set for the variable */
                  if (a_fab(i,j,k,c) != ZERO)
                  {

                  /* Check if a set constraint has been violated */
                  realtype a, x;
                  a = a_fab(i,j,k,c);
                  x = x_fab(i,j,k,c);
                  booleantype test = (SUNRabs(a) > ONEPT5 && x*a <= ZERO) ||
                     (SUNRabs(a) > HALF   && x*a <  ZERO);
                  if (test) {
                     m_fab(i,j,k,c) = ONE;
                  }
                  }
      });
   }

   temp = mf_m->norm1();
   /* Return false if any constraint was violated */
   amrex::ParallelDescriptor::ReduceRealMax(temp);

   return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_MultiFab(N_Vector num, N_Vector denom)
{
   using namespace amrex;

   amrex::MultiFab *mf_num = amrex::sundials::getMFptr(num);
   amrex::MultiFab *mf_denom = amrex::sundials::getMFptr(denom);
   sunindextype ncomp = mf_num->nComp();

   // ghost cells not included
   int nghost = 0;
   Real min = amrex::ReduceMin(*mf_num, *mf_denom, nghost,
	       [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& num_fab, Array4<Real const> const& denom_fab) -> Real
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    Real min_loc = std::numeric_limits<Real>::max();
#else
    Real min_loc = Real(1.e37);
#endif
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  if (denom_fab(i,j,k,c) != ZERO)
                  {
                     realtype num = num_fab(i,j,k,c);
                     realtype denom = denom_fab(i,j,k,c);
                     min_loc = SUNMIN(min_loc, num / denom);
                  }
               }
            }
         }
      }
      return min_loc;
    });

   amrex::ParallelDescriptor::ReduceRealMin(min);

   return min;
}

}
}
