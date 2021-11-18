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

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>

#include "NVector_Multifab.h"

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

N_Vector N_VNewEmpty_Multifab(sunindextype length)
{
   N_Vector v;
   N_Vector_Ops ops;
   N_VectorContent_Multifab content;

   /* Create vector */
   v = NULL;
   v = (N_Vector) malloc(sizeof *v);
   if (v == NULL) return(NULL);

   /* Create vector operation structure */
   ops = NULL;
   ops = (N_Vector_Ops) malloc(sizeof *ops);
   if (ops == NULL) { free(v); return(NULL); }

   ops->nvgetvectorid     = NULL;
   ops->nvclone           = N_VClone_Multifab;
   ops->nvcloneempty      = N_VCloneEmpty_Multifab;
   ops->nvdestroy         = N_VDestroy_Multifab;
   ops->nvspace           = N_VSpace_Multifab;
   ops->nvgetarraypointer = NULL;
   ops->nvsetarraypointer = NULL;

   /* standard vector operations */
   ops->nvlinearsum    = N_VLinearSum_Multifab;
   ops->nvconst        = N_VConst_Multifab;
   ops->nvprod         = N_VProd_Multifab;
   ops->nvdiv          = N_VDiv_Multifab;
   ops->nvscale        = N_VScale_Multifab;
   ops->nvabs          = N_VAbs_Multifab;
   ops->nvinv          = N_VInv_Multifab;
   ops->nvaddconst     = N_VAddConst_Multifab;
   ops->nvdotprod      = N_VDotProd_Multifab;
   ops->nvmaxnorm      = N_VMaxNorm_Multifab;
   ops->nvwrmsnormmask = N_VWrmsNormMask_Multifab;
   ops->nvwrmsnorm     = N_VWrmsNorm_Multifab;
   ops->nvmin          = N_VMin_Multifab;
   ops->nvwl2norm      = N_VWL2Norm_Multifab;
   ops->nvl1norm       = N_VL1Norm_Multifab;
   ops->nvcompare      = N_VCompare_Multifab;
   ops->nvinvtest      = N_VInvTest_Multifab;
   ops->nvconstrmask   = N_VConstrMask_Multifab;
   ops->nvminquotient  = N_VMinQuotient_Multifab;

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
   content = (N_VectorContent_Multifab) malloc(sizeof *content);
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

N_Vector N_VNew_Multifab(sunindextype length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost)
{
   N_Vector v;

   v = NULL;
   v = N_VNewEmpty_Multifab(length);
   if (v == NULL) return(NULL);

   // Create and attach new MultiFab
   if (length > 0)
   {
      amrex::MultiFab *mf_v = new amrex::MultiFab(ba, dm, nComp, nGhost);
      NV_OWN_MF_M(v) = SUNTRUE;
      NV_MFAB(v)     = mf_v;
   }

   return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a MultiFab N_Vector with user-specific MultiFab
 */

N_Vector N_VMake_Multifab(sunindextype length, amrex::MultiFab *v_mf)
{
   N_Vector v;

   v = NULL;
   v = N_VNewEmpty_Multifab(length);
   if (v == NULL) return(NULL);

   if (length > 0)
   {
      // Attach MultiFab
      NV_OWN_MF_M(v) = SUNFALSE;
      NV_MFAB(v)     = v_mf;
   }

   return(v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Multifab(N_Vector v)
{
   return NV_LENGTH_M(v);
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Multifab(N_Vector w)
{
   N_Vector v;
   N_Vector_Ops ops;
   N_VectorContent_Multifab content;

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
   content = (N_VectorContent_Multifab) malloc(sizeof *content);
   if (content == NULL) { free(ops); free(v); return(NULL); }

   content->length = NV_LENGTH_M(w);
   content->own_mf = SUNFALSE;
   content->mf     = NULL;

   /* Attach content and ops */
   v->content = content;
   v->ops     = ops;

   return(v);
}

N_Vector N_VClone_Multifab(N_Vector w)
{
   N_Vector v;
   sunindextype length;

   v = NULL;
   v = N_VCloneEmpty_Multifab(w);
   if (v == NULL) return(NULL);

   length = NV_LENGTH_M(w);

   if (length > 0)
   {
      // Copy the multifab
      amrex::MultiFab *mf_w = NV_MFAB(w);
      const amrex::BoxArray &ba = mf_w->boxArray();
      const amrex::DistributionMapping &dm = mf_w->DistributionMap();
      int nComp = mf_w->nComp();
      int nGhost = mf_w->nGrow();  // same number of ghost cells in the clone
      amrex::MultiFab *mf_v = new amrex::MultiFab(ba, dm, nComp, nGhost);

      // Attach multifab
      NV_OWN_MF_M(v) = SUNTRUE;
      NV_MFAB(v)     = mf_v;
   }

   return(v);
}

void N_VDestroy_Multifab(N_Vector v)
{
   if (NV_OWN_MF_M(v) == SUNTRUE)
   {
      delete NV_MFAB(v);
      NV_MFAB(v) = NULL;
   }
   free(v->content); v->content = NULL;
   free(v->ops); v->ops = NULL;
   free(v); v = NULL;

   return;
}

void N_VSpace_Multifab(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
   *lrw = NV_LENGTH_M(v);
   *liw = 1;

   return;
}

void N_VLinearSum_Multifab(realtype a, N_Vector x, realtype b, N_Vector y,
                           N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_y = NV_MFAB(y);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::LinComb(*mf_z, a, *mf_x, 0, b, *mf_y, 0, 0, ncomp, nghost);
}

void N_VConst_Multifab(realtype c, N_Vector z)
{
   sunindextype i, N;
   amrex::MultiFab *mf_z = NV_MFAB(z);
   *mf_z = c;
}

void N_VProd_Multifab(N_Vector x, N_Vector y, N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_y = NV_MFAB(y);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   amrex::MultiFab::Multiply(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VDiv_Multifab(N_Vector x, N_Vector y, N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_y = NV_MFAB(y);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   amrex::MultiFab::Divide(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VScale_Multifab(realtype c, N_Vector x, N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->mult(c, 0, ncomp, nghost);
}

void N_VAbs_Multifab(N_Vector x, N_Vector z)
{
   using namespace amrex;

   MultiFab *mf_x = NV_MFAB(x);
   MultiFab *mf_z = NV_MFAB(z);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& z_fab = mf_z->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  z_fab(i,j,k,c) = SUNRabs(x_fab(i,j,k,c));
               }
            }
         }
      }
   }
}

void N_VInv_Multifab(N_Vector x, N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->invert(1.0, 0, ncomp, nghost);
}

void N_VAddConst_Multifab(N_Vector x, realtype b, N_Vector z)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_z = NV_MFAB(z);

   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells

   amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
   mf_z->plus(b, nghost);
}

realtype N_VDotProd_Multifab(N_Vector x, N_Vector y)
{
   using namespace amrex;

   MultiFab *mf_x = NV_MFAB(x);
   MultiFab *mf_y = NV_MFAB(y);
   sunindextype ncomp = mf_x->nComp();
   //sunindextype nghost = mf_x->nGrow();
   sunindextype nghost = 0;  // do not include ghost cells in dot product

   realtype dotproduct = amrex::MultiFab::Dot(*mf_x, 0, *mf_y, 0, ncomp, nghost);

   return dotproduct;
}

realtype N_VMaxNorm_Multifab(N_Vector x)
{
   using namespace amrex;

   MultiFab *mf_x = NV_MFAB(x);
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

realtype N_VWrmsNorm_Multifab(N_Vector x, N_Vector w)
{
   using namespace amrex;

   MultiFab *mf_x = NV_MFAB(x);
   MultiFab *mf_w = NV_MFAB(w);
   sunindextype ncomp = mf_x->nComp();
   sunindextype N = NV_LENGTH_M(x);
   realtype sum = ZERO;
   realtype prodi;

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& w_fab = mf_w->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  prodi = x_fab(i,j,k,c) * w_fab(i,j,k,c);
                  sum += SUNSQR(prodi);
               }
            }
         }
      }
   }

   ParallelDescriptor::ReduceRealSum(sum);

   return SUNRsqrt(sum/N);
}

realtype N_VWrmsNormMask_Multifab(N_Vector x, N_Vector w, N_Vector id)
{
   using namespace amrex;

   MultiFab *mf_x = NV_MFAB(x);
   MultiFab *mf_w = NV_MFAB(w);
   MultiFab *mf_id = NV_MFAB(id);
   sunindextype ncomp = mf_x->nComp();
   sunindextype N = NV_LENGTH_M(x);
   realtype sum = ZERO;
   realtype prodi;

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& w_fab = mf_w->array(mfi);
      Array4<Real> const& id_fab = mf_id->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  if (id_fab(i,j,k,c) > ZERO)
                  {
                     prodi = x_fab(i,j,k,c) * w_fab(i,j,k,c);
                     sum += SUNSQR(prodi);
                  }
               }
            }
         }
      }
   }

   ParallelDescriptor::ReduceRealSum(sum);

   return SUNRsqrt(sum/N);
}

realtype N_VMin_Multifab(N_Vector x)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
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

realtype N_VWL2Norm_Multifab(N_Vector x, N_Vector w)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_w = NV_MFAB(w);
   sunindextype ncomp = mf_x->nComp();

   // do not include ghost cells in norm
   realtype sum = ZERO;
   realtype prodi;
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& w_fab = mf_w->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  prodi = x_fab(i,j,k,c) * w_fab(i,j,k,c);
                  sum += SUNSQR(prodi);
               }
            }
         }
      }
   }

   amrex::ParallelDescriptor::ReduceRealSum(sum);

   return SUNRsqrt(sum);
}

realtype N_VL1Norm_Multifab(N_Vector x)
{
   amrex::MultiFab *mf_x = NV_MFAB(x);
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

void N_VCompare_Multifab(realtype a, N_Vector x, N_Vector z)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_z = NV_MFAB(z);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& z_fab = mf_z->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  z_fab(i,j,k,c) = (SUNRabs(x_fab(i,j,k,c)) >= a) ? ONE : ZERO;
               }
            }
         }
      }
   }

   return;
}

booleantype N_VInvTest_Multifab(N_Vector x, N_Vector z)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_z = NV_MFAB(z);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   realtype val = ONE;
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& z_fab = mf_z->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  if(x_fab(i,j,k,c) == ZERO)
                  {
                     val = ZERO;
                  }
                  else
                  {
                     z_fab(i,j,k,c) = ONE / x_fab(i,j,k,c);
                  }
               }
            }
         }
      }
   }

   amrex::ParallelDescriptor::ReduceRealMin(val);

   if (val == ZERO)
   {
      return SUNFALSE;
   }
   else
   {
      return SUNTRUE;
   }
}

booleantype N_VConstrMask_Multifab(N_Vector a, N_Vector x, N_Vector m)
{
   using namespace amrex;

   amrex::MultiFab *mf_x = NV_MFAB(x);
   amrex::MultiFab *mf_a = NV_MFAB(a);
   amrex::MultiFab *mf_m = NV_MFAB(m);
   sunindextype ncomp = mf_x->nComp();

   // ghost cells not included
   realtype temp = ZERO;
   for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& x_fab = mf_x->array(mfi);
      Array4<Real> const& a_fab = mf_a->array(mfi);
      Array4<Real> const& m_fab = mf_a->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  m_fab(i,j,k,c) = ZERO;

                  /* Continue if no constraints were set for the variable */
                  if (a_fab(i,j,k,c) == ZERO)
                     continue;

                  /* Check if a set constraint has been violated */
                  realtype a, x;
                  a = a_fab(i,j,k,c);
                  x = x_fab(i,j,k,c);
                  booleantype test = (SUNRabs(a) > ONEPT5 && x*a <= ZERO) ||
                     (SUNRabs(a) > HALF   && x*a <  ZERO);
                  if (test) {
                     temp = m_fab(i,j,k,c) = ONE;
                  }
               }
            }
         }
      }
   }

   /* Return false if any constraint was violated */
   amrex::ParallelDescriptor::ReduceRealMax(temp);

   return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_Multifab(N_Vector num, N_Vector denom)
{
   using namespace amrex;

   amrex::MultiFab *mf_num = NV_MFAB(num);
   amrex::MultiFab *mf_denom = NV_MFAB(denom);
   sunindextype ncomp = mf_num->nComp();

   // ghost cells not included
   realtype min = BIG_REAL;
   booleantype notEvenOnce = SUNTRUE;
   for (MFIter mfi(*mf_num); mfi.isValid(); ++mfi)
   {
      const amrex::Box& bx = mfi.validbox();
      Array4<Real> const& num_fab = mf_num->array(mfi);
      Array4<Real> const& denom_fab = mf_denom->array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int c = 0; c < ncomp; ++c) {
         for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  if (denom_fab(i,j,k,c) == ZERO)
                  {
                     continue;
                  }
                  else
                  {
                     realtype num = num_fab(i,j,k,c);
                     realtype denom = denom_fab(i,j,k,c);
                     if (!notEvenOnce)
                     {
                        min = SUNMIN(min, num / denom);
                     }
                     else
                     {
                        min = num / denom;
                        notEvenOnce = SUNFALSE;
                     }
                  }
               }
            }
         }
      }
   }

   amrex::ParallelDescriptor::ReduceRealMin(min);

   return min;
}
