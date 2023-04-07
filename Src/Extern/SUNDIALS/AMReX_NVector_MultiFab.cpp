/*------------------------------------------------------------------------------
  Implementation file for N_Vector wrap of an AMReX 'MultiFab'. Based on example
  codes for the 2019 Argonne Training Program in Extreme-Scale Computing with
  SUNDIALS and AMReX.

  Authors (alphabetical):
    David Gardner (gardner48@llnl.gov)
    John Loffeld (loffeld1@llnl.gov)
    Daniel Reynolds (reynolds@smu.edu)
    Donald Willcox (dewillcox@lbl.gov)
  ----------------------------------------------------------------------------*/
#include "AMReX_NVector_MultiFab.H"
#include <type_traits>

namespace amrex::sundials {

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Function to create a new empty multifab vector
 */

N_Vector N_VNewEmpty_MultiFab(sunindextype length, ::sundials::Context* sunctx)
{
    /* Create vector */
    N_Vector v = N_VNewEmpty(*sunctx);
    if (v == nullptr) return(nullptr);

    v->ops->nvclone      = N_VClone_MultiFab;
    v->ops->nvcloneempty = N_VCloneEmpty_MultiFab;
    v->ops->nvdestroy    = N_VDestroy_MultiFab;
    v->ops->nvspace      = N_VSpace_MultiFab;
    v->ops->nvgetlength  = N_VGetLength_MultiFab;

    /* standard vector operations */
    v->ops->nvlinearsum    = N_VLinearSum_MultiFab;
    v->ops->nvconst        = N_VConst_MultiFab;
    v->ops->nvprod         = N_VProd_MultiFab;
    v->ops->nvdiv          = N_VDiv_MultiFab;
    v->ops->nvscale        = N_VScale_MultiFab;
    v->ops->nvabs          = N_VAbs_MultiFab;
    v->ops->nvinv          = N_VInv_MultiFab;
    v->ops->nvaddconst     = N_VAddConst_MultiFab;
    v->ops->nvdotprod      = N_VDotProd_MultiFab;
    v->ops->nvmaxnorm      = N_VMaxNorm_MultiFab;
    v->ops->nvwrmsnormmask = N_VWrmsNormMask_MultiFab;
    v->ops->nvwrmsnorm     = N_VWrmsNorm_MultiFab;
    v->ops->nvmin          = N_VMin_MultiFab;
    v->ops->nvwl2norm      = N_VWL2Norm_MultiFab;
    v->ops->nvl1norm       = N_VL1Norm_MultiFab;
    v->ops->nvcompare      = N_VCompare_MultiFab;
    v->ops->nvinvtest      = N_VInvTest_MultiFab;
    v->ops->nvconstrmask   = N_VConstrMask_MultiFab;
    v->ops->nvminquotient  = N_VMinQuotient_MultiFab;

    /* Create content */
    auto* content = (N_VectorContent_MultiFab)
        std::malloc(sizeof(std::remove_pointer_t<N_VectorContent_MultiFab>));
    if (content == nullptr) { N_VFreeEmpty(v); return(nullptr); }

    content->length = length;
    content->own_mf = SUNFALSE;
    content->mf     = nullptr;

    /* Attach content */
    v->content = content;

    return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new MultiFab vector
 */

N_Vector N_VNew_MultiFab(sunindextype length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost,
                         ::sundials::Context* sunctx)
{
    N_Vector v = N_VNewEmpty_MultiFab(length, sunctx);
    if (v == nullptr) return(nullptr);

    // Create and attach new MultiFab
    if (length > 0)
    {
         auto *mf_v = new amrex::MultiFab(ba, dm, static_cast<int>(nComp),
                                          static_cast<int>(nGhost));
         amrex::sundials::N_VSetOwnMF_MultiFab(v, SUNTRUE);
         amrex::sundials::getMFptr(v)     = mf_v;
    }

    return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a MultiFab N_Vector with user-specific MultiFab
 */

N_Vector N_VMake_MultiFab(sunindextype length, amrex::MultiFab *v_mf,
                          ::sundials::Context* sunctx)
{
    N_Vector v = N_VNewEmpty_MultiFab(length, sunctx);
    if (v == nullptr) return(nullptr);

    if (length > 0)
    {
         // Attach MultiFab
         amrex::sundials::N_VSetOwnMF_MultiFab(v, SUNFALSE);
         amrex::sundials::getMFptr(v)     = v_mf;
    }

    return(v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_MultiFab(N_Vector v)
{
    const auto* content = (amrex::sundials::N_VectorContent_MultiFab)(v->content);

    return content->length;
}

/* ----------------------------------------------------------------------------
 * Function to return if v owns the MultiFab*
 */
int N_VGetOwnMF_MultiFab(N_Vector v)
{
    const auto* content = (amrex::sundials::N_VectorContent_MultiFab)(v->content);

    return content->own_mf;
}

/* ----------------------------------------------------------------------------
 * Function to set if v owns the MultiFab*
 */
  void N_VSetOwnMF_MultiFab(N_Vector v, int own_mf_in)
{
    auto* content = (amrex::sundials::N_VectorContent_MultiFab)(v->content);

    content->own_mf = own_mf_in;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_MultiFab(N_Vector w)
{
    if (w == nullptr) return(nullptr);

    /* Create vector and copy operations */
    N_Vector v = N_VNewEmpty(w->sunctx);
    if (v == nullptr) return(nullptr);
    N_VCopyOps(w, v);

    /* Create content */
    auto* content = (N_VectorContent_MultiFab)
        std::malloc(sizeof(std::remove_pointer_t<N_VectorContent_MultiFab>));
    if (content == nullptr) { N_VFreeEmpty(v); return(nullptr); }

    content->length = amrex::sundials::N_VGetLength_MultiFab(w);
    content->own_mf = SUNFALSE;
    content->mf     = nullptr;

    /* Attach content */
    v->content = content;

    return(v);
}

N_Vector N_VClone_MultiFab(N_Vector w)
{
    N_Vector v = N_VCloneEmpty_MultiFab(w);
    if (v == nullptr) return(nullptr);

    sunindextype length = amrex::sundials::N_VGetLength_MultiFab(w);

    if (length > 0)
    {
         // Copy the multifab
         amrex::MultiFab *mf_w = amrex::sundials::getMFptr(w);
         const amrex::BoxArray &ba = mf_w->boxArray();
         const amrex::DistributionMapping &dm = mf_w->DistributionMap();
         int nComp = mf_w->nComp();
         int nGhost = mf_w->nGrow();  // same number of ghost cells in the clone
         auto *mf_v = new amrex::MultiFab(ba, dm, nComp, nGhost);

         // Attach multifab
         amrex::sundials::N_VSetOwnMF_MultiFab(v, SUNTRUE);
         amrex::sundials::getMFptr(v)     = mf_v;
    }

    return(v);
}

void N_VDestroy_MultiFab(N_Vector v)
{
    if (amrex::sundials::N_VGetOwnMF_MultiFab(v) == SUNTRUE)
    {
         delete amrex::sundials::getMFptr(v);
         amrex::sundials::getMFptr(v) = nullptr;
    }
    N_VFreeEmpty(v);
}

void N_VSpace_MultiFab(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
    *lrw = amrex::sundials::N_VGetLength_MultiFab(v);
    *liw = 1;
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
     return amrex::MultiFab(*((N_VectorContent_MultiFab)(v->content) )->mf,
                            amrex::make_alias, 0,
                            static_cast<int>((((N_VectorContent_MultiFab)(v->content) )->mf)->nComp()));
}

void N_VLinearSum_MultiFab(amrex::Real a, N_Vector x, amrex::Real b, N_Vector y,
                              N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::LinComb(*mf_z, a, *mf_x, 0, b, *mf_y, 0, 0, ncomp, nghost);
}

void N_VConst_MultiFab(amrex::Real c, N_Vector z)
{
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);
    *mf_z = c;
}

void N_VProd_MultiFab(N_Vector x, N_Vector y, N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
    amrex::MultiFab::Multiply(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VDiv_MultiFab(N_Vector x, N_Vector y, N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_y = amrex::sundials::getMFptr(y);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
    amrex::MultiFab::Divide(*mf_z, *mf_y, 0, 0, ncomp, nghost);
}

void N_VScale_MultiFab(amrex::Real c, N_Vector x, N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
    mf_z->mult(c, 0, ncomp, nghost);
}

void N_VAbs_MultiFab(N_Vector x, N_Vector z)
{
    using namespace amrex;

    MultiFab *mf_x = amrex::sundials::getMFptr(x);
    MultiFab *mf_z = amrex::sundials::getMFptr(z);
    int ncomp = mf_x->nComp();

    // ghost cells not included
    for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
    {
         const amrex::Box& bx = mfi.validbox();
         Array4<Real> const& x_fab = mf_x->array(mfi);
         Array4<Real> const& z_fab = mf_z->array(mfi);

         amrex::ParallelFor(bx, ncomp,
         [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
         {
             z_fab(i,j,k,c) = std::abs(x_fab(i,j,k,c));
         });
    }
}

void N_VInv_MultiFab(N_Vector x, N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
    mf_z->invert(1.0, 0, ncomp, nghost);
}

void N_VAddConst_MultiFab(N_Vector x, amrex::Real b, N_Vector z)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells

    amrex::MultiFab::Copy(*mf_z, *mf_x, 0, 0, ncomp, nghost);
    mf_z->plus(b, nghost);
}

amrex::Real N_VDotProd_MultiFab(N_Vector x, N_Vector y)
{
    using namespace amrex;

    MultiFab *mf_x = amrex::sundials::getMFptr(x);
    MultiFab *mf_y = amrex::sundials::getMFptr(y);
    int ncomp = mf_x->nComp();
    //int nghost = mf_x->nGrow();
    int nghost = 0;  // do not include ghost cells in dot product

    amrex::Real dotproduct = amrex::MultiFab::Dot(*mf_x, 0, *mf_y, 0, ncomp, nghost);

    return dotproduct;
}

amrex::Real N_VMaxNorm_MultiFab(N_Vector x)
{
    using namespace amrex;

    MultiFab *mf_x = amrex::sundials::getMFptr(x);
    int ncomp = mf_x->nComp();
    int startComp = 0;
    int nghost = 0;  // do not include ghost cells in the norm

    amrex::Real max = mf_x->max(startComp, nghost);

    // continue with rest of comps
    for (int c = 1; c < ncomp; ++c)
    {
         amrex::Real comp_max = mf_x->max(c, nghost); // comp c, no ghost zones
         if (comp_max > max)
         {
            max = comp_max;
         }
    }

    // no reduction needed, done in multifab
    return max;
}

amrex::Real N_VWrmsNorm_MultiFab(N_Vector x, N_Vector w)
{
    auto N = amrex::sundials::N_VGetLength_MultiFab(x);
    return N_VWL2Norm_MultiFab(x, w)*std::sqrt(1.0_rt/Real(N));
}

amrex::Real N_VWrmsNormMask_MultiFab(N_Vector x, N_Vector w, N_Vector id)
{
    return NormHelper_NVector_MultiFab(x, w, id, true, true);
}

amrex::Real N_VMin_MultiFab(N_Vector x)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    int ncomp = mf_x->nComp();

    int startComp = 0;
    int nghost = 0;  // ghost zones not included in min

    amrex::Real min = mf_x->min(startComp, nghost);

    // continue with rest of comps
    for (int c = 1; c < ncomp; ++c)
    {
         amrex::Real comp_min = mf_x->min(c, nghost); // comp c, no ghost zones
         if (comp_min < min)
         {
            min = comp_min;
         }
    }

    // no reduction needed, done in multifab
    return min;
}

amrex::Real NormHelper_NVector_MultiFab(N_Vector a_x, N_Vector a_w, N_Vector id, int use_id, bool rms)
{
    using namespace amrex;

    MultiFab *mf_x = amrex::sundials::getMFptr(a_x);
    //Call mf_y=mf_w is the MultiFab* from the N_Vector w
    MultiFab *mf_y = amrex::sundials::getMFptr(a_w);
    MultiFab *mf_id = use_id ? amrex::sundials::getMFptr(id) : nullptr;
    int numcomp = mf_x->nComp();
    sunindextype N = amrex::sundials::N_VGetLength_MultiFab(a_x);
    bool local = true;
    IntVect nghost = amrex::IntVect::TheZeroVector();
    Real sum = 0;
    int xcomp = 0;
    int ycomp = 0;

    // ghost cells not included
    if(use_id) {
        sum = amrex::NormHelper(*mf_id,
                                *mf_x, xcomp,
                                *mf_y, ycomp,
                                [=] AMREX_GPU_HOST_DEVICE (amrex::Real m) -> bool { return m > amrex::Real(0.0); },
                                [=] AMREX_GPU_HOST_DEVICE (amrex::Real x, amrex::Real y) -> amrex::Real { return x*x*y*y; },
                                numcomp, nghost, local);
    } else {
        sum = amrex::NormHelper(*mf_x, xcomp,
                                *mf_y, ycomp,
                                [=] AMREX_GPU_HOST_DEVICE (amrex::Real x, amrex::Real y) -> amrex::Real { return x*x*y*y; },
                                numcomp, nghost, local);
    }
    ParallelDescriptor::ReduceRealSum(sum);

    return rms ? std::sqrt(sum/Real(N)) : std::sqrt(sum);
}

amrex::Real N_VWL2Norm_MultiFab(N_Vector x, N_Vector w)
{
    using namespace amrex;

    return NormHelper_NVector_MultiFab(x, w, N_VCloneEmpty_MultiFab(x), false, false);
}

amrex::Real N_VL1Norm_MultiFab(N_Vector x)
{
    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    int ncomp = mf_x->nComp();

    int startComp = 0;
    int nghost = 0;  // ghost zones not included in norm

    amrex::Real sum = mf_x->norm1(startComp, nghost);

    // continue with rest of comps
    for (int c = 1; c < ncomp; ++c)
    {
         amrex::Real comp_sum = mf_x->norm1(c, nghost);
         sum += comp_sum;
    }

    // no reduction needed, it was done in multifab
    return sum;
}

void N_VCompare_MultiFab(amrex::Real a, N_Vector x, N_Vector z)
{
    using namespace amrex;

    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);
    int ncomp = mf_x->nComp();

    // ghost cells not included
    for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi)
    {
         const amrex::Box& bx = mfi.validbox();
         Array4<Real> const& x_fab = mf_x->array(mfi);
         Array4<Real> const& z_fab = mf_z->array(mfi);
         amrex::ParallelFor(bx, ncomp,
         [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
         {
           z_fab(i,j,k,c) = (std::abs(x_fab(i,j,k,c)) >= a) ? amrex::Real(1.0) : amrex::Real(0.0);
         });
    }
}

int N_VInvTest_MultiFab(N_Vector x, N_Vector z)
{
    using namespace amrex;

    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(x);
    amrex::MultiFab *mf_z = amrex::sundials::getMFptr(z);

    auto const& ma1 = mf_x->const_arrays();
    auto const& ma2 = mf_z->arrays();

    GpuTuple<bool> mm = ParReduce(TypeList<ReduceOpLogicalAnd>{},
                                    TypeList<bool>{},
                                    *mf_x,  amrex::IntVect::TheZeroVector(),
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
           -> GpuTuple<bool>
     {
         bool result = !(ma1[box_no](i,j,k) == amrex::Real(0.0));
         ma2[box_no](i,j,k) = result ? amrex::Real(1.0) / ma1[box_no](i,j,k) : 0.0;
         return { result };
     });

    bool val = amrex::get<0>(mm);

    if (val == false)
    {
         return SUNFALSE;
    }
    else
    {
         return SUNTRUE;
    }
}

int N_VConstrMask_MultiFab(N_Vector a_a, N_Vector a_x, N_Vector a_m)
{
    using namespace amrex;

    amrex::MultiFab *mf_x = amrex::sundials::getMFptr(a_x);
    amrex::MultiFab *mf_a = amrex::sundials::getMFptr(a_a);
    amrex::MultiFab *mf_m = amrex::sundials::getMFptr(a_m);
    int ncomp = mf_x->nComp();

    // ghost cells not included
    auto temp = amrex::Real(0.0);
    for (MFIter mfi(*mf_x); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        Array4<Real> const& x_fab = mf_x->array(mfi);
        Array4<Real> const& a_fab = mf_a->array(mfi);
        Array4<Real> const& m_fab = mf_a->array(mfi);
        //Changing continue to if check, temp calculation should be changed to reduction
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
            {
                m_fab(i,j,k,c) = amrex::Real(0.0);

                /* Continue if no constraints were set for the variable */
                if (a_fab(i,j,k,c) != amrex::Real(0.0)) {

                    /* Check if a set constraint has been violated */
                    amrex::Real a, x;
                    a = a_fab(i,j,k,c);
                    x = x_fab(i,j,k,c);
                    int test = (std::abs(a) > amrex::Real(1.5) && x*a <= amrex::Real(0.0)) ||
                        (std::abs(a) > amrex::Real(0.5)   && x*a <  amrex::Real(0.0));
                    if (test) {
                      m_fab(i,j,k,c) = amrex::Real(1.0);
                    }
                }
            });
    }

    temp = mf_m->norm1();
    /* Return false if any constraint was violated */
    amrex::ParallelDescriptor::ReduceRealMax(temp);

    return (temp == amrex::Real(1.0)) ? SUNFALSE : SUNTRUE;
}

amrex::Real N_VMinQuotient_MultiFab(N_Vector a_num, N_Vector a_denom)
{
    using namespace amrex;

    amrex::MultiFab *mf_num = amrex::sundials::getMFptr(a_num);
    amrex::MultiFab *mf_denom = amrex::sundials::getMFptr(a_denom);
    int ncomp = mf_num->nComp();

    // ghost cells not included
    int nghost = 0;
    Real min = amrex::ReduceMin(*mf_num, *mf_denom, nghost,
                 [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& num_fab, Array4<Real const> const& denom_fab) -> Real
    {
         Real min_loc = std::numeric_limits<Real>::max();
         const auto lo = lbound(bx);
         const auto hi = ubound(bx);

         for (int c = 0; c < ncomp; ++c) {
            for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
                  for (int i = lo.x; i <= hi.x; ++i) {
                     if (denom_fab(i,j,k,c) != amrex::Real(0.0))
                     {
                        amrex::Real num = num_fab(i,j,k,c);
                        amrex::Real denom = denom_fab(i,j,k,c);
                        min_loc = std::min(min_loc, num / denom);
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
