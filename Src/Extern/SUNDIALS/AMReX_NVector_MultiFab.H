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

/**
 * \file AMReX_NVector_MultiFab.H
 *
 * SUNDIALS N_Vector implementation that stores data inside AMReX MultiFabs,
 * plus the associated helper routines.
 */

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
  sunindextype length;   //!< Logical vector length (number of scalar entries).
  int own_mf;            //!< Non-zero if this wrapper owns the MultiFab.
  amrex::MultiFab *mf;   //!< Wrapped MultiFab storage.
};

using N_VectorContent_MultiFab = N_VectorContent_MultiFab_notptr*;

/* -----------------------------------------------------------------
 * Exported functions
 * -----------------------------------------------------------------*/

/**
 * \brief Create an empty MultiFab N_Vector shell of length \p vec_length.
 *
 * The returned vector owns no data until a MultiFab is attached.
 *
 * \param vec_length Number of scalar entries represented.
 * \param sunctx     Optional SUNDIALS context; defaults to amrex::sundials::The_Sundials_Context().
 * \return Newly allocated N_Vector shell.
 */
N_Vector N_VNewEmpty_MultiFab(sunindextype vec_length,
                              ::sundials::Context* sunctx = The_Sundials_Context());

/**
 * \brief Allocate a MultiFab-backed N_Vector of length \p vec_length.
 *
 * \param vec_length Number of scalar entries represented.
 * \param ba         Grid layout.
 * \param dm         Distribution map for \p ba.
 * \param nComp      Number of state components per cell.
 * \param nGhost     Number of grow cells.
 * \param sunctx     Optional SUNDIALS context used for bookkeeping.
 * \return Newly allocated N_Vector that owns its MultiFab data.
 */
N_Vector N_VNew_MultiFab(sunindextype vec_length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost,
                         ::sundials::Context* sunctx = The_Sundials_Context());

/**
 * \brief Wrap an existing MultiFab \p mf as an N_Vector without copying.
 *
 * \param vec_length Number of scalar entries represented.
 * \param mf         MultiFab that supplies the storage.
 * \param sunctx     Optional SUNDIALS context; defaults to The_Sundials_Context().
 * \return Non-owning N_Vector referencing \p mf.
 */
N_Vector N_VMake_MultiFab(sunindextype vec_length,
                          amrex::MultiFab *mf,
                          ::sundials::Context* sunctx = The_Sundials_Context());

/**
 * \brief Return the logical length of N_Vector \p v.
 *
 * \param v Vector to inspect.
 * \return Number of scalar entries represented by \p v.
 */
sunindextype N_VGetLength_MultiFab(N_Vector v);

/**
 * \brief Query whether \p v owns its MultiFab storage.
 *
 * \param v Vector to inspect.
 * \return Non-zero if the vector owns its MultiFab, zero otherwise.
 */
int N_VGetOwnMF_MultiFab(N_Vector v);

/**
 * \brief Set the ownership flag on \p v's MultiFab.
 *
 * \param v      Vector whose flag is modified.
 * \param own_mf Non-zero if the vector should free its MultiFab on destroy.
 */
void N_VSetOwnMF_MultiFab(N_Vector v, int own_mf);

/**
 * \brief Create an empty clone of \p w that owns no data.
 *
 * \param w Template vector describing the layout.
 * \return Empty clone matching \p w.
 */
N_Vector N_VCloneEmpty_MultiFab(N_Vector w);

/**
 * \brief Create a deep clone of \p w including its MultiFab data.
 *
 * \param w Template vector to clone.
 * \return New vector with copied storage.
 */
N_Vector N_VClone_MultiFab(N_Vector w);

/**
 * \brief Destroy \p v and free its resources (MultiFab if owned).
 *
 * \param v Vector to destroy.
 */
void N_VDestroy_MultiFab(N_Vector v);

/**
 * \brief Return workspace requirements for \p v.
 *
 * \param v   Vector to inspect.
 * \param lrw Set to the number of real workspace entries.
 * \param liw Set to the number of integer workspace entries.
 */
void N_VSpace_MultiFab(N_Vector v, sunindextype *lrw,
                       sunindextype *liw);

/**
 * \brief Access the implementation-specific content stored in \p v.
 *
 * \param v Vector whose content is returned.
 * \return Pointer to the MultiFab content structure.
 */
N_VectorContent_MultiFab N_VGetContent_MultiFab(N_Vector v);

/**
 * \brief Return the wrapped MultiFab pointer from \p v.
 *
 * \param v Vector whose MultiFab pointer is returned.
 * \return Pointer to the underlying MultiFab.
 */
MultiFab* N_VGetVectorPointer_MultiFab(N_Vector v);

/**
 * \brief Helper that computes weighted norms or RMS values.
 *
 * \param x      Input vector.
 * \param w      Weight vector (optional).
 * \param id     Mask vector (optional).
 * \param use_id Non-zero to restrict to entries where \p id is non-zero.
 * \param rms    True for RMS norm, false for weighted L2.
 * \return Computed norm value.
 */
amrex::Real NormHelper_NVector_MultiFab(N_Vector x, N_Vector w, N_Vector id, int use_id, bool rms);

/* standard vector operations */

/**
 * \brief Compute `z = a x + b y`.
 *
 * \param a Scalar multiplier for \p x.
 * \param x First operand.
 * \param b Scalar multiplier for \p y.
 * \param y Second operand.
 * \param z Output vector receiving the sum.
 */
void N_VLinearSum_MultiFab(amrex::Real a, N_Vector x,
                           amrex::Real b, N_Vector y, N_Vector z);

/**
 * \brief Fill \p z with the constant \p c.
 *
 * \param c Value assigned to every entry.
 * \param z Destination vector.
 */
void N_VConst_MultiFab(amrex::Real c, N_Vector z);

/**
 * \brief Compute the element-wise product `z = x * y`.
 *
 * \param x First operand.
 * \param y Second operand.
 * \param z Output vector receiving the product.
 */
void N_VProd_MultiFab(N_Vector x, N_Vector y, N_Vector z);

/**
 * \brief Compute the element-wise quotient `z = x / y`.
 *
 * \param x Numerator vector.
 * \param y Denominator vector.
 * \param z Output vector receiving the quotient.
 */
void N_VDiv_MultiFab(N_Vector x, N_Vector y, N_Vector z);

/**
 * \brief Compute the scalar product `z = c * x`.
 *
 * \param c Scalar multiplier.
 * \param x Input vector.
 * \param z Output vector receiving the scaled values.
 */
void N_VScale_MultiFab(amrex::Real c, N_Vector x, N_Vector z);

/**
 * \brief Compute the element-wise absolute value `z = |x|`.
 *
 * \param x Input vector.
 * \param z Output vector receiving the magnitudes.
 */
void N_VAbs_MultiFab(N_Vector x, N_Vector z);

/**
 * \brief Compute the element-wise inverse `z = 1 / x`.
 *
 * \param x Input vector.
 * \param z Output vector receiving the reciprocals.
 */
void N_VInv_MultiFab(N_Vector x, N_Vector z);

/**
 * \brief Compute `z = x + b`.
 *
 * \param x Input vector.
 * \param b Scalar shift applied to every entry.
 * \param z Output vector.
 */
void N_VAddConst_MultiFab(N_Vector x, amrex::Real b, N_Vector z);

/**
 * \brief Return the dot product `x · y`.
 *
 * \param x First operand.
 * \param y Second operand.
 * \return Scalar dot product.
 */
amrex::Real N_VDotProd_MultiFab(N_Vector x, N_Vector y);

/**
 * \brief Return the infinity norm of \p x.
 *
 * \param x Vector whose norm is computed.
 * \return Infinity norm.
 */
amrex::Real N_VMaxNorm_MultiFab(N_Vector x);

/**
 * \brief Return the weighted RMS norm of \p x using weights \p w.
 *
 * \param x Vector whose norm is computed.
 * \param w Weight vector.
 * \return Weighted RMS norm.
 */
amrex::Real N_VWrmsNorm_MultiFab(N_Vector x, N_Vector w);

/**
 * \brief Return the weighted RMS norm restricted by mask \p id.
 *
 * \param x  Vector whose norm is computed.
 * \param w  Weight vector.
 * \param id Mask vector indicating included entries.
 * \return Masked weighted RMS norm.
 */
amrex::Real N_VWrmsNormMask_MultiFab(N_Vector x, N_Vector w,
                                     N_Vector id);

/**
 * \brief Return the minimum entry of \p x.
 *
 * \param x Vector to inspect.
 * \return Minimum value.
 */
amrex::Real N_VMin_MultiFab(N_Vector x);

/**
 * \brief Return the weighted L2 norm of \p x using weights \p w.
 *
 * \param x Vector whose norm is computed.
 * \param w Weight vector.
 * \return Weighted L2 norm.
 */
amrex::Real N_VWL2Norm_MultiFab(N_Vector x, N_Vector w);

/**
 * \brief Return the L1 norm (sum of absolute values) of \p x.
 *
 * \param x Vector whose norm is computed.
 * \return L1 norm.
 */
amrex::Real N_VL1Norm_MultiFab(N_Vector x);

/**
 * \brief Build a mask \p z where entries of \p x exceed \p a.
 *
 * \param a Threshold value.
 * \param x Input vector.
 * \param z Output mask (1 where |x| >= a, 0 otherwise).
 */
void N_VCompare_MultiFab(amrex::Real a, N_Vector x, N_Vector z);

/**
 * \brief Compute the element-wise inverse test `z = 1 / x`, returning failure if zero encountered.
 *
 * \param x Input vector.
 * \param z Output vector receiving reciprocals.
 * \return 1 on success, 0 if a zero pivot was encountered.
 */
int N_VInvTest_MultiFab(N_Vector x, N_Vector z);

/**
 * \brief Apply constraint masking from \p a_a to produce mask \p m.
 *
 * \param a_a Constraint vector (per SUNDIALS conventions).
 * \param x   Solution vector tested against constraints.
 * \param m   Output mask updated in place.
 * \return Non-zero if all constraints satisfied, zero otherwise.
 */
int N_VConstrMask_MultiFab(N_Vector a_a, N_Vector x,
                           N_Vector m);

/**
 * \brief Return the minimum element-wise quotient `num / denom`.
 *
 * \param num   Numerator vector.
 * \param denom Denominator vector.
 * \return Minimum quotient over all entries.
 */
amrex::Real N_VMinQuotient_MultiFab(N_Vector num, N_Vector denom);

#ifdef __cplusplus
} // extern "C"

/**
 * \brief Access the MultiFab pointer stored inside \p v (non-const).
 *
 * \param v Vector whose pointer is returned.
 * \return Reference to the MultiFab pointer for direct manipulation.
 */
MultiFab*& getMFptr(N_Vector v);

/**
 * \brief Return a shallow alias of the MultiFab stored in \p v.
 *
 * \param v Vector whose data are aliased.
 * \return MultiFab alias referencing the same storage as \p v.
 */
MultiFab N_VGetVectorAlias_MultiFab(N_Vector v);

}//namespace amrex::sundials
#endif

#endif
