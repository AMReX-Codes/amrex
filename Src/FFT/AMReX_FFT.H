#ifndef AMREX_FFT_H_
#define AMREX_FFT_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_FFT_LocalR2C.H>
#include <AMReX_FFT_OpenBCSolver.H>
#include <AMReX_FFT_R2C.H>
#include <AMReX_FFT_R2X.H>

/**
 * \file AMReX_FFT.H
 *
 * Thin convenience layer that exposes FFT initialization/finalization helpers.
 */

namespace amrex
{
#ifdef AMREX_USE_MPI
    /**
     * \brief Initialize FFT helpers without invoking the full AMReX runtime.
     *
     * This is needed only when the user wants to use amrex::FFT, but does
     * not want to call amrex::Initialize to initialize the full version of
     * AMReX. Note that one usually only needs to call Init_FFT and
     * Finalize_FFT once in the entire program.
     *
     * \param comm MPI communicator passed to amrex::Init_minimal.
     */
    inline void Init_FFT (MPI_Comm comm ) { amrex::Init_minimal(comm); }
#else
    /**
     * \brief Initialize FFT helpers without invoking the full AMReX runtime.
     *
     * This overload is used in serial builds and has the same semantics as
     * the MPI version minus the communicator argument.
     */
    inline void Init_FFT () { amrex::Init_minimal(); }
#endif

    /**
     * \brief Tear down the minimal FFT runtime previously started via Init_FFT().
     */
    inline void Finalize_FFT () { amrex::Finalize_minimal(); }
}

/// \cond DOXYGEN_IGNORE
//! Not to be called by the user
namespace amrex::FFT::detail
{
    void Initialize ();
    void Finalize ();
    void Clear ();
}
/// \endcond

#endif
