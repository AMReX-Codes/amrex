#ifndef AMREX_FFT_H_
#define AMREX_FFT_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_FFT_LocalR2C.H>
#include <AMReX_FFT_OpenBCSolver.H>
#include <AMReX_FFT_R2C.H>
#include <AMReX_FFT_R2X.H>

namespace amrex
{
    /**
     * \brief Initialize FFT
     *
     * This is needed only when the user wants to use amrex::FFT, but does
     * not want to call amrex::Initialize to initialize the full version of
     * AMReX. Note that one usually only needs to call Init_FFT and
     * Finalize_FFT once in the entire program.
     */
#ifdef AMREX_USE_MPI
    inline void Init_FFT (MPI_Comm comm ) { amrex::Init_minimal(comm); }
#else
    inline void Init_FFT () { amrex::Init_minimal(); }
#endif

    //! If Init_FFT is called, this should be called after all the FFT works
    //! are done.
    inline void Finalize_FFT () { amrex::Finalize_minimal(); }
}

namespace amrex::FFT
{
    void Initialize ();
    void Finalize ();
    void Clear ();
}

#endif
