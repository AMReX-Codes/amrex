#ifndef AMREX_FFT_H_
#define AMREX_FFT_H_
#include <AMReX_Config.H>

#include <AMReX_FFT_LocalR2C.H>
#include <AMReX_FFT_OpenBCSolver.H>
#include <AMReX_FFT_R2C.H>
#include <AMReX_FFT_R2X.H>

namespace amrex::FFT
{
    void Initialize ();
    void Finalize ();
    void Clear ();
}

#endif
