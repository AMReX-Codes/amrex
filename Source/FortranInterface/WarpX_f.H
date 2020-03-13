/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Jean-Luc Vay, Ligia Diana Amorim
 * Luca Fedeli, Mathieu Lobet, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_F_H_
#define WARPX_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus

#if AMREX_SPACEDIM==2
#define WARPX_ARLIM_ANYD(x) std::array<int,3>{(x)[0], 0, (x)[1]}.data()
#else
#define WARPX_ARLIM_ANYD(x) x
#endif

#define WARPX_TO_FORTRAN_BOX(x) WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_ANYD(x) (x).dataPtr(), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())

#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef WARPX_USE_PSATD
    void warpx_fft_mpi_init (int fcomm);
    void warpx_fft_domain_decomp (int* warpx_local_nz, int* warpx_local_z0,
                                  const int* global_lo, const int* global_hi);
    void warpx_fft_dataplan_init (const int* nox, const int* noy, const int* noz,
                                  void* fft_data, const int* ndata,
                                  const amrex_real* dx_w, const amrex_real* dt_w,
                                  const int* fftw_plan_measure, const int* do_nodal );
    void warpx_fft_nullify ();
    void warpx_fft_push_eb (amrex_real* ex_w, const int* exlo, const int* exhi,
                            amrex_real* ey_w, const int* eylo, const int* eyhi,
                            amrex_real* ez_w, const int* ezlo, const int* ezhi,
                            amrex_real* bx_w, const int* bxlo, const int* bxhi,
                            amrex_real* by_w, const int* bylo, const int* byhi,
                            amrex_real* bz_w, const int* bzlo, const int* bzhi,
                            amrex_real* jx_w, const int* jxlo, const int* jxhi,
                            amrex_real* jy_w, const int* jylo, const int* jyhi,
                            amrex_real* jz_w, const int* jzlo, const int* jzhi,
                            amrex_real* rhoold_w, const int* r1lo, const int* r1hi,
                            amrex_real* rho_w, const int* r2lo, const int* r2hi);

#endif

#ifdef __cplusplus
}
#endif

#endif //WARPX_F_H_
