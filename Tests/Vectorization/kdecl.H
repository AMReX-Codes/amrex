#ifndef KDECL_H_
#define KDECL_H_

#include <AMReX_FArrayBox.H>

void ctoprim_c_simd
    (amrex::Box const& bx, amrex::FArrayBox const& ufab, amrex::FArrayBox & qfab);

void ctoprim_c_nosimd
    (amrex::Box const& bx, amrex::FArrayBox const& ufab, amrex::FArrayBox & qfab);

void flux_to_dudt_c_simd
    (amrex::Box const& bx,
     amrex::FArrayBox& dudtfab,
     amrex::FArrayBox const& fxfab,
     amrex::FArrayBox const& fyfab,
     amrex::FArrayBox const& fzfab,
     amrex::Array<amrex::Real,AMREX_SPACEDIM> const& dxinv);

void flux_to_dudt_c_nosimd
    (amrex::Box const& bx,
     amrex::FArrayBox& dudtfab,
     amrex::FArrayBox const& fxfab,
     amrex::FArrayBox const& fyfab,
     amrex::FArrayBox const& fzfab,
     amrex::Array<amrex::Real,AMREX_SPACEDIM> const& dxinv);

extern "C" {
    void ctoprim_f (const int* lo, const int* hi,
                    amrex::Real* u, const int* ulo, const int* uhi,
                    const amrex::Real* q, const int* qlo, const int* qhi);

    void flux_to_dudt_f (const int* lo, const int* hi,
                         amrex::Real* dudt, const int* ulo, const int* uhi,
                         const amrex::Real* fx, const int* xlo, const int* xhi,
                         const amrex::Real* fy, const int* ylo, const int* yhi,
                         const amrex::Real* fz, const int* zlo, const int* zhi,
                         const amrex::Real* dxinv, const int* ncomp);
}

#endif
