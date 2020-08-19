#ifndef _Adv_K_H_
#define _Adv_K_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void conservative(int i, int j, int k, 
                  amrex::Array4<amrex::Real> const& phi_out,
                  amrex::Array4<amrex::Real const> const& phi_in,
                  AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                               amrex::Array4<amrex::Real const> const& flxy,
                               amrex::Array4<amrex::Real const> const& flxz),
                  AMREX_D_DECL(amrex::Real dtdx, amrex::Real dtdy, amrex::Real dtdz))
{
    phi_out(i,j,k) = phi_in(i,j,k) + 
                ( AMREX_D_TERM( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx,
                              + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdy,
                              + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdz ) );
}

#endif
