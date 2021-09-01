#ifndef Adv_K_H_
#define Adv_K_H_

#include <AMReX_Array4.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void conservative(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& phi_out,
                  amrex::Array4<amrex::Real const> const& phi_in,
                  AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                               amrex::Array4<amrex::Real const> const& flxy,
                               amrex::Array4<amrex::Real const> const& flxz),
                  AMREX_D_DECL(const amrex::Real dtdx,
                               const amrex::Real dtdy,
                               const amrex::Real dtdz))
{
    phi_out(i,j,k) = phi_in(i,j,k) +
                ( AMREX_D_TERM( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx,
                              + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdy,
                              + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdz ) );
}

#endif /* Adv_K_H_ */
