#ifndef AMREX_EB_FLUXREGISTER_3D_C_H_
#define AMREX_EB_FLUXREGISTER_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_crseadd_va(int i, int j, int k, int n, Array4<Real> const& d,
                            Array4<int const> const& flag, Array4<Real const> const& fx,
                            Array4<Real const> const& fy, Array4<Real const> const& fz,
                            Array4<Real const> const& vfrac, Array4<Real const> const& ax,
                            Array4<Real const> const& ay, Array4<Real const> const& az,
                            Real dtdx, Real dtdy, Real dtdz)
{
    if (flag(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell
        and vfrac(i,j,k) > 1.e-14)
    {
        Real volinv = 1.0/vfrac(i,j,k);

        if (flag(i-1,j,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) -= dtdx*fx(i,j,k,n)*(ax(i,j,k)*volinv);
        } else if (flag(i+1,j,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) +=  dtdx*fx(i+1,j,k,n)*(ax(i+1,j,k)*volinv);
        }

        if (flag(i,j-1,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) -= dtdy*fy(i,j,k,n)*(ay(i,j,k)*volinv);
        } else if (flag(i,j+1,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) += dtdy*fy(i,j+1,k,n)*(ay(i,j+1,k)*volinv);
        }

        if (flag(i,j,k-1) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) -= dtdz*fz(i,j,k,n)*(az(i,j,k)*volinv);
        } else if (flag(i,j,k+1) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) += dtdz*fz(i,j,k+1,n)*(az(i,j,k+1)*volinv);
        }
    }
}

}

#endif
