#ifndef AMREX_HABEC_2D_H_
#define AMREX_HABEC_2D_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_mat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, int i, int j, int k,
                Dim3 const& boxlo, Dim3 const& boxhi,
                Real sa, Array4<Real const> const& a,
                Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int bho,
                GpuArray<Array4<int const>, AMREX_SPACEDIM*2> const& msk,
                Array4<Real> const& diaginv)
{
    sten[1] = -(sb / (dx[0]*dx[0])) * b[0](i,j,k);
    sten[2] = -(sb / (dx[0]*dx[0])) * b[0](i+1,j,k);
    sten[0] = -(sten[1] + sten[2]);
    if (sa != Real(0.0)) {
        sten[0] += sa * a(i,j,k);
    }

    // xlo
    if (i == boxlo.x) {
        int cdir = Orientation(Direction::x, Orientation::low);
        if (msk[cdir](i-1,j,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
        }
    }

    // xhi
    if (i == boxhi.x) {
        int cdir = Orientation(Direction::x, Orientation::high);
        if (msk[cdir](i+1,j,k) > 0) {
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
        }
    }
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_ijmat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, Array4<Int> const& ncols,
                  Array4<Real> const& diaginv, int i, int j, int k,
                  Array4<Int const> const& cell_id,
                  Real sa, Array4<Real const> const& a,
                  Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                  GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                  GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int bho,
                  Array4<int const> const& osm)
{
    if (!osm || osm(i,j,k) != 0) {
        sten[1] = -(sb / (dx[0]*dx[0])) * b[0](i,j,k);
        sten[2] = -(sb / (dx[0]*dx[0])) * b[0](i+1,j,k);
        sten[0] = -(sten[1] + sten[2]);
        if (sa != Real(0.0)) {
            sten[0] += sa * a(i,j,k);
        }

        // xlo
        if (cell_id(i-1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::low);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
        }

        // xhi
        if (cell_id(i+1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::high);
            Real bf1, bf2;
            detail::comp_bf(bf1, bf2, sb, dx[0], bctype[cdir], bcl[cdir], bho);
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
        }
    } else {
        sten[0] = Real(1.0);
        for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
            sten[m] = Real(0.0);
        }
    }

    diaginv(i,j,k) = Real(1.0) / sten[0];
    sten[0] = Real(1.0);
    for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
        sten[m] *= diaginv(i,j,k);
    }

    ncols(i,j,k) = 0;
    for (int m = 0; m < 2*AMREX_SPACEDIM+1; ++m) {
        ncols(i,j,k) += (sten[m] != Real(0.0));
    }
}

template <typename Int>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_cols (GpuArray<Int,2*AMREX_SPACEDIM+1>& sten, int i, int j, int /*k*/,
                 Array4<Int const> const& cell_id)
{
    sten[0] = cell_id(i  ,j  ,0);
    sten[1] = cell_id(i-1,j  ,0);
    sten[2] = cell_id(i+1,j  ,0);
}

}
#endif
