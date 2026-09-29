#ifndef AMREX_HABEC_1D_H_
#define AMREX_HABEC_1D_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_mat (GpuArray<Real,2*AMREX_SPACEDIM+1>& sten, int i, int j, int k,
                Dim3 const& boxlo, Dim3 const& boxhi,
                Real sa, Array4<Real const> const& a,
                Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int maxorder,
                GpuArray<Array4<int const>, AMREX_SPACEDIM*2> const& msk)
{
    maxorder = amrex::min(maxorder,3); // the structured stencil has room for order 3
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
            Real bf1, bf2, bf3;
            detail::comp_bf(bf1, bf2, bf3, sb, dx[0], bctype[cdir], bcl[cdir],
                            amrex::min(boxhi.x-boxlo.x+2, maxorder));
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
        }
    }

    // xhi
    if (i == boxhi.x) {
        int cdir = Orientation(Direction::x, Orientation::high);
        if (msk[cdir](i+1,j,k) > 0) {
            Real bf1, bf2, bf3;
            detail::comp_bf(bf1, bf2, bf3, sb, dx[0], bctype[cdir], bcl[cdir],
                            amrex::min(boxhi.x-boxlo.x+2, maxorder));
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
        }
    }
}

// NS is 2*AMREX_SPACEDIM+1, or 4*AMREX_SPACEDIM+1 for maxorder 4, whose extra
// entries couple to the cells two away (-x, +x, -y, +y, -z, +z).
template <typename AlgInt, unsigned int NS>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_ijmat (GpuArray<Real,NS>& sten, Array4<AlgInt> const& ncols,
                  int i, int j, int k, Box const& vbx,
                  Array4<AlgInt const> const& cell_id,
                  Real sa, Array4<Real const> const& a,
                  Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                  GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                  GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int maxorder,
                  Array4<int const> const& osm)
{
    constexpr unsigned int NR = 2*AMREX_SPACEDIM+1;
    static_assert(NS == NR || NS == 2*NR-1);
    int const mo = (NS > NR) ? maxorder : amrex::min(maxorder,3);
    for (unsigned int m = NR; m < NS; ++m) {
        sten[m] = Real(0.0);
    }

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
            Real bf1, bf2, bf3;
            detail::comp_bf(bf1, bf2, bf3, sb, dx[0], bctype[cdir], bcl[cdir],
                            amrex::min(vbx.length(0)+1, mo));
            sten[0] += bf1 * b[0](i,j,k);
            sten[1] = Real(0.0);
            sten[2] += bf2 * b[0](i,j,k);
            if constexpr (NS > NR) { sten[NR+1] += bf3 * b[0](i,j,k); }
        }

        // xhi
        if (cell_id(i+1,j,k) < 0) {
            int cdir = Orientation(Direction::x, Orientation::high);
            Real bf1, bf2, bf3;
            detail::comp_bf(bf1, bf2, bf3, sb, dx[0], bctype[cdir], bcl[cdir],
                            amrex::min(vbx.length(0)+1, mo));
            sten[0] += bf1 * b[0](i+1,j,k);
            sten[1] += bf2 * b[0](i+1,j,k);
            sten[2] = Real(0.0);
            if constexpr (NS > NR) { sten[NR+0] += bf3 * b[0](i+1,j,k); }
        }
    } else {
        sten[0] = Real(1.0);
        for (unsigned int m = 1; m < NS; ++m) {
            sten[m] = Real(0.0);
        }
    }

    ncols(i,j,k) = 0;
    for (unsigned int m = 0; m < NS; ++m) {
        ncols(i,j,k) += (sten[m] != Real(0.0));
    }
}

template <typename AlgInt, unsigned int NS>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void habec_cols (GpuArray<AlgInt,NS>& sten, int i, int j, int /*k*/, Box const& vbx,
                 Array4<AlgInt const> const& cell_id)
{
    sten[0] = cell_id(i  ,j  ,0);
    sten[1] = cell_id(i-1,j  ,0);
    sten[2] = cell_id(i+1,j  ,0);
    constexpr unsigned int NR = 2*AMREX_SPACEDIM+1;
    if constexpr (NS > NR) {
        sten[NR+0] = (i-2 >= vbx.smallEnd(0)) ? cell_id(i-2,j,0) : AlgInt(-1);
        sten[NR+1] = (i+2 <= vbx.bigEnd(0)) ? cell_id(i+2,j,0) : AlgInt(-1);
    }
}

}
#endif
