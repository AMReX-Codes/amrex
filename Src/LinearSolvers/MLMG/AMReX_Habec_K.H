#ifndef AMREX_Habec_K_H_
#define AMREX_Habec_K_H_
#include <AMReX_Config.H>

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LOUtil_K.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_BaseFabUtility.H>

/// \cond DOXYGEN_IGNORE
namespace amrex::detail {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void comp_bf (Real& bf1, Real& bf2, Real sb, Real h, int bct, Real bcl, int bho)
    {
        Real fac = sb / (h*h);
        if (bct == AMREX_LO_DIRICHLET) {
            Real h2 = Real(0.5)*h;
            if (bho >= 1) {
                Real h3 = Real(3.0)*h2;
                bf1 = fac*((h3 - bcl)/(bcl + h2) - Real(1.0));
                bf2 = fac*(bcl - h2)/(bcl + h3);
            } else {
                bf1 = fac*( h/(bcl + h2) - Real(1.0));
                bf2 = Real(0.0);
            }
        } else if (bct == AMREX_LO_NEUMANN) {
            bf1 = -fac;
            bf2 = Real(0.0);
        } else if (bct == AMREX_LO_REFLECT_ODD) {
            bf1 = fac; // ghost = -interior
            bf2 = Real(0.0);
        } else {
            bf1 = bf2 = Real(0.0);
            AMREX_ASSERT_WITH_MESSAGE(false,"hpmat: unsupported boundary type");
        }
    }

    // Boundary stencil of the given order (2 to 4); bf3 couples to the third cell.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void comp_bf (Real& bf1, Real& bf2, Real& bf3, Real sb, Real h, int bct, Real bcl,
                  int order)
    {
        if (order == 4 && bct == AMREX_LO_DIRICHLET) {
            Real const fac = sb / (h*h);
            Real x[4] = {-bcl/h, Real(0.5), Real(1.5), Real(2.5)};
            Real c[4];
            poly_interp_coeff(-Real(0.5), x, 4, c);
            bf1 = -fac*c[1];
            bf2 = -fac*c[2];
            bf3 = -fac*c[3];
        } else {
            comp_bf(bf1, bf2, sb, h, bct, bcl, amrex::min(order,3)-2);
            bf3 = Real(0.0);
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void comp_bflo (Real& bf1, Real& bf2, Real& bflo, Real sb, Real h, int bct, Real bcl, int bho)
    {
        Real fac = sb / (h*h);
        if (bct == AMREX_LO_DIRICHLET) {
            Real h2 = Real(0.5)*h;
            bflo = fac * ( h / (bcl + h2) - Real(1.0));
            if (bho >= 1) {
                Real h3 = Real(3.0)*h2;
                bf1 = fac*((h3 - bcl)/(bcl + h2) - Real(1.0));
                bf2 = fac*(bcl - h2)/(bcl + h3);
            } else {
                bf1 = bflo;
                bf2 = Real(0.0);
            }
        } else if (bct == AMREX_LO_NEUMANN) {
            bflo = -fac;
            bf1 = -fac;
            bf2 = Real(0.0);
        } else if (bct == AMREX_LO_REFLECT_ODD) {
            bflo = fac; // ghost = -interior
            bf1 = fac;
            bf2 = Real(0.0);
        } else {
            bf1 = bf2 = bflo = Real(0.0);
            AMREX_ASSERT_WITH_MESSAGE(false,"hpmat: unsupported boundary type");
        }
    }
}
/// \endcond

#if (AMREX_SPACEDIM == 1)
#include <AMReX_Habec_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_Habec_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_Habec_3D_K.H>
#endif

namespace amrex {

/**
 * \brief Fill the matrix and column arrays of the regular IJ stencil on box \p bx.
 *
 * \p mat and \p cols hold NS entries per cell. With \p drop_known, the
 * couplings into cells whose overset mask \p osm is 0 are removed; \p osm
 * then needs one ghost cell.
 */
template <unsigned int NS, typename AlgInt>
void habec_ij_fill (Box const& bx,
                    Real* mat, // NOLINT(readability-non-const-parameter)
                    AlgInt* cols, Array4<AlgInt> const& ncols,
                    Array4<AlgInt const> const& cell_id,
                    Real sa, Array4<Real const> const& a,
                    Real sb, GpuArray<Real,AMREX_SPACEDIM> const& dx,
                    GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& b,
                    GpuArray<int,AMREX_SPACEDIM*2> const& bctype,
                    GpuArray<Real,AMREX_SPACEDIM*2> const& bcl, int maxorder,
                    Array4<int const> const& osm, bool drop_known)
{
    BaseFab<GpuArray<Real,NS>> matfab(bx, 1, (GpuArray<Real,NS>*)mat);
    amrex::fill(matfab,
    [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,NS>& sten, int i, int j, int k)
    {
        habec_ijmat(sten, ncols, i, j, k, bx, cell_id, sa, a, sb, dx, b,
                    bctype, bcl, maxorder, osm);
        if (drop_known && osm && osm(i,j,k) != 0) {
            AMREX_D_TERM(if (osm(i-1,j,k) == 0) { sten[1] = Real(0.0); }
                         if (osm(i+1,j,k) == 0) { sten[2] = Real(0.0); },
                         if (osm(i,j-1,k) == 0) { sten[3] = Real(0.0); }
                         if (osm(i,j+1,k) == 0) { sten[4] = Real(0.0); },
                         if (osm(i,j,k-1) == 0) { sten[5] = Real(0.0); }
                         if (osm(i,j,k+1) == 0) { sten[6] = Real(0.0); })
            constexpr unsigned int NR = 2*AMREX_SPACEDIM+1;
            if constexpr (NS > NR) { // nonzero only inside bx, away from the boundary
                AMREX_D_TERM(int const sx = (cell_id(i-1,j,k) < 0) ? 2 : -2;
                             if (sten[NR+0] != Real(0.0) && osm(i+sx,j,k) == 0) { sten[NR+0] = Real(0.0); },
                             int const sy = (cell_id(i,j-1,k) < 0) ? 2 : -2;
                             if (sten[NR+1] != Real(0.0) && osm(i,j+sy,k) == 0) { sten[NR+1] = Real(0.0); },
                             int const sz = (cell_id(i,j,k-1) < 0) ? 2 : -2;
                             if (sten[NR+2] != Real(0.0) && osm(i,j,k+sz) == 0) { sten[NR+2] = Real(0.0); })
            }
            AlgInt n = 0;
            for (unsigned int m = 0; m < NS; ++m) { n += (sten[m] != Real(0.0)); }
            ncols(i,j,k) = n;
        }
    });

    BaseFab<GpuArray<AlgInt,NS>> colfab(bx, 1, (GpuArray<AlgInt,NS>*)cols);
    amrex::fill(colfab,
    [=] AMREX_GPU_HOST_DEVICE (GpuArray<AlgInt,NS>& sten, int i, int j, int k)
    {
        habec_cols(sten, i, j, k, bx, cell_id);
    });
}

}

#endif
