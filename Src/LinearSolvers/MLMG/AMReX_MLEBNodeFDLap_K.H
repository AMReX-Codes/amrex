#ifndef AMREX_MLEBNODEFDLAP_K_H_
#define AMREX_MLEBNODEFDLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLEBNodeFDLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBNodeFDLap_2D_K.H>
#else
#include <AMReX_MLEBNodeFDLap_3D_K.H>
#endif

namespace amrex {

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x_doit (int i, int j, int k, Array4<Real> const& px,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx, F const& phieb, Real dxi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i+1,j,k) >= 0) {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i+1,j,k) < 0) {
        px(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        // A centroid of 1 marks an uncut edge.  That happens when the EB is
        // exactly aligned with the grid, so the node itself is on the boundary
        // and the full edge length is used.
        Real const h = (ecx(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) - Real(2.0) * ecx(i,j,k);
        px(i,j,k) = dxi * (p(i+1,j,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = (ecx(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) + Real(2.0) * ecx(i,j,k);
        px(i,j,k) = dxi * (phieb(i+1,j,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y_doit (int i, int j, int k, Array4<Real> const& py,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecy, F const& phieb, Real dyi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j+1,k) >= 0) {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j+1,k) < 0) {
        py(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        // A centroid of 1 marks an uncut edge.  That happens when the EB is
        // exactly aligned with the grid, so the node itself is on the boundary
        // and the full edge length is used.
        Real const h = (ecy(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) - Real(2.0) * ecy(i,j,k);
        py(i,j,k) = dyi * (p(i,j+1,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = (ecy(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) + Real(2.0) * ecy(i,j,k);
        py(i,j,k) = dyi * (phieb(i,j+1,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

#if (AMREX_SPACEDIM > 2)
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z_doit (int i, int j, int k, Array4<Real> const& pz,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecz, F const& phieb, Real dzi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j,k+1) >= 0) {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j,k+1) < 0) {
        pz(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        // A centroid of 1 marks an uncut edge.  That happens when the EB is
        // exactly aligned with the grid, so the node itself is on the boundary
        // and the full edge length is used.
        Real const h = (ecz(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) - Real(2.0) * ecz(i,j,k);
        pz(i,j,k) = dzi * (p(i,j,k+1) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = (ecz(i,j,k) == Real(1.0)) ? Real(1.0)
                     : Real(1.0) + Real(2.0) * ecz(i,j,k);
        pz(i,j,k) = dzi * (phieb(i,j,k+1) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                         Real phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ecx,
                                [=] (int, int, int) -> Real { return phieb; },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                         Array4<Real const> const& phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ecx,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecy,
                         Real phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ecy,
                                [=] (int, int, int) -> Real { return phieb; },
                                dyi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecy,
                         Array4<Real const> const& phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ecy,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dyi);
    });
}

#if (AMREX_SPACEDIM > 2)

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecz,
                         Real phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ecz,
                                [=] (int, int, int) -> Real { return phieb; },
                                dzi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecz,
                         Array4<Real const> const& phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ecz,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dzi);
    });
}

#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    });
}

#if (AMREX_SPACEDIM > 1)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    });
}
#endif

#if (AMREX_SPACEDIM > 2)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    });
}
#endif


#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)

namespace ebnodefdlap_detail {

//! One assembled matrix row.  Entry 0 is the diagonal, the rest are neighbors.
struct Row
{
    static constexpr int nmax = 2*AMREX_SPACEDIM+1;
    int n = 0;
    Real val[nmax] = {};
    Dim3 node[nmax] = {};
};

/**
 * \brief sigma averaged onto the face shared by a node and one of its neighbors.
 *
 * The face is made up of the 2^(AMREX_SPACEDIM-1) cells that touch both nodes.
 * With EB the average is volume weighted, matching Fapply.
 *
 * \param idim face normal direction
 * \param side 0 for the low-side neighbor, 1 for the high-side neighbor
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_sigma (int i, int j, int k, int idim, int side, bool has_eb,
                 Array4<Real const> const& sig, Array4<Real const> const& vfrc) noexcept
{
    Real num = Real(0.0);
    Real den = Real(0.0);
    for (int c = 0; c < (1 << (AMREX_SPACEDIM-1)); ++c) {
        int off[3] = {0,0,0};
        int bit = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d == idim) {
                off[d] = side - 1;
            } else {
                off[d] = ((c >> bit) & 1) - 1;
                ++bit;
            }
        }
        int const ii = i + off[0];
        int const jj = j + off[1];
        int const kk = k + off[2];
        Real const w = has_eb ? vfrc(ii,jj,kk) : Real(1.0);
        num += sig(ii,jj,kk) * w;
        den += w;
    }
    return num / den;
}

}

/**
 * \brief Assemble one row of the bottom-level MLEBNodeFDLaplacian matrix.
 *
 * The row must reproduce MLEBNodeFDLaplacian::Fapply with homogeneous boundary
 * data, because MLMG's bottom solve works on the residual correction.  Nodes
 * that are not part of the linear system (Dirichlet, EB covered, or owned by
 * another box) have an invalid \p gid and their columns are dropped, which is
 * the same as using a zero value there.  A node outside a Neumann boundary is
 * folded onto its mirror image, matching what MLNodeLinOp::applyBC puts into
 * those ghost cells.
 *
 * \param gid         global node ids, invalid (max) outside the linear system
 * \param has_eb      whether this box has cut cells
 * \param has_sig     whether sigma is given as a MultiFab
 * \param bcoef       sigma_d / dx_d^2 per direction
 * \param ndlo, ndhi  bounds of the nodal domain
 * \param reflect_lo  whether the low-side domain BC mirrors ghost nodes
 * \param reflect_hi  whether the high-side domain BC mirrors ghost nodes
 */
template <typename AtomicInt>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ebnodefdlap_detail::Row
mlebndfdlap_ijmat_row (int i, int j, int k,
                       Array4<AtomicInt const> const& gid,
                       bool has_eb, bool has_sig,
                       GpuArray<Real,AMREX_SPACEDIM> const& bcoef,
                       Array4<Real const> const& sig,
                       Array4<Real const> const& vfrc,
                       Array4<Real const> const& levset,
                       GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& ec,
                       Dim3 const& ndlo, Dim3 const& ndhi,
                       GpuArray<bool,AMREX_SPACEDIM> const& reflect_lo,
                       GpuArray<bool,AMREX_SPACEDIM> const& reflect_hi) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AtomicInt>::max();

    Real diag = Real(0.0);
    Real scale = Real(1.0);
    Real hm_arr[AMREX_SPACEDIM] = {};
    Real hp_arr[AMREX_SPACEDIM] = {};
    bool cut_lo_arr[AMREX_SPACEDIM] = {};
    bool cut_hi_arr[AMREX_SPACEDIM] = {};
    Real coef[2*AMREX_SPACEDIM] = {};
    bool used[2*AMREX_SPACEDIM] = {};

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const di = (idim == 0) ? 1 : 0;
        int const dj = (idim == 1) ? 1 : 0;
        int const dk = (idim == 2) ? 1 : 0;

        Real hp = Real(1.0);
        Real hm = Real(1.0);
        if (has_eb) {
            Real const ecp = ec[idim](i,j,k);
            Real const ecm = ec[idim](i-di,j-dj,k-dk);
            hp = (ecp == Real(1.0)) ? Real(1.0) : (Real(1.0) + Real(2.0)*ecp);
            hm = (ecm == Real(1.0)) ? Real(1.0) : (Real(1.0) - Real(2.0)*ecm);
            cut_hi_arr[idim] = (levset(i+di,j+dj,k+dk) >= Real(0.0));
            cut_lo_arr[idim] = (levset(i-di,j-dj,k-dk) >= Real(0.0));
            scale = amrex::min(scale, hm, hp);
        }
        hp_arr[idim] = hp;
        hm_arr[idim] = hm;
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const di = (idim == 0) ? 1 : 0;
        int const dj = (idim == 1) ? 1 : 0;
        int const dk = (idim == 2) ? 1 : 0;
        Real const hp = hp_arr[idim];
        Real const hm = hm_arr[idim];
        Real const fac = bcoef[idim] * Real(2.0) / (hp+hm);

        int const nd = (idim == 0) ? i : ((idim == 1) ? j : k);
        int const dlo = (idim == 0) ? ndlo.x : ((idim == 1) ? ndlo.y : ndlo.z);
        int const dhi = (idim == 0) ? ndhi.x : ((idim == 1) ? ndhi.y : ndhi.z);

        for (int side = 0; side < 2; ++side) {
            Real const s = has_sig
                ? ebnodefdlap_detail::face_sigma(i,j,k,idim,side,has_eb,sig,vfrc)
                : Real(1.0);
            bool const cut = (side == 0) ? cut_lo_arr[idim] : cut_hi_arr[idim];
            Real const h = (side == 0) ? hm : hp;
            if (cut) {
                // The neighbor is covered, so it only contributes a Dirichlet
                // value, which is zero here.
                diag -= s * fac * mlebndfdlap_scaled_hinv(scale, h);
            } else {
                diag -= s * fac * scale;
                int ii = (side == 0) ? i-di : i+di;
                int jj = (side == 0) ? j-dj : j+dj;
                int kk = (side == 0) ? k-dk : k+dk;
                int slot = 2*idim + side;
                bool const outside = (side == 0) ? (nd-1 < dlo) : (nd+1 > dhi);
                bool const reflect = (side == 0) ? reflect_lo[idim] : reflect_hi[idim];
                if (outside && reflect) {
                    ii = (side == 0) ? i+di : i-di;
                    jj = (side == 0) ? j+dj : j-dj;
                    kk = (side == 0) ? k+dk : k-dk;
                    slot = 2*idim + (1-side);
                }
                if (gid(ii,jj,kk) < gidmax) {
                    coef[slot] += s * fac * scale;
                    used[slot] = true;
                }
            }
        }
    }

    ebnodefdlap_detail::Row row;
    row.val[0] = diag;
    row.node[0] = Dim3{.x=i, .y=j, .z=k};
    row.n = 1;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const di = (idim == 0) ? 1 : 0;
        int const dj = (idim == 1) ? 1 : 0;
        int const dk = (idim == 2) ? 1 : 0;
        for (int side = 0; side < 2; ++side) {
            int const slot = 2*idim + side;
            if (used[slot]) {
                row.val[row.n] = coef[slot];
                row.node[row.n] = (side == 0) ? Dim3{.x=i-di, .y=j-dj, .z=k-dk}
                                              : Dim3{.x=i+di, .y=j+dj, .z=k+dk};
                ++row.n;
            }
        }
    }
    return row;
}

#endif
}

#endif
