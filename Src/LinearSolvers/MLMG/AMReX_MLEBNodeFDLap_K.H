#ifndef AMREX_MLEBNODEFDLAP_K_H_
#define AMREX_MLEBNODEFDLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

namespace amrex {
//! Open length of an edge measured from its high node.  See mlebndfdlap_edge_len_from_cent.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_hm (Real p) noexcept
{
    return (p == Real(1.0)) ? Real(1.0) : Real(1.0) - p;
}
}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLEBNodeFDLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBNodeFDLap_2D_K.H>
#else
#include <AMReX_MLEBNodeFDLap_3D_K.H>
#endif

namespace amrex {

//! Constant sigma, already folded into the stencil coefficients.
struct EBNodeFDConstSigma
{
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real x (int, int, int) noexcept { return Real(1.0); }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real y (int, int, int) noexcept { return Real(1.0); }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real z (int, int, int) noexcept { return Real(1.0); }
};

//! Constant sigma for RZ: sr in the radial direction and 1 in the axial direction.
struct EBNodeFDRZConstSigma
{
    Real sr = Real(1.0);
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real x (int, int, int) const noexcept { return sr; }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real y (int, int, int) noexcept { return Real(1.0); }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real z (int, int, int) noexcept { return Real(1.0); }
};

//! Edge-centered sigma.  s[d](i,j,k) is the edge from node (i,j,k) in direction d.
struct EBNodeFDEdgeSigma
{
    GpuArray<Array4<Real const>,AMREX_SPACEDIM> s;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real x (int i, int j, int k) const noexcept { return s[0](i,j,k); }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real y (int i, int j, int k) const noexcept { // NOLINT(readability-convert-member-functions-to-static)
#if (AMREX_SPACEDIM > 1)
        return s[1](i,j,k);
#else
        amrex::ignore_unused(i,j,k); return Real(1.0);
#endif
    }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real z (int i, int j, int k) const noexcept { // NOLINT(readability-convert-member-functions-to-static)
#if (AMREX_SPACEDIM > 2)
        return s[2](i,j,k);
#else
        amrex::ignore_unused(i,j,k); return Real(1.0);
#endif
    }
};

/**
 * \brief Convert an EB edge centroid into the edge length used by the kernels.
 *
 * The stored value p is the position of the EB on the edge measured from the
 * low node in units of dx, with 1 meaning the edge is not cut.  The open
 * length is p from the low node and mlebndfdlap_hm(p) from the high node.
 * A centroid of 1 is an uncut edge and -1 a covered edge; otherwise it is
 * the centroid of the open part, which lies next to the uncovered node.
 *
 * \param off unit vector in the edge direction
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_edge_len_from_cent (int i, int j, int k, Array4<Real> const& el,
                                     Array4<Real const> const& ec,
                                     Array4<Real const> const& levset,
                                     Dim3 const& off) noexcept
{
    Real const c = ec(i,j,k);
    if (c == Real(1.0)) {
        el(i,j,k) = Real(1.0);
    } else if (c == Real(-1.0)) {
        el(i,j,k) = Real(0.0);
    } else {
        bool const lo_covered = levset(i,j,k) >= Real(0.0);
        bool const hi_covered = levset(i+off.x,j+off.y,k+off.z) >= Real(0.0);
        Real p;
        if (lo_covered && hi_covered) {
            p = Real(0.0);
        } else if (lo_covered) {
            p = Real(2.0)*c;
        } else {
            p = Real(1.0) + Real(2.0)*c;
        }
        el(i,j,k) = amrex::min(Real(1.0), amrex::max(Real(0.0), p));
    }
}

/**
 * \brief Edge length of a coarse edge from the two fine edges it contains.
 *
 * When both coarse nodes are uncovered but the fine node between them is
 * covered, the coarse edge is blocked and the EB is placed at the middle of
 * the covered part so that the edge stays blocked from both sides.
 *
 * \param off   unit vector in the edge direction
 * \param ratio coarsening ratio; the fine edge is simply copied when the
 *              ratio in the edge direction is 1
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_coarsen_edge_len (int ic, int jc, int kc, Array4<Real> const& cel,
                                   Array4<Real const> const& fel,
                                   Array4<Real const> const& flevset,
                                   Dim3 const& off, Dim3 const& ratio) noexcept
{
    int const i = ic*ratio.x;
    int const j = jc*ratio.y;
    int const k = kc*ratio.z;
    if (off.x*ratio.x + off.y*ratio.y + off.z*ratio.z == 1) {
        cel(ic,jc,kc) = fel(i,j,k);
    } else {
        int const im = i+off.x;
        int const jm = j+off.y;
        int const km = k+off.z;
        int const ib = im+off.x;
        int const jb = jm+off.y;
        int const kb = km+off.z;
        bool const a_covered = flevset(i ,j ,k ) >= Real(0.0);
        bool const m_covered = flevset(im,jm,km) >= Real(0.0);
        bool const b_covered = flevset(ib,jb,kb) >= Real(0.0);
        Real const p1 = fel(i ,j ,k );
        Real const p2 = fel(im,jm,km);
        Real const hp1 = p1;                  // from a toward m
        Real const hm1 = mlebndfdlap_hm(p1);  // from m toward a
        Real const hp2 = p2;                  // from m toward b
        Real const hm2 = mlebndfdlap_hm(p2);  // from b toward m
        // open lengths from a and from b, in fine units
        Real const la = (!m_covered && hp1 == Real(1.0)) ? Real(1.0)+hp2 : hp1;
        Real const lb = (!m_covered && hm2 == Real(1.0)) ? Real(1.0)+hm1 : hm2;
        Real p;
        if (a_covered && b_covered) {
            p = Real(0.0);
        } else if (b_covered) {
            p = Real(0.5)*la;
        } else if (a_covered) {
            p = Real(1.0) - Real(0.5)*lb;
        } else if (la == Real(2.0)) {
            p = Real(1.0);
        } else {
            p = Real(0.5) + Real(0.25)*(la - lb);
        }
        cel(ic,jc,kc) = p;
    }
}

/**
 * \brief Edge sigma from cell-centered sigma.
 *
 * The 2^(D-1) cells touching the edge are averaged with weight
 * vfrc (1 if vfrc is empty) times the cell radius for RZ axial edges.
 *
 * \param idim edge direction
 * \param rz   whether the geometry is RZ; then direction 0 is r
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_cc_to_edge_sigma (int i, int j, int k, Array4<Real> const& se,
                                   Array4<Real const> const& sc,
                                   Array4<Real const> const& vfrc, int idim,
                                   bool rz, Real dr, Real rlo) noexcept
{
    Real num = Real(0.0);
    Real den = Real(0.0);
    for (int c = 0; c < (1 << (AMREX_SPACEDIM-1)); ++c) {
        int off[3] = {0,0,0};
        int bit = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != idim) {
                off[d] = ((c >> bit) & 1) - 1;
                ++bit;
            }
        }
        int const ii = i + off[0];
        int const jj = j + off[1];
        int const kk = k + off[2];
        Real w = vfrc ? vfrc(ii,jj,kk) : Real(1.0);
        if (rz && idim == 1) {
            w *= amrex::Math::abs(rlo + (Real(ii)+Real(0.5))*dr);
        }
        num += sc(ii,jj,kk) * w;
        den += w;
    }
    se(i,j,k) = (den > Real(0.0)) ? num/den : Real(0.0);
}

//! Total open length of a fine edge, used as a weight when coarsening sigma.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_open_len (int i, int j, int k, Array4<Real const> const& el,
                           Array4<Real const> const& levset, Dim3 const& off) noexcept
{
    if (!el) { return Real(1.0); }
    bool const lo_covered = levset(i,j,k) >= Real(0.0);
    bool const hi_covered = levset(i+off.x,j+off.y,k+off.z) >= Real(0.0);
    Real const p = el(i,j,k);
    if (lo_covered && hi_covered) {
        return Real(0.0);
    } else if (hi_covered) {
        return p;
    } else if (lo_covered) {
        return mlebndfdlap_hm(p);
    } else {
        return Real(1.0);
    }
}

/**
 * \brief Coarsen edge sigma: the open-length weighted harmonic mean of the two
 * fine edges making up the coarse edge.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_coarsen_edge_sigma (int ic, int jc, int kc, Array4<Real> const& cs,
                                     Array4<Real const> const& fs,
                                     Array4<Real const> const& fel,
                                     Array4<Real const> const& flevset,
                                     Dim3 const& off, Dim3 const& ratio) noexcept
{
    int const i = ic*ratio.x;
    int const j = jc*ratio.y;
    int const k = kc*ratio.z;
    if (off.x*ratio.x + off.y*ratio.y + off.z*ratio.z == 1) {
        cs(ic,jc,kc) = fs(i,j,k);
    } else {
        int const im = i+off.x;
        int const jm = j+off.y;
        int const km = k+off.z;
        Real const s1 = fs(i,j,k);
        Real const s2 = fs(im,jm,km);
        Real const h1 = mlebndfdlap_open_len(i ,j ,k , fel, flevset, off);
        Real const h2 = mlebndfdlap_open_len(im,jm,km, fel, flevset, off);
        if (h1 == Real(0.0)) {
            cs(ic,jc,kc) = s2;
        } else if (h2 == Real(0.0)) {
            cs(ic,jc,kc) = s1;
        } else {
            Real const den = h1*s2 + h2*s1;
            cs(ic,jc,kc) = (den > Real(0.0)) ? (h1+h2)*s1*s2/den : Real(0.0);
        }
    }
}

/**
 * \brief Fill ghost cells outside the domain by reflection.
 *
 * Works for any index type: a direction in which the data are nodal is
 * mirrored about the boundary node, a cell-centered direction about the
 * boundary face.  Periodic directions are skipped.
 *
 * \param flip_dir when the ghost cell is outside the domain in this direction,
 *                 an edge length p is mirrored to mlebndfdlap_hm(p); -1 disables this.
 */
inline
void mlebndfdlap_fill_domain_ghost (Box const& gbx, Array4<Real> const& a,
                                    Box const& domain,
                                    GpuArray<bool,AMREX_SPACEDIM> const& is_periodic,
                                    int flip_dir) noexcept
{
    if (domain.contains(gbx)) { return; }
    auto const dlo = amrex::lbound(domain);
    auto const dhi = amrex::ubound(domain);
    GpuArray<int,AMREX_SPACEDIM> nodal;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        nodal[d] = domain.type(d) == IndexType::NODE ? 1 : 0;
    }
    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int idx[3];
        idx[0] = i; idx[1] = j; idx[2] = k;
        int lo[3];
        lo[0] = dlo.x; lo[1] = dlo.y; lo[2] = dlo.z;
        int hi[3];
        hi[0] = dhi.x; hi[1] = dhi.y; hi[2] = dhi.z;
        bool outside = false;
        bool flip = false;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!is_periodic[d]) {
                if (idx[d] < lo[d]) {
                    idx[d] = 2*lo[d] - idx[d] - (1-nodal[d]);
                    outside = true;
                    flip = flip || (d == flip_dir);
                } else if (idx[d] > hi[d]) {
                    idx[d] = 2*hi[d] - idx[d] + (1-nodal[d]);
                    outside = true;
                    flip = flip || (d == flip_dir);
                }
            }
        }
        if (outside) {
            Real const v = a(idx[0],idx[1],idx[2]);
            a(i,j,k) = flip ? mlebndfdlap_hm(v) : v;
        }
    });
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x_doit (int i, int j, int k, Array4<Real> const& px,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& elx, F const& phieb, Real dxi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i+1,j,k) >= 0) {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i+1,j,k) < 0) {
        px(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(elx(i,j,k));
        px(i,j,k) = dxi * (p(i+1,j,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = elx(i,j,k);
        px(i,j,k) = dxi * (phieb(i+1,j,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y_doit (int i, int j, int k, Array4<Real> const& py,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ely, F const& phieb, Real dyi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j+1,k) >= 0) {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j+1,k) < 0) {
        py(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(ely(i,j,k));
        py(i,j,k) = dyi * (p(i,j+1,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = ely(i,j,k);
        py(i,j,k) = dyi * (phieb(i,j+1,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

#if (AMREX_SPACEDIM > 2)
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z_doit (int i, int j, int k, Array4<Real> const& pz,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& elz, F const& phieb, Real dzi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j,k+1) >= 0) {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j,k+1) < 0) {
        pz(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(elz(i,j,k));
        pz(i,j,k) = dzi * (p(i,j,k+1) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = elz(i,j,k);
        pz(i,j,k) = dzi * (phieb(i,j,k+1) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& elx,
                         Real phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, elx,
                                [=] (int, int, int) -> Real { return phieb; },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& elx,
                         Array4<Real const> const& phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, elx,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ely,
                         Real phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ely,
                                [=] (int, int, int) -> Real { return phieb; },
                                dyi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ely,
                         Array4<Real const> const& phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ely,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dyi);
    });
}

#if (AMREX_SPACEDIM > 2)

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& elz,
                         Real phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, elz,
                                [=] (int, int, int) -> Real { return phieb; },
                                dzi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& elz,
                         Array4<Real const> const& phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, elz,
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
 * \brief sigma on the edge between a node and one of its neighbors.
 *
 * \param idim edge direction
 * \param side 0 for the low-side neighbor, 1 for the high-side neighbor
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real edge_sigma (int i, int j, int k, int idim, int side, bool has_sig,
                 GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& sig) noexcept
{
    if (!has_sig) { return Real(1.0); }
    int const di = (idim == 0 && side == 0) ? 1 : 0;
    int const dj = (idim == 1 && side == 0) ? 1 : 0;
    int const dk = (idim == 2 && side == 0) ? 1 : 0;
    return sig[idim](i-di,j-dj,k-dk);
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
 * \param sig         edge sigma per direction
 * \param el          open edge lengths per direction
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
                       GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& sig,
                       Array4<Real const> const& levset,
                       GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& el,
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
            hp = el[idim](i,j,k);
            hm = mlebndfdlap_hm(el[idim](i-di,j-dj,k-dk));
            cut_hi_arr[idim] = (levset(i+di,j+dj,k+dk) >= Real(0.0)) || (hp != Real(1.0));
            cut_lo_arr[idim] = (levset(i-di,j-dj,k-dk) >= Real(0.0)) || (hm != Real(1.0));
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
            Real const s = ebnodefdlap_detail::edge_sigma(i,j,k,idim,side,has_sig,sig);
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
