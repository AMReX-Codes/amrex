#ifndef AMREX_MLEBNODEFDLAP_K_H_
#define AMREX_MLEBNODEFDLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#include <limits>

namespace amrex {
//! Open length of an edge measured from its high node.  See mlebndfdlap_eb_pos_from_cent.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_hm (Real p) noexcept
{
    return (p == Real(1.0)) ? Real(1.0) : Real(1.0) - p;
}

//! Largest EB position below 1.  It puts the EB on the high node, since 1 means uncut.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_pmax () noexcept
{
    return Real(1.0) - std::numeric_limits<Real>::epsilon();
}

/**
 * \brief Weight of a matrix row in RZ: r/dr, or 1/8 on the axis.
 *
 * Multiplying the RZ operator rows by it makes the radial couplings
 * symmetric.  It is dimensionless, so absolute tolerances keep the units
 * of the residual.  The matrix and fillRHS both apply it.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_rz_row_weight (int i, Real dr, Real rlo) noexcept
{
    Real const r = rlo + Real(i)*dr;
    return (r == Real(0.0)) ? Real(0.125) : r/dr;
}

//! Restriction for an MG level coarsened in direction \p cdir only.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_line_restriction (int i, int j, int k, Array4<Real> const& crse,
                                   Array4<Real const> const& fine,
                                   Array4<int const> const& msk, int cdir) noexcept
{
    Dim3 const off{.x = int(cdir == 0), .y = int(cdir == 1), .z = int(cdir == 2)};
    int const ii = i*(1+off.x);
    int const jj = j*(1+off.y);
    int const kk = k*(1+off.z);
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = Real(0.0);
    } else {
        crse(i,j,k) = Real(0.25)*(fine(ii-off.x,jj-off.y,kk-off.z) +
                                  fine(ii+off.x,jj+off.y,kk+off.z))
            +         Real(0.5) * fine(ii,jj,kk);
    }
}

/**
 * \brief Restriction of a residual whose rows are scaled near the EB.
 *
 * Full weighting in the coarsened directions (ratio 2; 1 means not
 * coarsened).  Each fine residual is multiplied by min(1, sc/sf), which
 * brings it to the scale of the coarse row without enlarging it.
 *
 * \param sf row scales of the fine nodes
 * \param sc row scale of the coarse node
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_eb_restriction (int ic, int jc, int kc, Array4<Real> const& crse,
                                 Array4<Real const> const& fine,
                                 Array4<int const> const& msk,
                                 Array4<Real const> const& sf, Real sc,
                                 Dim3 const& ratio) noexcept
{
    int const i = ic*ratio.x;
    int const j = jc*ratio.y;
    int const k = kc*ratio.z;
    if (msk(i,j,k)) {
        crse(ic,jc,kc) = Real(0.0);
        return;
    }
    int const nx = ratio.x-1;
    int const ny = ratio.y-1;
    int const nz = ratio.z-1;
    Real r = Real(0.0);
    for (int dk = -nz; dk <= nz; ++dk) {
        Real const wz = nz ? ((dk == 0) ? Real(0.5) : Real(0.25)) : Real(1.0);
        for (int dj = -ny; dj <= ny; ++dj) {
            Real const wy = ny ? ((dj == 0) ? Real(0.5) : Real(0.25)) : Real(1.0);
            for (int di = -nx; di <= nx; ++di) {
                Real const wx = nx ? ((di == 0) ? Real(0.5) : Real(0.25)) : Real(1.0);
                Real const s = sf(i+di,j+dj,k+dk);
                Real const f = (s > sc) ? sc/s : Real(1.0);
                r += wx*wy*wz*f*fine(i+di,j+dj,k+dk);
            }
        }
    }
    crse(ic,jc,kc) = r;
}

//! Linear interpolation for an MG level coarsened in direction \p cdir only.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_line_interpadd (int i, int j, int k, Array4<Real> const& fine,
                                 Array4<Real const> const& crse,
                                 Array4<int const> const& msk, int cdir) noexcept
{
    if (!msk(i,j,k)) {
        Dim3 const off{.x = int(cdir == 0), .y = int(cdir == 1), .z = int(cdir == 2)};
        int const ic = off.x ? amrex::coarsen(i,2) : i;
        int const jc = off.y ? amrex::coarsen(j,2) : j;
        int const kc = off.z ? amrex::coarsen(k,2) : k;
        bool const odd = (ic*(1+off.x) != i) || (jc*(1+off.y) != j) || (kc*(1+off.z) != k);
        if (odd) {
            fine(i,j,k) += Real(0.5)*(crse(ic,jc,kc) + crse(ic+off.x,jc+off.y,kc+off.z));
        } else {
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
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
 * \brief Convert an EB edge centroid into the EB position used by the kernels.
 *
 * The stored value p is the position of the EB on the edge measured from the
 * low node in units of dx, with 1 meaning the edge is not cut and
 * mlebndfdlap_pmax() that the EB is on the high node.  The open length is p
 * from the low node and mlebndfdlap_hm(p) from the high node.  A centroid of
 * 1 is an uncut edge and -1 a covered edge; otherwise it is the centroid of
 * the open part, which lies next to the uncovered node.
 *
 * \param off unit vector in the edge direction
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_eb_pos_from_cent (int i, int j, int k, Array4<Real> const& ebp,
                                     Array4<Real const> const& ec,
                                     Array4<Real const> const& levset,
                                     Dim3 const& off) noexcept
{
    Real const c = ec(i,j,k);
    if (c == Real(1.0)) {
        ebp(i,j,k) = Real(1.0);
    } else {
        bool const lo_covered = levset(i,j,k) >= Real(0.0);
        bool const hi_covered = levset(i+off.x,j+off.y,k+off.z) >= Real(0.0);
        Real p;
        if (lo_covered && hi_covered) {
            p = Real(0.0);
        } else if (lo_covered) {
            // A covered edge (c == -1) has the EB on the uncovered high node.
            p = (c == Real(-1.0)) ? Real(1.0) : Real(2.0)*c;
        } else {
            p = (c == Real(-1.0)) ? Real(0.0) : Real(1.0) + Real(2.0)*c;
        }
        ebp(i,j,k) = amrex::min(mlebndfdlap_pmax(), amrex::max(Real(0.0), p));
    }
}

/**
 * \brief EB position on a coarse edge from the two fine edges it contains.
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
void mlebndfdlap_coarsen_eb_pos (int ic, int jc, int kc, Array4<Real> const& cebp,
                                   Array4<Real const> const& febp,
                                   Array4<Real const> const& flevset,
                                   Dim3 const& off, Dim3 const& ratio) noexcept
{
    int const i = ic*ratio.x;
    int const j = jc*ratio.y;
    int const k = kc*ratio.z;
    if (off.x*ratio.x + off.y*ratio.y + off.z*ratio.z == 1) {
        cebp(ic,jc,kc) = febp(i,j,k);
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
        Real const p1 = febp(i ,j ,k );
        Real const p2 = febp(im,jm,km);
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
        // Only an uncut edge may be 1.
        bool const uncut = !a_covered && !b_covered && la == Real(2.0);
        cebp(ic,jc,kc) = uncut ? p : amrex::min(p, mlebndfdlap_pmax());
    }
}

/**
 * \brief Whether an EB feature is lost inside a coarse cell.
 *
 * True when a fine edge with an EB lies on no edge line of the coarse cell
 * while the cell's nodes are all open and its edges all uncut.  Multigrid
 * diverges near such a feature.
 *
 * \param ratio coarsening ratio, 1 in directions that are not coarsened
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool mlebndfdlap_eb_lost (int ic, int jc, int kc, Array4<Real const> const& clevset,
                          GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& cebp,
                          GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& febp,
                          Dim3 const& ratio) noexcept
{
#if (AMREX_SPACEDIM == 3)
    constexpr int cmax = 1;
#else
    constexpr int cmax = 0;
#endif
    for (int c = 0; c <= cmax; ++c) {
    for (int b = 0; b <= 1; ++b) {
    for (int a = 0; a <= 1; ++a) {
        if (clevset(ic+a,jc+b,kc+c) >= Real(0.0)) { return false; }
    }}}
    for (int c = 0; c <= cmax; ++c) {
    for (int b = 0; b <= 1; ++b) {
        if (cebp[0](ic,jc+b,kc+c) != Real(1.0)) { return false; }
        if (cebp[1](ic+b,jc,kc+c) != Real(1.0)) { return false; }
    }}
#if (AMREX_SPACEDIM == 3)
    for (int b = 0; b <= 1; ++b) {
    for (int a = 0; a <= 1; ++a) {
        if (cebp[2](ic+a,jc+b,kc) != Real(1.0)) { return false; }
    }}
#endif
    int const i0 = ic*ratio.x;
    int const j0 = jc*ratio.y;
    int const k0 = kc*ratio.z;
    // A fine edge is on a coarse edge line if its transverse positions are
    // at the cell boundary.
    for (int c = 0; c <= ratio.z*cmax; ++c) {
        bool const cb = (c == 0 || c == ratio.z);
        for (int b = 0; b <= ratio.y; ++b) {
            bool const bb = (b == 0 || b == ratio.y);
            for (int a = 0; a <= ratio.x; ++a) {
                bool const ab = (a == 0 || a == ratio.x);
                if (a < ratio.x && !(bb && cb) &&
                    febp[0](i0+a,j0+b,k0+c) != Real(1.0)) { return true; }
                if (b < ratio.y && !(ab && cb) &&
                    febp[1](i0+a,j0+b,k0+c) != Real(1.0)) { return true; }
#if (AMREX_SPACEDIM == 3)
                if (c < ratio.z && !(ab && bb) &&
                    febp[2](i0+a,j0+b,k0+c) != Real(1.0)) { return true; }
#endif
            }
        }
    }
    return false;
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
Real mlebndfdlap_open_len (int i, int j, int k, Array4<Real const> const& ebp,
                           Array4<Real const> const& levset, Dim3 const& off) noexcept
{
    if (!ebp) { return Real(1.0); }
    bool const lo_covered = levset(i,j,k) >= Real(0.0);
    bool const hi_covered = levset(i+off.x,j+off.y,k+off.z) >= Real(0.0);
    Real const p = ebp(i,j,k);
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
                                     Array4<Real const> const& febp,
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
        Real const h1 = mlebndfdlap_open_len(i ,j ,k , febp, flevset, off);
        Real const h2 = mlebndfdlap_open_len(im,jm,km, febp, flevset, off);
        if (h1 == Real(0.0)) {
            cs(ic,jc,kc) = s2;
        } else if (h2 == Real(0.0)) {
            cs(ic,jc,kc) = s1;
        } else {
            // Avoid s1*s2, which can overflow or underflow.
            bool const lo1 = s1 <= s2;
            Real const slo = lo1 ? s1 : s2;
            Real const shi = lo1 ? s2 : s1;
            Real const hlo = lo1 ? h1 : h2;
            Real const hhi = lo1 ? h2 : h1;
            Real const den = (shi > Real(0.0)) ? hlo + hhi*(slo/shi) : Real(0.0);
            cs(ic,jc,kc) = (den > Real(0.0)) ? slo*((h1+h2)/den) : Real(0.0);
        }
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x_doit (int i, int j, int k, Array4<Real> const& px,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ebpx, F const& phieb, Real dxi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i+1,j,k) >= 0) {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i+1,j,k) < 0) {
        px(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(ebpx(i,j,k));
        px(i,j,k) = dxi * (p(i+1,j,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = ebpx(i,j,k);
        px(i,j,k) = dxi * (phieb(i+1,j,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y_doit (int i, int j, int k, Array4<Real> const& py,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ebpy, F const& phieb, Real dyi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j+1,k) >= 0) {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j+1,k) < 0) {
        py(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(ebpy(i,j,k));
        py(i,j,k) = dyi * (p(i,j+1,k) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = ebpy(i,j,k);
        py(i,j,k) = dyi * (phieb(i,j+1,k) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}

#if (AMREX_SPACEDIM > 2)
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z_doit (int i, int j, int k, Array4<Real> const& pz,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ebpz, F const& phieb, Real dzi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j,k+1) >= 0) {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j,k+1) < 0) {
        pz(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        Real const h = mlebndfdlap_hm(ebpz(i,j,k));
        pz(i,j,k) = dzi * (p(i,j,k+1) - phieb(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    } else {
        Real const h = ebpz(i,j,k);
        pz(i,j,k) = dzi * (phieb(i,j,k+1) - p(i,j,k)) * mlebndfdlap_scaled_hinv(Real(1.0), h);
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpx,
                         Real phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ebpx,
                                [=] (int, int, int) -> Real { return phieb; },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpx,
                         Array4<Real const> const& phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ebpx,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpy,
                         Real phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ebpy,
                                [=] (int, int, int) -> Real { return phieb; },
                                dyi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpy,
                         Array4<Real const> const& phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ebpy,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dyi);
    });
}

#if (AMREX_SPACEDIM > 2)

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpz,
                         Real phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ebpz,
                                [=] (int, int, int) -> Real { return phieb; },
                                dzi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ebpz,
                         Array4<Real const> const& phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ebpz,
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


#if (AMREX_SPACEDIM > 1)

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
 * \param ebp          EB positions on the edges per direction
 * \param ndlo, ndhi  bounds of the nodal domain
 * \param reflect_lo  whether the low-side domain BC mirrors ghost nodes
 * \param reflect_hi  whether the high-side domain BC mirrors ghost nodes
 */
template <typename AlgGid>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ebnodefdlap_detail::Row
mlebndfdlap_ijmat_row (int i, int j, int k,
                       Array4<AlgGid const> const& gid,
                       bool has_eb, bool has_sig,
                       GpuArray<Real,AMREX_SPACEDIM> const& bcoef,
                       GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& sig,
                       Array4<Real const> const& levset,
                       GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& ebp,
                       Dim3 const& ndlo, Dim3 const& ndhi,
                       GpuArray<bool,AMREX_SPACEDIM> const& reflect_lo,
                       GpuArray<bool,AMREX_SPACEDIM> const& reflect_hi) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();

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
            hp = ebp[idim](i,j,k);
            hm = mlebndfdlap_hm(ebp[idim](i-di,j-dj,k-dk));
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

#if (AMREX_SPACEDIM == 2)
/**
 * \brief Assemble one row of the MLEBNodeFDLaplacian matrix in RZ.
 *
 * Same contract as mlebndfdlap_ijmat_row, reproducing the RZ kernels
 * mlebndfdlap_adotx_rz and mlebndfdlap_adotx_rz_eb.  The row is multiplied
 * by mlebndfdlap_rz_row_weight; fillRHS applies the same weight.
 *
 * \param sig0   constant radial sigma, used when \p has_sig is false
 * \param dr     radial cell size
 * \param dz     axial cell size divided by sqrt(sigma_z)
 * \param rlo    r of node 0
 * \param alpha  coefficient of the alpha/r^2 term
 */
template <typename AlgGid>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ebnodefdlap_detail::Row
mlebndfdlap_ijmat_row_rz (int i, int j, int k,
                          Array4<AlgGid const> const& gid,
                          bool has_eb, bool has_sig, Real sig0,
                          GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& sig,
                          Array4<Real const> const& levset,
                          GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& ebp,
                          Real dr, Real dz, Real rlo, Real alpha,
                          Dim3 const& ndlo, Dim3 const& ndhi,
                          GpuArray<bool,AMREX_SPACEDIM> const& reflect_lo,
                          GpuArray<bool,AMREX_SPACEDIM> const& reflect_hi) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();

    Real const r = rlo + Real(i)*dr;
    bool const axis = (r == Real(0.0));

    // Open length toward each neighbor (1 if uncut) and whether the EB
    // blocks that side.  Indexed [idim][side].
    Real h[2][2] = {{Real(1.0), Real(1.0)}, {Real(1.0), Real(1.0)}};
    bool cut[2][2] = {};
    Real scale = Real(1.0);
    if (has_eb) {
        h[0][0] = mlebndfdlap_hm(ebp[0](i-1,j,k));
        h[0][1] = ebp[0](i,j,k);
        h[1][0] = mlebndfdlap_hm(ebp[1](i,j-1,k));
        h[1][1] = ebp[1](i,j,k);
        cut[0][0] = (levset(i-1,j,k) >= Real(0.0)) || (h[0][0] != Real(1.0));
        cut[0][1] = (levset(i+1,j,k) >= Real(0.0)) || (h[0][1] != Real(1.0));
        cut[1][0] = (levset(i,j-1,k) >= Real(0.0)) || (h[1][0] != Real(1.0));
        cut[1][1] = (levset(i,j+1,k) >= Real(0.0)) || (h[1][1] != Real(1.0));
        scale = amrex::min(h[0][1], h[1][0], h[1][1]);
        if (!axis) { scale = amrex::min(scale, h[0][0]); }
    }

    Real diag = Real(0.0);
    Real coef[4] = {};
    bool used[4] = {};

    for (int idim = 0; idim < 2; ++idim) {
        for (int side = 0; side < 2; ++side) {
            // The axis stencil has no node inside the axis.
            if (idim == 0 && axis && side == 0) { continue; }

            Real const s = (idim == 0 && !has_sig)
                ? sig0 : ebnodefdlap_detail::edge_sigma(i,j,k,idim,side,has_sig,sig);
            Real const hh = h[idim][side];
            bool const blocked = cut[idim][side];
            // Coefficient of the neighbor, or of the EB value if blocked.
            Real c;
            if (idim == 1) {
                Real const fac = Real(2.0) / ((h[1][0]+h[1][1]) * dz*dz);
                c = s * fac * (blocked ? mlebndfdlap_scaled_hinv(scale, hh) : scale);
            } else if (axis) {
                Real const fac = Real(4.0) / (dr*dr);
                c = s * fac * (blocked ? mlebndfdlap_scaled_h2inv(scale, hh) : scale);
            } else {
                Real const fac = Real(2.0) / ((h[0][0]+h[0][1]) * r * dr*dr);
                Real const sgn = side ? Real(1.0) : Real(-1.0);
                c = blocked
                    ? s * fac * mlebndfdlap_scaled_hinv(scale, hh) * (r + sgn*Real(0.5)*hh*dr)
                    : s * fac * scale * (r + sgn*Real(0.5)*dr);
            }
            diag -= c;
            if (blocked) { continue; } // Dirichlet value, zero here

            int const nd = (idim == 0) ? i : j;
            int const dlo = (idim == 0) ? ndlo.x : ndlo.y;
            int const dhi = (idim == 0) ? ndhi.x : ndhi.y;
            bool const outside = (side == 0) ? (nd-1 < dlo) : (nd+1 > dhi);
            bool const reflect = (side == 0) ? reflect_lo[idim] : reflect_hi[idim];
            // A ghost node outside a Neumann boundary is its mirror image.
            int const ss = (outside && reflect) ? 1-side : side;
            int const ii = (idim == 0) ? i + (ss ? 1 : -1) : i;
            int const jj = (idim == 1) ? j + (ss ? 1 : -1) : j;
            if (gid(ii,jj,k) < gidmax) {
                coef[2*idim+ss] += c;
                used[2*idim+ss] = true;
            }
        }
    }
    if (!axis) { diag -= scale * alpha / (r*r); }

    Real const w = mlebndfdlap_rz_row_weight(i, dr, rlo);
    ebnodefdlap_detail::Row row;
    row.val[0] = diag * w;
    row.node[0] = Dim3{.x=i, .y=j, .z=k};
    row.n = 1;
    for (int idim = 0; idim < 2; ++idim) {
        for (int side = 0; side < 2; ++side) {
            int const slot = 2*idim + side;
            if (used[slot]) {
                int const d = side ? 1 : -1;
                row.val[row.n] = coef[slot] * w;
                row.node[row.n] = (idim == 0) ? Dim3{.x=i+d, .y=j, .z=k}
                                              : Dim3{.x=i, .y=j+d, .z=k};
                ++row.n;
            }
        }
    }
    return row;
}
#endif

#endif
}

#endif
