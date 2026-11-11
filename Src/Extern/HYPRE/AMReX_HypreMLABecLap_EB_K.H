#ifndef AMREX_HYPRE_ML_ABECLAP_EB_K_H_
#define AMREX_HYPRE_ML_ABECLAP_EB_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_BndryData.H>
#include <AMReX_Box.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LOUtil_K.H>
#include <AMReX_Loop.H>
#include <AMReX_MLLinOp_K.H>
#include <AMReX_Orientation.H>

#include <cmath>

namespace amrex {

/**
 * \brief The data of a box for the rows of HypreMLABecLap's composite EB
 * operator (hypmlabeclap_eb_row).  Null EB arrays mean regular cells, null
 * coefficients a = 0 and b = 1.
 */
struct HypMLABecLapEB
{
    Box vbx;
    bool irregular = false; // ghost cells are filled as in an EB fab
    IntVect rr{1};          // ratio to the coarser level
    GpuArray<Real,AMREX_SPACEDIM> dxinv{};
    Real sa = 0.0, sb = 0.0;
    IntVect fixed_pt = IntVect::TheMaxVector();
    Array4<Long const> cid;      // row ids, with the same-level neighbors in the ghost cells
    Array4<Long const> crse_id;  // row ids of the coarser level, on the coarsened box (null on level 0)
    Array4<int const> fine_mask; // covered by the finer level
    Array4<Real const> a, beb, phieb;
    GpuArray<Array4<Real const>,AMREX_SPACEDIM> b;
    Array4<EBCellFlag const> flag;
    Array4<Real const> vfrac, barea, bcent;
    Array4<Real const> scale;      // of the unknowns (null: 1)
    Array4<Real const> crse_scale; // of the coarser level's, on the coarsened box
    GpuArray<Array4<Real const>,AMREX_SPACEDIM> apert, fcent;
    GpuArray<int,2*AMREX_SPACEDIM> bct{};
    GpuArray<Real,2*AMREX_SPACEDIM> bcl{};
    GpuArray<bool,2*AMREX_SPACEDIM> cf_face{};   // boundary values from the coarser level's unknowns
    GpuArray<bool,2*AMREX_SPACEDIM> phys_face{}; // on the domain boundary (not periodic)
    GpuArray<Array4<int const>,2*AMREX_SPACEDIM> bmask;
    GpuArray<Array4<Real const>,2*AMREX_SPACEDIM> bval;
    //! Set to 1 when a row reaches a coarse cell off the coarser level's
    //! grids (the levels are not properly nested), to 2 on an internal error;
    //! the host checks it after the pass
    int* error = nullptr;
};

//! A ghost cell's value: c[0] times the boundary value plus the sum of
//! c[m]*phi(g+m*s) for m = 1..n-1, where s points from \p face into the box
//! (Neumann: c = {0, 1}).  n = 0: not filled, the value is zero.
struct HypMLABecLapEBFill
{
    int n = 0;
    Orientation face;
    Real c[3] = {};
};

/**
 * \brief How MLEBABecLap::applyBC (MLCellLinOp's for a regular fab, maxorder
 * 3) fills ghost cell g of the box.
 *
 * MLEBABecLap does not fill a cell outside two faces whose inner cell in
 * either direction is off the level's grids (a convex corner of the grids
 * at the coarse/fine boundary, an edge or a vertex in 3D; outside a physical
 * boundary and the coarse/fine boundary; a domain corner); MLMG then uses
 * the solution's ghost value, whatever it holds.  The composite operator
 * fills such a cell at the lowest order: outside two coarse/fine faces with
 * the coarse interpolation at the face, and at a physical boundary through
 * the physical face with the boundary value, or for Neumann the inner cell,
 * itself a ghost cell (c[1] applies to it); the physical face has priority
 * (a coarse/fine fill from outside the domain would reach coarse cells off
 * the grids).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
HypMLABecLapEBFill hypmlabeclap_eb_fill (HypMLABecLapEB const& e, IntVect const& g)
{
    HypMLABecLapEBFill f;
    // No break: as in MLEBABecLap::applyBC, a corner ghost cell in the boxes
    // of two faces gets the fill of the last direction
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        IntVect const ed = IntVect::TheDimensionVector(d);
        for (int side = 0; side < 2; ++side) {
            int const s = 1-2*side;
            Box gbx = (side == 0) ? amrex::adjCellLo(e.vbx, d) : amrex::adjCellHi(e.vbx, d);
            if (e.irregular) { gbx.grow(IntVect(1)-ed); }
            if (!gbx.contains(g)) { continue; }
            Orientation const face(d, (side == 0) ? Orientation::low : Orientation::high);
            if (e.cid(g+s*ed) < 0) {
                if (e.phys_face[face]) {
                    // The corner case above, through this physical face
                    f.face = face;
                    if (e.bct[face] == AMREX_LO_NEUMANN) {
                        f.c[0] = Real(0.0);
                        f.c[1] = Real(1.0);
                        f.n = 2;
                    } else {
                        f.c[0] = Real(1.0);
                        f.n = 1;
                    }
                } else if (e.cf_face[face] && f.n == 0) {
                    // Outside two coarse/fine faces: the coarse interpolation
                    // at this face (not overriding a physical face's fill)
                    f.face = face;
                    f.c[0] = Real(1.0);
                    f.n = 1;
                }
                continue;
            }
            f.face = face;
            if (e.bct[f.face] == AMREX_LO_NEUMANN) {
                f.c[0] = Real(0.0);
                f.c[1] = Real(1.0);
                f.n = 2;
                continue;
            }
            int n = amrex::min(e.vbx.length(d)+1, 3);
            if (e.irregular) {
                // A covered face along the normal lowers the order, a cut one
                // limits it to 2
                int order = 1;
                bool cut = false;
                for (int r = 0; r < n-1; ++r) {
                    IntVect fc = g;
                    fc[d] += (1-side) + s*r;
                    Real const ap = e.apert[d] ? e.apert[d](fc) : Real(1.0);
                    if (ap == Real(0.0)) { break; }
                    ++order;
                    cut = cut || (ap < Real(1.0));
                }
                n = cut ? amrex::min(order, 2) : order;
            }
            // Fallback MLMG does not have (it would read another box's
            // ghost cell): the inner cells of a corner ghost cell may be off
            // the grids
            for (int m = 1; m < n; ++m) {
                if (e.cid(g+(m*s)*ed) < 0) { n = m; }
            }
            Real x[3] = {-e.bcl[f.face]*e.dxinv[d], Real(0.5), Real(1.5)};
            poly_interp_coeff(Real(-0.5), x, n, f.c);
            f.n = n;
        }
    }
    return f;
}

/**
 * \brief The coarse cells (on the coarsened box) and weights of the
 * coarse/fine boundary value at ghost cell g of \p face, as InterpBndryData
 * interpolates it (order 3, half width 1).  Returns the number of cells.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int hypmlabeclap_eb_crse_interp (HypMLABecLapEB const& e, IntVect const& g, Orientation face,
                                 IntVect* civ, Real* w)
{
    auto const& msk = e.bmask[face];
    int const d = face.coordDir();
    IntVect const cg = amrex::coarsen(g, e.rr);
    civ[0] = cg;
    w[0] = Real(1.0);
    int n = 1;
    GpuArray<Real,AMREX_SPACEDIM> x{};
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (t == d) { continue; }
        IntVect const et = IntVect::TheDimensionVector(t);
        bool const lo = msk(g-e.rr[t]*et) == BndryData::not_covered;
        bool const hi = msk(g+e.rr[t]*et) == BndryData::not_covered;
        x[t] = Real(-0.5) + (Real(g[t]-cg[t]*e.rr[t])+Real(0.5))/Real(e.rr[t]);
        if (lo && hi) {
            w[0] -= x[t]*x[t];
            civ[n] = cg-et; w[n++] = Real(0.5)*x[t]*(x[t]-Real(1.0));
            civ[n] = cg+et; w[n++] = Real(0.5)*x[t]*(x[t]+Real(1.0));
        } else if (lo) {
            w[0] += x[t];
            civ[n] = cg-et; w[n++] = -x[t];
        } else if (hi) {
            w[0] -= x[t];
            civ[n] = cg+et; w[n++] = x[t];
        }
    }
#if (AMREX_SPACEDIM == 3)
    IntVect const e1 = IntVect::TheDimensionVector((d+1)%3);
    IntVect const e2 = IntVect::TheDimensionVector((d+2)%3);
    IntVect const r1 = e.rr[(d+1)%3]*e1;
    IntVect const r2 = e.rr[(d+2)%3]*e2;
    if (msk(g-r1-r2) == BndryData::not_covered && msk(g+r1-r2) == BndryData::not_covered &&
        msk(g-r1+r2) == BndryData::not_covered && msk(g+r1+r2) == BndryData::not_covered)
    {
        Real const f = Real(0.25)*x[(d+1)%3]*x[(d+2)%3];
        civ[n] = cg+e1+e2; w[n++] =  f;
        civ[n] = cg-e1+e2; w[n++] = -f;
        civ[n] = cg-e1-e2; w[n++] =  f;
        civ[n] = cg+e1-e2; w[n++] = -f;
    }
#endif
    AMREX_ASSERT(n <= 9);
    return n;
}

/**
 * \brief The entries of the Dirichlet flux on the EB of cut cell c (b_eb on
 * the EB, see MLEBABecLap's mlebabeclap_adotx): the gradient between the EB
 * value and phi interpolated at a point inside the cell along the EB normal,
 * passed to put(cell, value); the EB value's term goes to \p bcrhs.
 */
template <typename P>
AMREX_GPU_HOST_DEVICE
void hypmlabeclap_eb_dirichlet (HypMLABecLapEB const& e, IntVect const& c, Real kappa,
                                P const& put, Real& bcrhs)
{
    GpuArray<Real,AMREX_SPACEDIM> nrm;
    Real anorm = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        IntVect const ed = IntVect::TheDimensionVector(d);
        nrm[d] = e.apert[d](c) - e.apert[d](c+ed);
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            if (t != d) { nrm[d] *= e.dxinv[0]/e.dxinv[t]; }
        }
        anorm += nrm[d]*nrm[d];
    }
    anorm = std::sqrt(anorm);
    Real mmax = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        nrm[d] /= anorm;
        mmax = amrex::max(mmax, std::abs(nrm[d]*e.dxinv[d]));
    }
    Real const dg = get_dx_eb(kappa)/mmax;
#if (AMREX_SPACEDIM == 2)
    Real const area = get_eb_area(e.barea(c), nrm[0], nrm[1], e.dxinv[0]/e.dxinv[1]);
    Real const fb = e.sb*e.dxinv[1]*area*e.beb(c)/dg;
#else
    Real const area = get_eb_area(e.barea(c), nrm[0], nrm[1], nrm[2],
                                  e.dxinv[0]/e.dxinv[1], e.dxinv[0]/e.dxinv[2]);
    Real const fb = e.sb*e.dxinv[1]*(e.dxinv[2]/e.dxinv[0])*area*e.beb(c)/dg;
#endif
    // Multilinear interpolation at the point: the 2^D cells around it
    for (int m = 0; m < (1 << AMREX_SPACEDIM); ++m) {
        Real v = fb;
        IntVect q = c;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real const s = std::copysign(Real(1.0), nrm[d]);
            Real const gs = s*(e.bcent(c,d) - dg*nrm[d]*e.dxinv[d]);
            if (m & (1 << d)) {
                v *= -gs;
                q[d] -= int(s);
            } else {
                v *= Real(1.0)+gs;
            }
        }
        put(q, v);
    }
    if (e.phieb) { bcrhs += fb*e.phieb(c); }
}

/**
 * \brief The composite EB row of cell c and the cell's entries in coarser
 * rows, passed to put(row, col, value).
 *
 * The row is MLEBABecLap's operator (cell-centered phi, face-centered b)
 * times the volume fraction.  The fluxes across coarse/fine faces are the
 * finer level's: a fine cell adds them to its coarse neighbor's row, which
 * skips that face.  Ghost cells are filled as MLEBABecLap fills them, with
 * coarse/fine boundary values interpolated from the coarser level's unknowns
 * and other boundary values moved to the right-hand side (\p bcrhs).  A cell
 * covered by the finer level has the average-down row of EB_average_down;
 * the fine cells add their weights to it.  The same-level entries of the
 * cell's own row are summed first.  With the scaling s of the unknowns, an
 * entry (i,j) is multiplied by s_i*s_j.
 *
 * In order: the cell's entry in its parent's average-down row; the trivial
 * rows (covered cells); the face fluxes (also into coarse neighbors' rows);
 * the EB flux; the cell's own stencil.
 */
template <typename F>
AMREX_GPU_HOST_DEVICE
void hypmlabeclap_eb_row (HypMLABecLapEB const& e, IntVect const& c, F const& put, Real& bcrhs)
{
    constexpr int nsten = AMREX_D_TERM(3,*3,*3);
    Long const row = e.cid(c);
    bool const pinned = (c == e.fixed_pt);
    bcrhs = Real(0.0);

    auto scale = [&] (Array4<Real const> const& s, IntVect const& q)
    {
        return s ? s(q) : Real(1.0);
    };
    // The coarser level must cover the coarse cells the row reaches
    // (proper nesting); the host checks the flag after the pass
    auto crse_row = [&] (IntVect const& cq)
    {
        Long const id = e.crse_id(cq);
        if (id < 0 && e.error) { *e.error = 1; }
        return id;
    };

    Real const kappa = e.vfrac(c);
    Real const sc = scale(e.scale, c);
    if (e.crse_id && kappa > Real(0.0)) {
        // The grids are coarsenable, so the fine cells of pc are in the box
        IntVect const pc = amrex::coarsen(c, e.rr);
        Real vtot = Real(0.0);
        amrex::Loop(Box(pc*e.rr, pc*e.rr+e.rr-1), [&] (int i, int j, int k)
        {
            vtot += e.vfrac(i,j,k);
        });
        put(crse_row(pc), row, -kappa/vtot*scale(e.crse_scale, pc)*sc);
    }

    if (kappa == Real(0.0) || (e.fine_mask && e.fine_mask(c))) {
        put(row, row, sc*sc);
        return;
    }

    // A row with id and scaling s.  The cell's own row sums its same-level
    // entries in sten.  put_crse adds a coarse cell's entry except to the
    // off-diagonals of a pinned row, put_cell maps a cell of this level to
    // sten or to another row, put_phi also expands ghost cells.
    struct Row { Long id; Real s; };
    Row const own{row, sc};
    Real sten[nsten] = {};
    bool used[nsten] = {};
    auto put_crse = [&] (Row const& t, IntVect const& cq, Real v)
    {
        if (t.id != row || !pinned) {
            put(t.id, crse_row(cq), t.s*v*scale(e.crse_scale, cq));
        }
    };
    auto put_cell = [&] (Row const& t, IntVect const& q, Real v)
    {
        if (t.id == row) {
            IntVect const o = q-c+1;
            int const m = AMREX_D_TERM(o[0], +3*o[1], +9*o[2]);
            sten[m] += v;
            used[m] = true;
        } else {
            put(t.id, e.cid(q), t.s*v*scale(e.scale, q));
        }
    };
    auto put_phi = [&] (Row const& t, IntVect const& q0, Real v0)
    {
        // The fill of a ghost cell may use another ghost cell (a corner
        // ghost cell outside a Neumann wall at the coarse/fine boundary,
        // see hypmlabeclap_eb_fill), whose own fill uses cells on the grids
        // or the coarser level: a work list of at most two cells
        IntVect qs[2];
        Real vs[2];
        int ns = 0;
        qs[ns] = q0;
        vs[ns++] = v0;
        while (ns > 0) {
            IntVect const q = qs[--ns];
            Real const v = vs[ns];
            if (e.cid(q) >= 0) {
                put_cell(t, q, v);
                continue;
            }
            auto const f = hypmlabeclap_eb_fill(e, q);
            if (f.n == 0) { continue; }
            int const d = f.face.coordDir();
            int const s = f.face.isLow() ? 1 : -1;
            for (int m = 1; m < f.n; ++m) {
                IntVect qm = q;
                qm[d] += m*s;
                if (e.cid(qm) >= 0) {
                    put_cell(t, qm, v*f.c[m]);
                } else {
                    AMREX_ASSERT(ns < 2);
                    qs[ns] = qm;
                    vs[ns++] = v*f.c[m];
                }
            }
            if (e.bct[f.face] == AMREX_LO_NEUMANN) { continue; }
            if (e.cf_face[f.face]) {
                IntVect civ[9]; // 1 + 2 per transverse direction + 4 cross terms
                Real w[9];
                int const nc = hypmlabeclap_eb_crse_interp(e, q, f.face, civ, w);
                for (int m = 0; m < nc; ++m) {
                    put_crse(t, civ[m], v*f.c[0]*w[m]);
                }
            } else if (t.id == row) {
                bcrhs -= v*f.c[0]*e.bval[f.face](q);
            } else if (e.error) {
                // A coarse neighbor's row never reaches a physical boundary:
                // the centroid shift of a cut face toward a ghost cell off
                // the grids is excluded below (as MLEBABecLap's cc_mask)
                *e.error = 2;
            }
        }
    };

    if (e.a) { put_cell(own, c, kappa*e.sa*e.a(c)); }

    Real const rrvol = Real(AMREX_D_TERM(e.rr[0],*e.rr[1],*e.rr[2]));
    bool open = false; // a face with a nonzero aperture
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        IntVect const ed = IntVect::TheDimensionVector(d);
        for (int side = 0; side < 2; ++side) {
            Orientation const face(d, (side == 0) ? Orientation::low : Orientation::high);
            IntVect const f = (side == 0) ? c : c+ed;
            IntVect const nb = (side == 0) ? c-ed : c+ed;
            Real const ap = e.apert[d] ? e.apert[d](f) : Real(1.0);
            if (ap == Real(0.0)) { continue; }
            open = true;
            // The flux across a coarse/fine face is the fine cells'
            if (e.fine_mask && e.fine_mask(nb)) { continue; }
            // The coarse neighbor across a coarse/fine face
            IntVect const cnb = amrex::coarsen(nb, e.rr);
            Row const crow = (e.cid(nb) < 0 && e.cf_face[face])
                ? Row{crse_row(cnb), scale(e.crse_scale, cnb)} : Row{-1, Real(0.0)};
            // b*(phi(f)-phi(f-ed)) on a cut face is interpolated toward its
            // centroid from the neighboring faces; as MLEBABecLap's cc_mask,
            // not toward a face whose two cells are both off the level's grids
            IntVect off[AMREX_SPACEDIM-1];
            Real frac[AMREX_SPACEDIM-1];
            int nt = 0;
            if (ap < Real(1.0)) {
                for (int t = 0, it = 0; t < AMREX_SPACEDIM; ++t) {
                    if (t == d) { continue; }
                    Real const fc = e.fcent[d](f,it++);
                    IntVect const o = IntVect::TheDimensionVector(t)
                        * int(std::copysign(Real(1.0), fc));
                    if (fc != Real(0.0) && (e.cid(f+o-ed) >= 0 || e.cid(f+o) >= 0)) {
                        off[nt] = o;
                        frac[nt++] = std::abs(fc);
                    }
                }
            }
            Real const fac = ((side == 0) ? Real(1.0) : Real(-1.0))
                * e.sb*e.dxinv[d]*e.dxinv[d]*ap;
            for (int m = 0; m < (1 << nt); ++m) {
                Real v = fac;
                IntVect fp = f;
                for (int t = 0; t < nt; ++t) {
                    if (m & (1 << t)) {
                        v *= frac[t];
                        fp += off[t];
                    } else {
                        v *= Real(1.0)-frac[t];
                    }
                }
                if (e.b[d]) { v *= e.b[d](fp); }
                put_phi(own, fp, v);
                put_phi(own, fp-ed, -v);
                if (crow.id >= 0) {
                    put_phi(crow, fp, -v/rrvol);
                    put_phi(crow, fp-ed, v/rrvol);
                }
            }
        }
    }

    if (e.beb && e.flag(c).isSingleValued()) {
        hypmlabeclap_eb_dirichlet(e, c, kappa, [&] (IntVect const& q, Real v) { put_phi(own, q, v); },
                                  bcrhs);
    }

    // A fluid cell with all faces covered and no EB Dirichlet or a term (an
    // isolated cell) has no entries: give it the identity row.  Not a cell
    // whose open faces all border the finer level: its entries are the fine
    // cells' fluxes, added to its row by the fine cells
    if (!open && !used[nsten/2]) { put(row, row, sc*sc); }
    for (int m = 0; m < nsten; ++m) {
        if (used[m] && (!pinned || m == nsten/2)) {
            IntVect const o(AMREX_D_DECL(m%3-1, (m/3)%3-1, m/9-1));
            put(row, e.cid(c+o), sc*sten[m]*scale(e.scale, c+o));
        }
    }
}

/**
 * \brief The volume-fraction-weighted average of the fine cells of coarse
 * cell (i,j,k) of the coarsened fine box: the average-down row of a cell
 * covered by the finer level, solved exactly.  Zero if all the fine cells are
 * covered by the EB.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_eb_avgdown (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine,
                              Array4<Real const> const& vfrac, IntVect const& rr)
{
    IntVect const lo = IntVect(AMREX_D_DECL(i,j,k))*rr;
    Real sum = Real(0.0);
    Real vtot = Real(0.0);
    amrex::Loop(Box(lo, lo+rr-1), [&] (int ii, int jj, int kk)
    {
        sum += vfrac(ii,jj,kk)*fine(ii,jj,kk);
        vtot += vfrac(ii,jj,kk);
    });
    crse(i,j,k) = (vtot > Real(0.0)) ? sum/vtot : Real(0.0);
}

}

#endif
