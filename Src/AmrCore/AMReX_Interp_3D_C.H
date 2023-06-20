#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pcinterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp);
                }
            }
        }
    }
}

namespace {
    constexpr int ix   = 0;
    constexpr int iy   = 1;
    constexpr int iz   = 2;
    constexpr int ixy  = 3;
    constexpr int ixz  = 4;
    constexpr int iyz  = 5;
    constexpr int ixyz = 6;
}

// Let's keep these nodal functions even though amrex is no longer using
// them, because they are used by other codes.

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/Real(ratio[0]);
    const Real ry = Real(1.)/Real(ratio[1]);
    const Real rz = Real(1.)/Real(ratio[2]);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    T dx00 = u(i+1,j,k,nu) - u(i,j,k,nu);
                    T d0x0 = u(i,j+1,k,nu) - u(i,j,k,nu);
                    T d00x = u(i,j,k+1,nu) - u(i,j,k,nu);

                    T dx10 = u(i+1,j+1,k,nu) - u(i,j+1,k,nu);
                    T dx01 = u(i+1,j,k+1,nu) - u(i,j,k+1,nu);
                    T d0x1 = u(i,j+1,k+1,nu) - u(i,j,k+1,nu);

                    T dx11 = u(i+1,j+1,k+1,nu) - u(i,j+1,k+1,nu);

                    slope(i,j,k,n+ncomp*ix  ) = rx*dx00;
                    slope(i,j,k,n+ncomp*iy  ) = ry*d0x0;
                    slope(i,j,k,n+ncomp*iz  ) = rz*d00x;
                    slope(i,j,k,n+ncomp*ixy ) = rx*ry*(dx10 - dx00);
                    slope(i,j,k,n+ncomp*ixz ) = rx*rz*(dx01 - dx00);
                    slope(i,j,k,n+ncomp*iyz ) = ry*rz*(d0x1 - d0x0);
                    slope(i,j,k,n+ncomp*ixyz) = rx*ry*rz*(dx11 - dx01 - dx10 + dx00);
                }
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_interp (Box const& bx, Array4<T> const& fine, const int fcomp, const int ncomp,
                  Array4<T const> const& slope, Array4<T const> const& crse,
                  const int ccomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto chi = amrex::ubound(slope);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::min(amrex::coarsen(k,ratio[2]),chi.z);
            const Real fz = k - kc*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
                const Real fy = j - jc*ratio[1];
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                    const Real fx = i - ic*ratio[0];
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp)
                        + fx*slope(ic,jc,kc,n+ncomp*ix)
                        + fy*slope(ic,jc,kc,n+ncomp*iy)
                        + fz*slope(ic,jc,kc,n+ncomp*iz)
                        + fx*fy*slope(ic,jc,kc,n+ncomp*ixy)
                        + fx*fz*slope(ic,jc,kc,n+ncomp*ixz)
                        + fy*fz*slope(ic,jc,kc,n+ncomp*iyz)
                        + fx*fy*fz*slope(ic,jc,kc,n+ncomp*ixyz);
                }
            }
        }
    }
}


template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_face_interp (int ci, int cj, int ck,
                     int nc, int nf, int idir,
                     Array4<T const> const& crse,
                     Array4<T> const& fine,
                     Array4<int const> const& mask,
                     IntVect const& ratio) noexcept
{

    if (mask) {
        if (!mask(ci, cj, ck, nc))
            { return; }
    }

    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];
    const int fk = ck*ratio[2];

    switch (idir) {
        case 0:
        {
            const Real ll = crse(ci, cj-1, ck-1, nc);
            const Real cl = crse(ci, cj-1, ck,   nc);
            const Real rl = crse(ci, cj-1, ck+1, nc);

            const Real lc = crse(ci, cj  , ck-1, nc);
            const Real cc = crse(ci, cj  , ck,   nc);
            const Real rc = crse(ci, cj  , ck+1, nc);

            const Real lu = crse(ci, cj+1, ck-1, nc);
            const Real cu = crse(ci, cj+1, ck,   nc);
            const Real ru = crse(ci, cj+1, ck+1, nc);

            fine(fi,fj  ,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi,fj  ,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi,fj+1,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi,fj+1,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
        case 1:
        {
            const Real ll = crse(ci-1, cj, ck-1, nc);
            const Real cl = crse(ci  , cj, ck-1, nc);
            const Real rl = crse(ci+1, cj, ck-1, nc);

            const Real lc = crse(ci-1, cj, ck  , nc);
            const Real cc = crse(ci  , cj, ck  , nc);
            const Real rc = crse(ci+1, cj, ck  , nc);

            const Real lu = crse(ci-1, cj, ck+1, nc);
            const Real cu = crse(ci  , cj, ck+1, nc);
            const Real ru = crse(ci+1, cj, ck+1, nc);

            fine(fi  ,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi+1,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi  ,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi+1,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
        case 2:
        {
            const Real ll = crse(ci-1, cj-1, ck, nc);
            const Real cl = crse(ci  , cj-1, ck, nc);
            const Real rl = crse(ci+1, cj-1, ck, nc);

            const Real lc = crse(ci-1, cj  , ck, nc);
            const Real cc = crse(ci  , cj  , ck, nc);
            const Real rc = crse(ci+1, cj  , ck, nc);

            const Real lu = crse(ci-1, cj+1, ck, nc);
            const Real cu = crse(ci  , cj+1, ck, nc);
            const Real ru = crse(ci+1, cj+1, ck, nc);

            fine(fi  ,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi+1,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi  ,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi+1,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int ci, int cj, int ck, int nf,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& fine,
             IntVect const& ratio,
             GpuArray<Real, AMREX_SPACEDIM> const& cellSize) noexcept
{
    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];
    const int fk = ck*ratio[1];

    // References to fine exterior values needed for interior calculation.
    const Real u000 = fine[0](fi,   fj  , fk  , nf);
    const Real u200 = fine[0](fi+2, fj  , fk  , nf);
    const Real u010 = fine[0](fi,   fj+1, fk  , nf);
    const Real u210 = fine[0](fi+2, fj+1, fk  , nf);
    const Real u001 = fine[0](fi,   fj  , fk+1, nf);
    const Real u201 = fine[0](fi+2, fj  , fk+1, nf);
    const Real u011 = fine[0](fi,   fj+1, fk+1, nf);
    const Real u211 = fine[0](fi+2, fj+1, fk+1, nf);

    const Real v000 = fine[1](fi  , fj  , fk  , nf);
    const Real v020 = fine[1](fi  , fj+2, fk  , nf);
    const Real v100 = fine[1](fi+1, fj  , fk  , nf);
    const Real v120 = fine[1](fi+1, fj+2, fk  , nf);
    const Real v001 = fine[1](fi  , fj  , fk+1, nf);
    const Real v021 = fine[1](fi  , fj+2, fk+1, nf);
    const Real v101 = fine[1](fi+1, fj  , fk+1, nf);
    const Real v121 = fine[1](fi+1, fj+2, fk+1, nf);

    const Real w000 = fine[2](fi  , fj  , fk  , nf);
    const Real w002 = fine[2](fi  , fj  , fk+2, nf);
    const Real w100 = fine[2](fi+1, fj  , fk  , nf);
    const Real w102 = fine[2](fi+1, fj  , fk+2, nf);
    const Real w010 = fine[2](fi  , fj+1, fk  , nf);
    const Real w012 = fine[2](fi  , fj+1, fk+2, nf);
    const Real w110 = fine[2](fi+1, fj+1, fk  , nf);
    const Real w112 = fine[2](fi+1, fj+1, fk+2, nf);

    const Real dx = cellSize[0];
    const Real dy = cellSize[1];
    const Real dz = cellSize[2];

    const Real dx3 = dx*dx*dx;
    const Real dy3 = dy*dy*dy;
    const Real dz3 = dz*dz*dz;

    // aspbs = a squared plus b squared
    const Real xspys = dx*dx + dy*dy;
    const Real yspzs = dy*dy + dz*dz;
    const Real zspxs = dz*dz + dx*dx;

    fine[0](fi+1, fj  , fk  , nf) = Real(0.5)*(u000+u200)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  +                dx3/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w000+w102-w002-w100)
                                  +                dx3/(8*dz*xspys)*(w010+w112-w012-w110);

    fine[0](fi+1, fj+1, fk  , nf) = Real(0.5)*(u010+u210)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  +                dx3/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w010+w112-w012-w110)
                                  +                dx3/(8*dz*xspys)*(w000+w102-w100-w002);

    fine[0](fi+1, fj  , fk+1, nf) = Real(0.5)*(u001+u201)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  +                dx3/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w000+w102-w002-w100)
                                  +                dx3/(8*dz*xspys)*(w010+w112-w012-w110);

    fine[0](fi+1, fj+1, fk+1, nf) = Real(0.5)*(u011+u211)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  +                dx3/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w010+w112-w012-w110)
                                  +                dx3/(8*dz*xspys)*(w000+w102-w002-w100);

    fine[1](fi  , fj+1, fk  , nf) = Real(0.5)*(v000+v020)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  +                dy3/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w000+w012-w002-w010)
                                  +                dy3/(8*dz*xspys)*(w100+w112-w102-w110);

    fine[1](fi+1, fj+1, fk  , nf) = Real(0.5)*(v001+v021)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  +                dy3/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w100+w112-w102-w110)
                                  +                dy3/(8*dz*xspys)*(w000+w012-w002-w010);

    fine[1](fi  , fj+1, fk+1, nf) = Real(0.5)*(v100+v120)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  +                dy3/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w000+w012-w002-w010)
                                  +                dy3/(8*dz*xspys)*(w100+w112-w102-w110);

    fine[1](fi+1, fj+1, fk+1, nf) = Real(0.5)*(v101+v121)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  +                dy3/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w100+w112-w102-w110)
                                  +                dy3/(8*dz*xspys)*(w000+w012-w002-w010);

    fine[2](fi  , fj  , fk+1, nf) = Real(0.5)*(w000+w002)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  +                dz3/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v000+v021-v001-v020)
                                  +                dz3/(8*dy*zspxs)*(v100+v121-v101-v120);

    fine[2](fi  , fj+1, fk+1, nf) = Real(0.5)*(w010+w012)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  +                dz3/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);

    fine[2](fi+1, fj  , fk+1, nf) = Real(0.5)*(w100+w102)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v000+v021-v001-v020)
                                  +                dz3/(8*dy*zspxs)*(v100+v121-v101-v120);

    fine[2](fi+1, fj+1, fk+1, nf) = Real(0.5)*(w110+w112)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccprotect_3d (int ic, int jc, int kc, int nvar,
                   Box const& fine_bx,
                   IntVect const& ratio,
                   Array4<T>       const& fine,
                   Array4<T const> const& fine_state) noexcept
{
    // Calculate bounds for interpolation
    Dim3 fnbxlo = lbound(fine_bx);
    Dim3 fnbxhi = ubound(fine_bx);
    int ilo = amrex::max(ratio[0]*ic,              fnbxlo.x);
    int ihi = amrex::min(ratio[0]*ic+(ratio[0]-1), fnbxhi.x);
    int jlo = amrex::max(ratio[1]*jc,              fnbxlo.y);
    int jhi = amrex::min(ratio[1]*jc+(ratio[1]-1), fnbxhi.y);
    int klo = amrex::max(ratio[2]*kc,              fnbxlo.z);
    int khi = amrex::min(ratio[2]*kc+(ratio[2]-1), fnbxhi.z);

    /*
     * Check if interpolation needs to be redone for derived components (n > 0)
     */
    for (int n = 1; n < nvar-1; ++n) {

        bool redo_me = false;
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    if ((fine_state(i,j,k,n) + fine(i,j,k,n)) < 0.0) {
                        redo_me = true;
                    }
                }
            }
        }

        /*
         * If all the fine values are non-negative after the original
         * interpolated correction, then we do nothing here.
         *
         * If any of the fine values are negative after the original
         * interpolated correction, then we do our best.
         */
        if (redo_me) {

            // Calculate number of fine cells
            int numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1);

            /*
             * First, calculate the following quantities:
             *
             * crseTot = volume-weighted sum of all interpolated values
             *           of the correction, which is equivalent to
             *           the total volume-weighted coarse correction
             *
             * SumN = volume-weighted sum of all negative values of fine_state
             *
             * SumP = volume-weighted sum of all positive values of fine_state
             */
            Real crseTot = 0.0;
            Real SumN = 0.0;
            Real SumP = 0.0;
            for         (int k = klo; k <= khi; ++k) {
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        crseTot += fine(i,j,k,n);
                        if (fine_state(i,j,k,n) <= 0.0) {
                            SumN += fine_state(i,j,k,n);
                        } else {
                            SumP += fine_state(i,j,k,n);
                        }
                    }
                }
            }

            if ( (crseTot > 0) && (crseTot > std::abs(SumN)) ) {

                /*
                 * Special case 1:
                 *
                 * Coarse correction > 0, and fine_state has some cells
                 * with negative values which will be filled before
                 * adding to the other cells.
                 *
                 * Use the correction to bring negative cells to zero,
                 * then distribute the remaining positive proportionally.
                 */
                for         (int k = klo; k <= khi; ++k) {
                    for     (int j = jlo; j <= jhi; ++j) {
                        for (int i = ilo; i <= ihi; ++i) {

                            // Fill in negative values first.
                            if (fine_state(i,j,k,n) < 0.0) {
                                fine(i,j,k,n) = -fine_state(i,j,k,n);
                            }

                            // Then, add the remaining positive proportionally.
                            if (SumP > 0.0) {
                                if (fine_state(i,j,k,n) > 0.0) {
                                    Real alpha = (crseTot - std::abs(SumN)) / SumP;
                                    fine(i,j,k,n) = alpha * fine_state(i,j,k,n);
                                }
                            } else { /* (SumP < 0) */
                                Real posVal = (crseTot - std::abs(SumN)) / (Real)numFineCells;
                                fine(i,j,k,n) += posVal;
                            }

                        }
                    }
                }

            } else if ( (crseTot > 0) && (crseTot < std::abs(SumN)) ) {

                /*
                 * Special case 2:
                 *
                 * Coarse correction > 0, and correction can not make
                 * them all positive.
                 *
                 * Add correction only to the negative cells
                 * in proportion to their magnitude, and
                 * don't add any correction to the states already positive.
                 */
                for         (int k = klo; k <= khi; ++k) {
                    for     (int j = jlo; j <= jhi; ++j) {
                        for (int i = ilo; i <= ihi; ++i) {
                            Real alpha = crseTot / std::abs(SumN);
                            if (fine_state(i,j,k,n) < 0.0) {
                                // Add correction to negative cells proportionally.
                                fine(i,j,k,n) = alpha * std::abs(fine_state(i,j,k,n));
                            } else {
                                // Don't touch the positive states.
                                fine(i,j,k,n) = 0.0;
                            }

                        }
                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) > SumP) ) {

                /*
                 * Special case 3:
                 *
                 * Coarse correction < 0, and fine_state DOES NOT have
                 * enough positive states to absorb it.
                 *
                 * Here we distribute the remaining negative amount
                 * in such a way as to make them all as close to the
                 * same negative value as possible.
                 */
                for         (int k = klo; k <= khi; ++k) {
                    for     (int j = jlo; j <= jhi; ++j) {
                        for (int i = ilo; i <= ihi; ++i) {

                            // We want to make them all as close to the same negative value.
                            Real negVal = (SumP + SumN + crseTot) / (Real)numFineCells;
                            fine(i,j,k,n) = negVal - fine_state(i,j,k,n);
                        }
                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) < SumP) &&
                        ((SumP+SumN+crseTot) > 0.0) )  {

                /*
                 * Special case 4:
                 *
                 * Coarse correction < 0, and fine_state has enough
                 * positive states to absorb all the negative
                 * correction *and* to redistribute to make
                 * negative cells positive.
                 */
                for         (int k = klo; k <= khi; ++k) {
                    for     (int j = jlo; j <= jhi; ++j) {
                        for (int i = ilo; i <= ihi; ++i) {

                            if (fine_state(i,j,k,n) < 0.0) {
                                // Absorb the negative correction
                                fine(i,j,k,n) = -fine_state(i,j,k,n);
                            } else {
                                // Redistribute the rest to make negative cells positive
                                Real alpha = (crseTot + SumN) / SumP;
                                fine(i,j,k,n) = alpha * fine_state(i,j,k,n);
                            }

                        }
                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) < SumP) &&
                        ((SumP+SumN+crseTot) < 0.0) )  {
                /*
                 * Special case 5:
                 *
                 * Coarse correction < 0, and fine_state has enough
                 * positive states to absorb all the negative
                 * correction, but not enough to fix the states
                 * already negative.
                 *
                 * Here we take a constant percentage away from each
                 * positive cell and don't touch the negatives.
                 */
                for         (int k = klo; k <= khi; ++k) {
                    for     (int j = jlo; j <= jhi; ++j) {
                        for (int i = ilo; i <= ihi; ++i) {

                            if (fine_state(i,j,k,n) > 0.0) {
                                // Don't touch the negatives
                                fine(i,j,k,n) = -fine_state(i,j,k,n);
                            } else {
                                // Take a constant percentage away from each positive cell
                                Real alpha = (crseTot + SumP) / SumN;
                                fine(i,j,k,n) = alpha * fine_state(i,j,k,n);
                            }

                        }
                    }
                }

            } // special cases
        } // redo_me
    } // n

    // Set sync for density (n=0) to sum of spec sync (1:nvar)
    for         (int k = klo; k <= khi; ++k) {
        for     (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) {
                fine(i,j,k,0) = 0.0;
                for (int n = 1; n < nvar-1; ++n) {
                    fine(i,j,k,0) += fine(i,j,k,n);
                }
            }
        }
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccquartic_interp (int i, int j, int k, int n,
                       Array4<Real const> const& crse,
                       Array4<Real>       const& fine) noexcept
{
    // Note: there are asserts in CellConservativeQuartic::interp()
    //       to check whether ratio is all equal to 2.

    constexpr Array1D<Real, -2, 2> cL = { -0.01171875_rt,  0.0859375_rt, 0.5_rt, -0.0859375_rt, 0.01171875_rt };

    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    int kc = amrex::coarsen(k,2);
    int irx = i - 2*ic; // = abs(i % 2);
    int jry = j - 2*jc; // = abs(j % 2);
    int krz = k - 2*kc; // = abs(k % 2);

    Array2D<Real, -2, 2, -2, 2> ctmp2;
    for     (int jj = -2; jj <= 2; ++jj) {
        for (int ii = -2; ii <= 2; ++ii) {
            ctmp2(ii,jj) = 2.0_rt * ( cL(-2)*crse(ic+ii,jc+jj,kc-2,n)
                                    + cL(-1)*crse(ic+ii,jc+jj,kc-1,n)
                                    + cL( 0)*crse(ic+ii,jc+jj,kc  ,n)
                                    + cL( 1)*crse(ic+ii,jc+jj,kc+1,n)
                                    + cL( 2)*crse(ic+ii,jc+jj,kc+2,n) );
            if (krz) {
                ctmp2(ii,jj) = 2.0_rt * crse(ic+ii,jc+jj,kc,n) - ctmp2(ii,jj);
            }
        } // ii
    } // jj

    Array1D<Real, -2, 2> ctmp;
    for (int ii = -2; ii <= 2; ++ii) {
        ctmp(ii) = 2.0_rt * ( cL(-2)*ctmp2(ii,-2)
                            + cL(-1)*ctmp2(ii,-1)
                            + cL( 0)*ctmp2(ii, 0)
                            + cL( 1)*ctmp2(ii, 1)
                            + cL( 2)*ctmp2(ii, 2) );
        if (jry) {
            ctmp(ii) = 2.0_rt * ctmp2(ii, 0) - ctmp(ii);
        }
    } // ii

    Real ftmp = 2.0_rt * ( cL(-2)*ctmp(-2)
                         + cL(-1)*ctmp(-1)
                         + cL( 0)*ctmp( 0)
                         + cL( 1)*ctmp( 1)
                         + cL( 2)*ctmp( 2) );
    if (irx) {
        ftmp = 2.0_rt * ctmp(0) - ftmp;
    }

    fine(i,j,k,n) = ftmp;

}

}  // namespace amrex


#endif
