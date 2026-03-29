#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>
#include <AMReX_LUSolver.H>

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

/// \cond DOXYGEN_IGNORE
namespace interp_detail {
    constexpr int ix   = 0;
    constexpr int iy   = 1;
    constexpr int iz   = 2;
    constexpr int ixy  = 3;
    constexpr int ixz  = 4;
    constexpr int iyz  = 5;
    constexpr int ixyz = 6;
}
/// \endcond

// Let's keep these nodal functions even though amrex is no longer using
// them, because they are used by other codes.

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    using namespace interp_detail;

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
    using namespace interp_detail;

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
//            const Real ll = crse(ci, cj-1, ck-1, nc);
//            const Real cl = crse(ci, cj-1, ck,   nc);
//            const Real rl = crse(ci, cj-1, ck+1, nc);
//
//            const Real lc = crse(ci, cj  , ck-1, nc);
//            const Real cc = crse(ci, cj  , ck,   nc);
//            const Real rc = crse(ci, cj  , ck+1, nc);

//            const Real lu = crse(ci, cj+1, ck-1, nc);
//            const Real cu = crse(ci, cj+1, ck,   nc);
//            const Real ru = crse(ci, cj+1, ck+1, nc);

//            fine(fi,fj  ,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
//            fine(fi,fj  ,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
//            fine(fi,fj+1,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
//            fine(fi,fj+1,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            const Real slpy  = Real(0.5)* (crse(ci,cj+1,ck  ,nc) - crse(ci,cj-1,ck  ,nc));
            const Real slpz  = Real(0.5)* (crse(ci,cj,  ck+1,nc) - crse(ci,cj,  ck-1,nc));
            const Real slpyz = Real(0.25)*(crse(ci,cj+1,ck+1,nc) - crse(ci,cj+1,ck-1,nc) - crse(ci,cj-1,ck+1,nc) + crse(ci,cj-1,ck-1,nc));

            const Real dys = Real(1.)/ratio[1];
            const Real dzs = Real(1.)/ratio[2];

            for (int k=0; k<ratio[2]; ++k){
            for (int j=0; j<ratio[1]; ++j){

              Real ysloc = (j+Real(0.5))*dys - Real(0.5);
              Real zsloc = (k+Real(0.5))*dzs - Real(0.5);
              fine(fi,fj+j,fk+k,nc) = crse(ci,cj,ck,nc) +  ysloc*slpy + zsloc*slpz + ysloc*zsloc*slpyz;

            }
            }

            break;
        }
        case 1:
        {
//            const Real ll = crse(ci-1, cj, ck-1, nc);
//            const Real cl = crse(ci  , cj, ck-1, nc);
//            const Real rl = crse(ci+1, cj, ck-1, nc);

//            const Real lc = crse(ci-1, cj, ck  , nc);
//            const Real cc = crse(ci  , cj, ck  , nc);
//            const Real rc = crse(ci+1, cj, ck  , nc);

//            const Real lu = crse(ci-1, cj, ck+1, nc);
//            const Real cu = crse(ci  , cj, ck+1, nc);
//            const Real ru = crse(ci+1, cj, ck+1, nc);

//            fine(fi  ,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
//            fine(fi+1,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
//            fine(fi  ,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
//            fine(fi+1,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            const Real slpx  = Real(0.5)* (crse(ci+1,cj,ck  ,nc) - crse(ci-1,cj,ck  ,nc));
            const Real slpz  = Real(0.5)* (crse(ci  ,cj,ck+1,nc) - crse(ci  ,cj,ck-1,nc));
            const Real slpxz = Real(0.25)*(crse(ci+1,cj,ck+1,nc) - crse(ci+1,cj,ck-1,nc) - crse(ci-1,cj,ck+1,nc) + crse(ci-1,cj,ck-1,nc));

            const Real dxs = Real(1.)/ratio[0];
            const Real dzs = Real(1.)/ratio[2];

            for (int k=0; k<ratio[2]; ++k){
            for (int i=0; i<ratio[0]; ++i){

              Real xsloc = (i+Real(0.5))*dxs - Real(0.5);
              Real zsloc = (k+Real(0.5))*dzs - Real(0.5);
              fine(fi+i,fj,fk+k,nc) = crse(ci,cj,ck,nc) + xsloc*slpx + zsloc*slpz + xsloc*zsloc*slpxz;

            }
            }


            break;
        }
        case 2:
        {
//            const Real ll = crse(ci-1, cj-1, ck, nc);
//            const Real cl = crse(ci  , cj-1, ck, nc);
//            const Real rl = crse(ci+1, cj-1, ck, nc);

//            const Real lc = crse(ci-1, cj  , ck, nc);
//            const Real cc = crse(ci  , cj  , ck, nc);
//            const Real rc = crse(ci+1, cj  , ck, nc);

//            const Real lu = crse(ci-1, cj+1, ck, nc);
//            const Real cu = crse(ci  , cj+1, ck, nc);
//            const Real ru = crse(ci+1, cj+1, ck, nc);

//            fine(fi  ,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
//            fine(fi+1,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
//            fine(fi  ,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
//            fine(fi+1,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            const Real slpx  = Real(0.5)* (crse(ci+1,cj,  ck,nc) - crse(ci-1,cj,  ck  ,nc));
            const Real slpy  = Real(0.5)* (crse(ci  ,cj+1,ck,nc) - crse(ci,  cj-1,ck  ,nc));
            const Real slpxy = Real(0.25)*(crse(ci+1,cj+1,ck,nc) - crse(ci+1,cj-1,ck,nc) - crse(ci-1,cj+1,ck,nc) + crse(ci-1,cj-1,ck,nc));

            const Real dxs = Real(1.)/ratio[0];
            const Real dys = Real(1.)/ratio[1];

            for (int j=0; j<ratio[2]; ++j){
            for (int i=0; i<ratio[0]; ++i){

              Real xsloc = (i+Real(0.5))*dxs - Real(0.5);
              Real ysloc = (j+Real(0.5))*dys - Real(0.5);
              fine(fi+i,fj+j,fk,nc) = crse(ci,cj,ck,nc) +  xsloc*slpx + ysloc*slpy + xsloc*ysloc*slpxy;

            }
            }


            break;
        }
        default : { break; }
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
    const int fk = ck*ratio[2];

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

    fine[1](fi+1, fj+1, fk  , nf) = Real(0.5)*(v100+v120)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  +                dy3/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w100+w112-w102-w110)
                                  +                dy3/(8*dz*xspys)*(w000+w012-w002-w010);

    fine[1](fi  , fj+1, fk+1, nf) = Real(0.5)*(v001+v021)
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
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v000+v021-v001-v020)
                                  +                dz3/(8*dy*zspxs)*(v100+v121-v101-v120);

    fine[2](fi+1, fj  , fk+1, nf) = Real(0.5)*(w100+w102)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  +                dz3/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);

    fine[2](fi+1, fj+1, fk+1, nf) = Real(0.5)*(w110+w112)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);
}


template <int i, std::enable_if_t<(i>=0 && i<8),int> = 0>
constexpr IntVect facediv_get_refratio ()
{
    if constexpr (i == 0) {
        return IntVect(2,2,2);
    } else if constexpr (i == 1) {
        return IntVect(4,2,2);
    } else if constexpr (i == 2) {
        return IntVect(2,4,2);
    } else if constexpr (i == 3) {
        return IntVect(4,4,2);
    } else if constexpr (i == 4) {
        return IntVect(2,2,4);
    } else if constexpr (i == 5) {
        return IntVect(4,2,4);
    } else if constexpr (i == 6) {
        return IntVect(2,4,4);
    } else {
        return IntVect(4,4,4);
    }
}

template <int N>
LUSolver<N,Real>
facediv_init_lusolver (IntVect const& ratio,
                       GpuArray<Real, AMREX_SPACEDIM> const& cellSize)
{
    //  set up matrix for local MAC projection
    Array2D<Real,0,N-1,0,N-1,Order::C> mat;

    for (int irow = 0; irow < N; irow++) {
        for (int icol = 0; icol < N; icol++) {
            mat(irow,icol) = Real(0.0);
        }
    }

    Real dxinv = Real(1.)/cellSize[0];
    Real dyinv = Real(1.)/cellSize[1];
    Real dzinv = Real(1.)/cellSize[2];

    Real ddx2 = dxinv*dxinv;
    Real ddy2 = dyinv*dyinv;
    Real ddz2 = dzinv*dzinv;

    for( int k=0; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

        int irow = i + j*ratio[0] +k*ratio[0]*ratio[1];

        if(i > 0){
            mat(irow,irow) += ddx2;
            mat(irow,irow-1) += -ddx2;
        }

        if(i < ratio[0]-1){
            mat(irow,irow) += ddx2;
            mat(irow,irow+1) += -ddx2;
        }

        if(j > 0){
            mat(irow,irow) += ddy2;
            mat(irow,irow-ratio[0]) += -ddy2;
        }

        if(j < ratio[1]-1){
            mat(irow,irow) += ddy2;
            mat(irow,irow+ratio[0]) += -ddy2;
        }

        if(k > 0){
            mat(irow,irow) += ddz2;
            mat(irow,irow-ratio[0]*ratio[1]) += -ddz2;
        }

        if(k < ratio[2]-1){
            mat(irow,irow) += ddz2;
            mat(irow,irow+ratio[0]*ratio[1]) += -ddz2;
        }
    }}}

    //  hcak matrxi to be nonsingular
    mat(N-1,N-1) +=ddx2;

    return LUSolver<N,Real>{mat};
}

template<typename T, int N>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
facediv_int_gen (int ci, int cj, int ck, int nf,
                 GpuArray<Array4<T>, AMREX_SPACEDIM> const& fine,
                 IntVect const& ratio,
                 GpuArray<Real, AMREX_SPACEDIM> const& cellSize,
                 LUSolver<N,Real> const& lusolver) noexcept
{
    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];
    const int fk = ck*ratio[2];

    const Real dx = cellSize[0];
    const Real dy = cellSize[1];
    const Real dz = cellSize[2];

    const Real dxinv = Real(1.)/cellSize[0];
    const Real dyinv = Real(1.)/cellSize[1];
    const Real dzinv = Real(1.)/cellSize[2];

// interpolate initial veloc

    for( int k=0; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=1; i < ratio[0]; ++i){

       fine[0](fi + i, fj + j, fk + k, nf) = ((ratio[0]-i)*fine[0](fi,   fj  , fk  , nf) + i* fine[0](fi+ratio[0],   fj  , fk  , nf) )/ratio[0];

    }
    }
    }

    for( int k=0; k < ratio[2]; ++k){
    for( int j=1; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

       fine[1](fi + i, fj + j, fk + k, nf) = ((ratio[1]-j)*fine[1](fi,   fj  , fk  , nf) + j* fine[1](fi,   fj +ratio[1] , fk  , nf) )/ratio[1];

    }
    }
    }

    for( int k=1; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

       fine[2](fi + i, fj + j, fk + k, nf) = ((ratio[2]-k)*fine[2](fi,   fj  , fk  , nf) + k* fine[2](fi,   fj , fk +ratio[2]  , nf) )/ratio[2];

    }
    }
    }

// compute target divergence

    Real uplus  = 0;
    Real uminus = 0;

    for( int k=0; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){

         uplus  += fine[0](fi+ratio[0],fj+j,fk+k);
         uminus += fine[0](fi         ,fj+j,fk+k);

    }
    }

    Real vplus  = 0;
    Real vminus = 0;

    for( int k=0; k < ratio[2]; ++k){
    for( int i=0; i < ratio[0]; ++i){

         vplus  += fine[1](fi+i,fj+ratio[1],fk+k);
         vminus += fine[1](fi+i,fj,         fk+k);

    }
    }

    Real wplus  = 0;
    Real wminus = 0;

    for( int j=0; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

         wplus  += fine[2](fi+i,fj+j,fk+ratio[2]);
         wminus += fine[2](fi+i,fj+j,fk);

    }
    }

    Real target_div = (uplus-uminus)/(N*dx) + (vplus-vminus)/(N*dy) + (wplus-wminus)/(N*dz);

//  set up matrix for local MAC projection

    Array2D<Real,0,N-1,0,N-1,Order::C> mat;



// form rhs
    Array1D<Real,0,N-1> rhs;
    Array1D<Real,0,N-1> phi;

    for( int k=0; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

       int irow = i + j*ratio[0] +k*ratio[0]*ratio[1];

       rhs(irow) = target_div - ((fine[0](fi+i+1,fj+j,  fk+k,  nf)  -fine[0](fi+i,fj+j,fk+k,nf))*dxinv
                              +  (fine[1](fi+i  ,fj+j+1,fk+k,  nf) - fine[1](fi+i,fj+j,fk+k,nf))*dyinv
                              +  (fine[2](fi+i,  fj+j,  fk+k+1,nf) - fine[2](fi+i,fj+j,fk+k,nf))*dzinv);
    }
    }
    }

//  solver the system

    lusolver(phi.begin(), rhs.begin());

//  this is for sanity check
//  can be put in if needed

//    Array3D<Real,0,ratio[0]-1,0,ratio[1]-1,0,ratio[2]-1> phi_reshape;

//    for( int i=0; i < ratio[0]  , ++j){
//    for( int j=0; j < ratio[1]  , ++j){
//    for( int k=0; k < ratio[2]  , ++k){

//       int irow = i + j*ratio[0] +k*ratio[0]*ratio[1];

//       phi_reshape(i,j,k) = phi(irow)

//    }
//    }
//    }


    for( int k=0; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=1; i < ratio[0]; ++i){

       int plus  = i     + j*ratio[0] +k*ratio[0]*ratio[1];
       int minus = i - 1 + j*ratio[0] +k*ratio[0]*ratio[1];
//       fine[0](fi + i, fj + j, fk + k, nf) += ( phi_reshape(i,j,k) - phi_reshape(i-1,j,k))*dxinv
       fine[0](fi + i, fj + j, fk + k, nf) += ( phi(plus) - phi(minus))*dxinv;

    }
    }
    }

    for( int k=0; k < ratio[2]; ++k){
    for( int j=1; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

       int plus  = i  +     j*ratio[0] +k*ratio[0]*ratio[1];
       int minus = i  + (j-1)*ratio[0] +k*ratio[0]*ratio[1];
//       fine[1](fi + i, fj + j, fk + k, nf) += ( phi_reshape(i,j,k) - phi_reshape(i,j-1,k))*dyinv
       fine[1](fi + i, fj + j, fk + k, nf) += ( phi(plus) - phi(minus))*dyinv;

    }
    }
    }

    for( int k=1; k < ratio[2]; ++k){
    for( int j=0; j < ratio[1]; ++j){
    for( int i=0; i < ratio[0]; ++i){

       int plus  = i  + j*ratio[0] + k    *ratio[0]*ratio[1];
       int minus = i  + j*ratio[0] + (k-1)*ratio[0]*ratio[1];
//       fine[2](fi + i, fj + j, fk + k, nf) += ( phi_reshape(i,j,k) - phi_reshape(i,j,k-1))*dzinv
       fine[2](fi + i, fj + j, fk + k, nf) += ( phi(plus) - phi(minus))*dzinv;

    }
    }
    }


}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (i-ii*ratio[0] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/Real(ratio[0]));
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii+1,jj,kk,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (j-jj*ratio[1] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/Real(ratio[1]));
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj+1,kk,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_z (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (k-kk*ratio[2] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(k-kk*ratio[2]) * (Real(1.)/Real(ratio[2]));
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj,kk+1,n);
    }
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
                    if ((fine_state(i,j,k,n) + fine(i,j,k,n)) < Real(0.0)) {
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
            Real crseTot = Real(0.0);
            Real SumN = Real(0.0);
            Real SumP = Real(0.0);
            for         (int k = klo; k <= khi; ++k) {
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        crseTot += fine(i,j,k,n);
                        if (fine_state(i,j,k,n) <= Real(0.0)) {
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
                            if (fine_state(i,j,k,n) < Real(0.0)) {
                                fine(i,j,k,n) = -fine_state(i,j,k,n);
                            }

                            // Then, add the remaining positive proportionally.
                            if (SumP > Real(0.0)) {
                                if (fine_state(i,j,k,n) > Real(0.0)) {
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
                            if (fine_state(i,j,k,n) < Real(0.0)) {
                                // Add correction to negative cells proportionally.
                                fine(i,j,k,n) = alpha * std::abs(fine_state(i,j,k,n));
                            } else {
                                // Don't touch the positive states.
                                fine(i,j,k,n) = Real(0.0);
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
                        ((SumP+SumN+crseTot) > Real(0.0)) )  {

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

                            if (fine_state(i,j,k,n) < Real(0.0)) {
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
                        ((SumP+SumN+crseTot) < Real(0.0)) )  {
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

                            if (fine_state(i,j,k,n) > Real(0.0)) {
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
                fine(i,j,k,0) = Real(0.0);
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
