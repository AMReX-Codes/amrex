#ifndef AMREX_EB_MULTIFAB_UTIL_3D_C_H_
#define AMREX_EB_MULTIFAB_UTIL_3D_C_H_
#include <AMReX_Config.H>
#include <AMReX_REAL.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
EB_interp_in_quad(Real xint,Real yint,Real v0,Real v1,Real v2,Real v3,
                  Real x0,Real y0, Real x1,Real y1,
                  Real x2,Real y2, Real x3,Real y3)
{
//
//      2D Stencil
//
//      v3 @ (x3,y3)    v2  @ (x2,y2)
//
//             (xint,yint)
//
//      v0 @ (x0,y0)    v1  @ (x1,y1)

        Real va = v3 - v0;
        Real vb = v2 - v0;
        Real vc = v1 - v0;

        Real xa = x3 - x0;
        Real xb = x2 - x0;
        Real xc = x1 - x0;

        Real ya = y3 - y0;
        Real yb = y2 - y0;
        Real yc = y1 - y0;

        Real xder  = vc*(-xa + xb)*ya*yb + vb*(xa - xc)*ya*yc + va*(-xb + xc)*yb*yc;
        Real yder  = vc*xa*xb*(ya - yb) + va*xb*xc*(yb - yc) + vb*xa*xc*(-ya + yc);
        Real xyder = -(vc*xb*ya) + vb*xc*ya + vc*xa*yb - va*xc*yb - vb*xa*yc + va*xb*yc ;

        Real det = -(xa*xc*ya*yb) + xb*xc*ya*yb + xa*xb*ya*yc - xb*xc*ya*yc - xa*xb*yb*yc + xa*xc*yb*yc ;

        xder  =  xder / det;
        yder  =  yder / det;
        xyder = xyder / det;

        Real value = v0 + xder*(xint-x0)+yder*(yint-y0)+ xyder*(xint-x0)*(yint-y0);

        return value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real v)
{
    if (f(i-1,j-1,k-1).isCovered() && f(i  ,j-1,k-1).isCovered() &&
        f(i-1,j  ,k-1).isCovered() && f(i  ,j  ,k-1).isCovered() &&
        f(i-1,j-1,k  ).isCovered() && f(i  ,j-1,k  ).isCovered() &&
        f(i-1,j  ,k  ).isCovered() && f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real const * AMREX_RESTRICT v)
{
    if (f(i-1,j-1,k-1).isCovered() && f(i  ,j-1,k-1).isCovered() &&
        f(i-1,j  ,k-1).isCovered() && f(i  ,j  ,k-1).isCovered() &&
        f(i-1,j-1,k  ).isCovered() && f(i  ,j-1,k  ).isCovered() &&
        f(i-1,j  ,k  ).isCovered() && f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_with_vol (int i, int j, int k,
                          Array4<Real const> const& fine, int fcomp,
                          Array4<Real> const& crse, int ccomp,
                          Array4<Real const> const& fv, Array4<Real const> const& vfrc,
                          Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real cv = Real(0.0);
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = fv(ii,jj,kk)*vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown (int i, int j, int k,
                 Array4<Real const> const& fine, int fcomp,
                 Array4<Real> const& crse, int ccomp,
                 Array4<Real const> const& vfrc,
                 Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real cv = Real(0.0);
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_x (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int ii = i*ratio.x;
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real a = Real(0.0);
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(ii,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_y (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int jj = j*ratio.y;
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real a = Real(0.0);
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,jj,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_z (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int kk = k*ratio.z;
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real a = Real(0.0);
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_boundaries (int i, int j, int k,
                            Array4<Real const> const& fine, int fcomp,
                            Array4<Real> const& crse, int ccomp,
                            Array4<Real const> const& ba,
                            Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = Real(0.0);
        Real cv = Real(0.0);
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = ba(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > Real(1.e-30)) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = Real(0.0);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_compute_divergence (int i, int j, int k, int n, Array4<Real> const& divu,
                            Array4<Real const> const& u, Array4<Real const> const& v,
                            Array4<Real const> const& w, Array4<int const> const& ccm,
                            Array4<EBCellFlag const> const& flag, Array4<Real const> const& vfrc,
                            Array4<Real const> const& apx, Array4<Real const> const& apy,
                            Array4<Real const> const& apz, Array4<Real const> const& fcx,
                            Array4<Real const> const& fcy, Array4<Real const> const& fcz,
                            GpuArray<Real,3> const& dxinv, bool already_on_centroids)
{
    if (flag(i,j,k).isCovered())
    {
        divu(i,j,k,n) = Real(0.0);
    }
    else if (flag(i,j,k).isRegular())
    {
        divu(i,j,k,n) = dxinv[0] * (u(i+1,j,k,n)-u(i,j,k,n))
            +           dxinv[1] * (v(i,j+1,k,n)-v(i,j,k,n))
            +           dxinv[2] * (w(i,j,k+1,n)-w(i,j,k,n));
    }
    else if (already_on_centroids)
    {
        divu(i,j,k,n) = (Real(1.0)/vfrc(i,j,k)) * (
                        dxinv[0] * (apx(i+1,j,k)*u(i+1,j,k,n)-apx(i,j,k)*u(i,j,k,n))
            +           dxinv[1] * (apy(i,j+1,k)*v(i,j+1,k,n)-apy(i,j,k)*v(i,j,k,n))
            +           dxinv[2] * (apz(i,j,k+1)*w(i,j,k+1,n)-apz(i,j,k)*w(i,j,k,n)) );
    }
    else
    {
        Real fxm = u(i,j,k,n);
        if (apx(i,j,k) != Real(0.0) && apx(i,j,k) != Real(1.0)) {
            int jj = j + static_cast<int>(std::copysign(Real(1.0), fcx(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(Real(1.0), fcx(i,j,k,1)));
            Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : Real(0.0);
            Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : Real(0.0);
            fxm = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxm
                +      fracy *(Real(1.0)-fracz)*u(i,jj,k ,n)
                +      fracz *(Real(1.0)-fracy)*u(i,j ,kk,n)
                +      fracy *     fracz *u(i,jj,kk,n);
        }

        Real fxp = u(i+1,j,k,n);
        if (apx(i+1,j,k) != Real(0.0) && apx(i+1,j,k) != Real(1.0)) {
            int jj = j + static_cast<int>(std::copysign(Real(1.0),fcx(i+1,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(Real(1.0),fcx(i+1,j,k,1)));
            Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k,0)) : Real(0.0);
            Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? std::abs(fcx(i+1,j,k,1)) : Real(0.0);
            fxp = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxp
                +      fracy *(Real(1.0)-fracz)*u(i+1,jj,k ,n)
                +      fracz *(Real(1.0)-fracy)*u(i+1,j ,kk,n)
                +      fracy *     fracz *u(i+1,jj,kk,n);
        }

        Real fym = v(i,j,k,n);
        if (apy(i,j,k) != Real(0.0) && apy(i,j,k) != Real(1.0)) {
            int ii = i + static_cast<int>(std::copysign(Real(1.0),fcy(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(Real(1.0),fcy(i,j,k,1)));
            Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : Real(0.0);
            Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : Real(0.0);
            fym = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fym
                +      fracx *(Real(1.0)-fracz)*v(ii,j,k ,n)
                +      fracz *(Real(1.0)-fracx)*v(i ,j,kk,n)
                +      fracx *     fracz *v(ii,j,kk,n);
        }

        Real fyp = v(i,j+1,k,n);
        if (apy(i,j+1,k) != Real(0.0) && apy(i,j+1,k) != Real(1.0)) {
            int ii = i + static_cast<int>(std::copysign(Real(1.0),fcy(i,j+1,k,0)));
            int kk = k + static_cast<int>(std::copysign(Real(1.0),fcy(i,j+1,k,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k,0)) : Real(0.0);
            Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? std::abs(fcy(i,j+1,k,1)) : Real(0.0);
            fyp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fyp
                +      fracx *(Real(1.0)-fracz)*v(ii,j+1,k ,n)
                +      fracz *(Real(1.0)-fracx)*v(i ,j+1,kk,n)
                +      fracx *     fracz *v(ii,j+1,kk,n);
        }

        Real fzm = w(i,j,k,n);
        if (apz(i,j,k) != Real(0.0) && apz(i,j,k) != Real(1.0)) {
            int ii = i + static_cast<int>(std::copysign(Real(1.0),fcz(i,j,k,0)));
            int jj = j + static_cast<int>(std::copysign(Real(1.0),fcz(i,j,k,1)));
            Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : Real(0.0);
            Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : Real(0.0);
            fzm = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzm
                +      fracx *(Real(1.0)-fracy)*w(ii,j ,k,n)
                +      fracy *(Real(1.0)-fracx)*w(i ,jj,k,n)
                +      fracx *     fracy *w(ii,jj,k,n);
        }

        Real fzp = w(i,j,k+1,n);
        if (apz(i,j,k+1) != Real(0.0) && apz(i,j,k+1) != Real(1.0)) {
            int ii = i + static_cast<int>(std::copysign(Real(1.0),fcz(i,j,k+1,0)));
            int jj = j + static_cast<int>(std::copysign(Real(1.0),fcz(i,j,k+1,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? std::abs(fcz(i,j,k+1,0)) : Real(0.0);
            Real fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? std::abs(fcz(i,j,k+1,1)) : Real(0.0);
            fzp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzp
                +      fracx *(Real(1.0)-fracy)*w(ii,j ,k+1,n)
                +      fracy *(Real(1.0)-fracx)*w(i ,jj,k+1,n)
                +      fracx *     fracy *w(ii,jj,k+1,n);
        }

        divu(i,j,k,n) = (Real(1.0)/vfrc(i,j,k)) *
            ( dxinv[0] * (apx(i+1,j,k)*fxp-apx(i,j,k)*fxm)
            + dxinv[1] * (apy(i,j+1,k)*fyp-apy(i,j,k)*fym)
            + dxinv[2] * (apz(i,j,k+1)*fzp-apz(i,j,k)*fzm) );
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avg_fc_to_cc (int i, int j, int k, int n, Array4<Real> const& cc,
                      Array4<Real const> const& fx, Array4<Real const> const& fy,
                      Array4<Real const> const& fz, Array4<Real const> const& ax,
                      Array4<Real const> const& ay, Array4<Real const> const& az,
                      Array4<EBCellFlag const> const& flag)
{
    if (flag(i,j,k).isCovered()) {
        cc(i,j,k,n+0) = Real(0.0);
        cc(i,j,k,n+1) = Real(0.0);
        cc(i,j,k,n+2) = Real(0.0);
    } else {
        if (ax(i,j,k) == Real(0.0)) {
            cc(i,j,k,n+0) = fx(i+1,j,k);
        } else if (ax(i+1,j,k) == Real(0.0)) {
            cc(i,j,k,n+0) = fx(i,j,k);
        } else {
            cc(i,j,k,n+0) = Real(0.5) * (fx(i,j,k) + fx(i+1,j,k));
        }

        if (ay(i,j,k) == Real(0.0)) {
            cc(i,j,k,n+1) = fy(i,j+1,k);
        } else if (ay(i,j+1,k) == Real(0.0)) {
            cc(i,j,k,n+1) = fy(i,j,k);
        } else {
            cc(i,j,k,n+1) = Real(0.5) * (fy(i,j,k) + fy(i,j+1,k));
        }

        if (az(i,j,k) == Real(0.0)) {
            cc(i,j,k,n+2) = fz(i,j,k+1);
        } else if (az(i,j,k+1) == Real(0.0)) {
            cc(i,j,k,n+2) = fz(i,j,k);
        } else {
            cc(i,j,k,n+2) = Real(0.5) * (fz(i,j,k) + fz(i,j,k+1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2cent (Box const& box,
                        const Array4<Real>& phicent,
                        Array4<Real const > const& phicc,
                        Array4<EBCellFlag const> const& flag,
                        Array4<Real const> const& cent,
                        int ncomp) noexcept
{
  amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (flag(i,j,k).isCovered())
    {
      phicent(i,j,k,n) = phicc(i,j,k,n); //do nothing
    }
    else
    {
      if (flag(i,j,k).isRegular())
      {
        phicent(i,j,k,n) = phicc(i,j,k,n);
      }
      else
      {
        Real gx = cent(i,j,k,0);
        Real gy = cent(i,j,k,1);
        Real gz = cent(i,j,k,2);
        int ii = (gx < Real(0.0)) ? i - 1 : i + 1;
        int jj = (gy < Real(0.0)) ? j - 1 : j + 1;
        int kk = (gz < Real(0.0)) ? k - 1 : k + 1;
        gx = std::abs(gx);
        gy = std::abs(gy);
        gz = std::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        phicent(i,j,k,n)
          = ( Real(1.0) - gx - gy - gz + gxy + gxz + gyz - gxyz) * phicc(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phicc(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phicc(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phicc(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phicc(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phicc(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phicc(ii,jj,k ,n)
          + (                                        gxyz) * phicc(ii,jj,kk,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2facecent_x (Box const& ubx,
                              Array4<Real const> const& phi,
                              Array4<Real const> const& apx,
                              Array4<Real const> const& fcx,
                              Array4<Real> const& edg_x,
                              int ncomp,
                              const Box& domain,
                              const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apx(i,j,k) == 0)
    {
#ifdef AMREX_USE_FLOAT
      edg_x(i,j,k,n) = Real(1e20);
#else
      edg_x(i,j,k,n) = 1e40;
#endif
    }
    else
    {
      if (apx(i,j,k) == Real(1.))
      {
        if ( (i == domlo.x) && (bc[n].lo(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
        }
        else if ( (i == domhi.x+1) && (bc[n].hi(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
        }
        else
        {
          edg_x(i,j,k,n) = Real(0.5) * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
        }
      }
      else
      {
        Real gx = Real(0.5);
        Real gy = fcx(i,j,k,0);
        Real gz = fcx(i,j,k,1);
        int ii = i - 1;
        int jj = (gy < Real(0.0)) ? j - 1 : j + 1;
        int kk = (gz < Real(0.0)) ? k - 1 : k + 1;
        gy = std::abs(gy);
        gz = std::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_x(i,j,k,n) =
          (   Real(1.0) - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2facecent_y (Box const& vbx,
                              Array4<Real const> const& phi,
                              Array4<Real const> const& apy,
                              Array4<Real const> const& fcy,
                              Array4<Real> const& edg_y,
                              int ncomp,
                              const Box& domain,
                              const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apy(i,j,k) == 0)
    {
#ifdef AMREX_USE_FLOAT
      edg_y(i,j,k,n) = Real(1e20);
#else
      edg_y(i,j,k,n) = 1e40;
#endif
    }
    else
    {
      if (apy(i,j,k) == Real(1.))
      {
        if ( (j == domlo.y) && (bc[n].lo(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) && (bc[n].hi(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else
        {
          edg_y(i,j,k,n) = Real(0.5) * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
      }
      else
      {
        Real gx = fcy(i,j,k,0);
        Real gy = Real(0.5);
        Real gz = fcy(i,j,k,1);
        int ii = (gx < Real(0.0)) ? i - 1 : i + 1;
        int jj = j - 1;
        int kk = (gz < Real(0.0)) ? k - 1 : k + 1;
        gx = std::abs(gx);
        gz = std::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_y(i,j,k,n) =
          (   Real(1.0) - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2facecent_z (Box const& wbx,
                              Array4<Real const> const& phi,
                              Array4<Real const> const& apz,
                              Array4<Real const> const& fcz,
                              Array4<Real> const& edg_z,
                              int ncomp,
                              const Box& domain,
                              const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apz(i,j,k) == 0)
    {
#ifdef AMREX_USE_FLOAT
      edg_z(i,j,k,n) = Real(1e20);
#else
      edg_z(i,j,k,n) = 1e40;
#endif
    }
    else
    {
      if (apz(i,j,k) == Real(1.))
      {
        if ( (k == domlo.z) && (bc[n].lo(2) == BCType::ext_dir) )
        {
          edg_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
        }
        else if ( (k == domhi.z+1) && (bc[n].hi(2) == BCType::ext_dir) )
        {
          edg_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
        }
        else
        {
          edg_z(i,j,k,n) = Real(0.5) * (phi(i,j  ,k,n) + phi(i,j,k-1,n));
        }
      }
      else
      {
        Real gx = fcz(i,j,k,0);
        Real gy = fcz(i,j,k,1);
        Real gz = Real(0.5);
        int ii = (gx < Real(0.0)) ? i - 1 : i + 1;
        int jj = (gy < Real(0.0)) ? j - 1 : j + 1;
        int kk = k -1;
        gx = std::abs(gx);
        gy = std::abs(gy);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_z(i,j,k,n) =
          (   Real(1.0) - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_centroid2facecent_x (Box const& ubx,
                                    Array4<Real const> const& phi,
                                    Array4<Real const> const& apx,
                                    Array4<Real const> const& cvol,
                                    Array4<Real const> const& ccent,
                                    Array4<Real const> const& fcx,
                                    Array4<Real> const& phi_x,
                                    int ncomp,
                                    const Box& domain,
                                    const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  // Note that ccent holds (x,y,z) of the cell centroids  as components (0/1/2)
  //           fcx   holds (  y,z) of the x-face centroid as components ( /0/1)
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
      if (apx(i,j,k) == 0)
      {
#ifdef AMREX_USE_FLOAT
          phi_x(i,j,k,n) = Real(1e20);
#else
          phi_x(i,j,k,n) = 1e40;
#endif
      }
      else if ( (i == domlo.x) && (bc[n].lo(0) == BCType::ext_dir) )
      {
          phi_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
      }
      else if ( (i == domhi.x+1) && (bc[n].hi(0) == BCType::ext_dir) )
      {
          phi_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
      }
      else if (apx(i,j,k) == Real(1.) && cvol(i,j,k) == Real(1.) && cvol(i-1,j,k) == Real(1.))
      {
          phi_x(i,j,k,n) = Real(0.5) * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
      }
      else if ( apx(i,j,k) == Real(1.) && (std::abs(ccent(i,j,k,1)-ccent(i-1,j,k,1)) < Real(1.e-8)) &&
                                   (std::abs(ccent(i,j,k,2)-ccent(i-1,j,k,2)) < Real(1.e-8)) )
      {
          Real d0 = Real(0.5) + ccent(i  ,j,k,0);                               // distance in x-dir from centroid(i  ,j) to face(i,j)
          Real d1 = Real(0.5) - ccent(i-1,j,k,0);                               // distance in x-dir from centroid(i-1,j) to face(i,j)
          Real a0 = d1 / (d0 + d1);                                       // coefficient of phi(i  ,j,k,n)
          Real a1 = d0 / (d0 + d1);                                       // coefficient of phi(i-1,j,k,n)
          phi_x(i,j,k,n) = a0*phi(i,j,k,n) + a1*phi(i-1,j,k,n);           // interpolation of phi in x-direction
      }
      else
      {
          // We must add additional tests to avoid the case where fcx is very close to zero, but j-1/j+1 or k-1/k+1
          //    might be covered cells -- this can happen when the EB is exactly aligned with the horizontal grid planes
          int jj,kk;
          if (std::abs(fcx(i,j,k,0)) > Real(1.e-8))
              jj = (fcx(i,j,k,0) < Real(0.0)) ? j - 1 : j + 1;
          else if (apx(i,j-1,k) > Real(0.))
              jj = j-1;
          else
              jj = j+1;

          if (std::abs(fcx(i,j,k,1)) > Real(1.e-8))
              kk = (fcx(i,j,k,1) < Real(0.0)) ? k - 1 : k + 1;
          else if (apx(i,j,k-1) > Real(0.))
              kk = k-1;
          else
              kk = k+1;

          // If any of these cells has zero volume we don't want to use this stencil
          Real test_zero =  cvol(i-1,jj,k  ) * cvol(i-1,j  ,kk) *  cvol(i-1,jj,kk) *
                            cvol(i  ,jj,k  ) * cvol(i  ,j  ,kk) *  cvol(i  ,jj,kk);

          if (test_zero ==  Real(0.))
          {
              int jalt = 2*j - jj;
              int kalt = 2*k - kk;
              Real test_zero_jalt  = cvol(i-1,jalt,k) * cvol(i-1,j,kk  ) *  cvol(i-1,jalt,kk  ) *
                                     cvol(i  ,jalt,k) * cvol(i  ,j,kk  ) *  cvol(i  ,jalt,kk  );

              Real test_zero_kalt  = cvol(i-1,jj  ,k) * cvol(i-1,j,kalt) *  cvol(i-1,jj  ,kalt) *
                                     cvol(i  ,jj  ,k) * cvol(i  ,j,kalt) *  cvol(i  ,jj  ,kalt);

              Real test_zero_jkalt = cvol(i-1,jalt,k) * cvol(i-1,j,kalt) *  cvol(i-1,jalt,kalt) *
                                     cvol(i  ,jalt,k) * cvol(i  ,j,kalt) *  cvol(i  ,jalt,kalt);

              if (test_zero_jalt > Real(0.))
                 jj = jalt;
              else if (test_zero_kalt > Real(0.))
                 kk = kalt;
              else if (test_zero_jkalt > Real(0.))
              {
                 jj = jalt;
                 kk = kalt;
              }
          }

          Real d0       = Real(0.5) + ccent(i  , j, k,0);                     // distance in x-dir from centroid(i  , j, k) to face(i,j,k)
          Real d1       = Real(0.5) - ccent(i-1, j, k,0);                     // distance in x-dir from centroid(i-1, j, k) to face(i,j,k)

          Real d0_jj    = Real(0.5) + ccent(i  ,jj, k,0);                     // distance in x-dir from centroid(i  ,jj, k) to face(i,jj,k)
          Real d1_jj    = Real(0.5) - ccent(i-1,jj, k,0);                     // distance in x-dir from centroid(i-1,jj, k) to face(i,jj,k)

          Real d0_kk    = Real(0.5) + ccent(i  , j,kk,0);                     // distance in x-dir from centroid(i  , j,kk) to face(i,j,kk)
          Real d1_kk    = Real(0.5) - ccent(i-1, j,kk,0);                     // distance in x-dir from centroid(i-1, j,kk) to face(i,j,kk)

          Real d0_jj_kk = Real(0.5) + ccent(i  ,jj,kk,0);                     // distance in x-dir from centroid(i  ,jj,kk) to face(i,jj,kk)
          Real d1_jj_kk = Real(0.5) - ccent(i-1,jj,kk,0);                     // distance in x-dir from centroid(i-1,jj,kk) to face(i,jj,kk)

          Real a0       = d1 / (d0 + d1);                             // coefficient of phi(i  , j, k,n)
          Real a1       = Real(1.0) - a0;                                   // coefficient of phi(i-1, j, k,n)

          Real a0_jj    = d1_jj / (d0_jj + d1_jj);                    // coefficient of phi(i  ,jj, k,n)
          Real a1_jj    = Real(1.0) - a0_jj;                                // coefficient of phi(i-1,jj, k,n)

          Real a0_kk    = d1_kk / (d0_kk + d1_kk);                    // coefficient of phi(i  , j,kk,n)
          Real a1_kk    = Real(1.0) - a0_kk;                                // coefficient of phi(i-1, j,kk,n)

          Real a0_jj_kk = d1_jj_kk / (d0_jj_kk + d1_jj_kk);           // coefficient of phi(i  ,jj,kk,n)
          Real a1_jj_kk = Real(1.0) - a0_jj_kk;                             // coefficient of phi(i-1,jj,kk,n)

          Real phi01       = a0       *   phi(i, j, k,n)  + a1       *   phi(i-1, j, k,n); // interpolation in x-dir of phi
          Real y01         = a0       * ccent(i, j, k,1)  + a1       * ccent(i-1, j, k,1); // interpolation in x-dir of y-loc
          Real z01         = a0       * ccent(i, j, k,2)  + a1       * ccent(i-1, j, k,2); // interpolation in x-dir of z-loc

          Real phi01_jj    = a0_jj    *   phi(i,jj, k,n)  + a1_jj    *   phi(i-1,jj, k,n); // interpolation in x-dir of phi
          Real y01_jj      = a0_jj    * ccent(i,jj, k,1)  + a1_jj    * ccent(i-1,jj, k,1); // interpolation in x-dir of y-loc
          Real z01_jj      = a0_jj    * ccent(i,jj, k,2)  + a1_jj    * ccent(i-1,jj, k,2); // interpolation in x-dir of z-loc

          Real phi01_kk    = a0_kk    *   phi(i, j,kk,n)  + a1_kk    *   phi(i-1, j,kk,n); // interpolation in x-dir of phi
          Real y01_kk      = a0_kk    * ccent(i, j,kk,1)  + a1_kk    * ccent(i-1, j,kk,1); // interpolation in x-dir of y-loc
          Real z01_kk      = a0_kk    * ccent(i, j,kk,2)  + a1_kk    * ccent(i-1, j,kk,2); // interpolation in x-dir of z-loc

          Real phi01_jj_kk = a0_jj_kk *   phi(i,jj,kk,n)  + a1_jj_kk *   phi(i-1,jj,kk,n); // interpolation in x-dir of phi
          Real y01_jj_kk   = a0_jj_kk * ccent(i,jj,kk,1)  + a1_jj_kk * ccent(i-1,jj,kk,1); // interpolation in x-dir of y-loc
          Real z01_jj_kk   = a0_jj_kk * ccent(i,jj,kk,2)  + a1_jj_kk * ccent(i-1,jj,kk,2); // interpolation in x-dir of z-loc

          // 2D interpolation on x-face from interpolated points at (y01,z01), (y01_jj,z01_jj), (y01_kk,z01_kk), (y01_jj_kk,z01_jj_kk)
          // to centroid of x-face(i,j,k) at (fcx(i,j,k,0), fcx(i,j,k,1))

          // We order these in order to form a set of four points on the x-face ...
          // (x0,y0) : lower left
          // (x1,y1) : lower right
          // (x2,y2) : upper right
          // (x3,y3) : upper left

          Real y,z;

          // This is the location of the face centroid relative to the central node
          // Recall fcx holds (y,z) of the x-face centroid as components ( /0/1)
          if (j < jj)
             y = Real(-0.5) + fcx(i,j,k,0);  // (j,k) is in lower half of stencil so y < 0
          else
             y =  Real(0.5) + fcx(i,j,k,0);  // (j,k) is in upper half of stencil so y > 0

          if (k < kk)
             z = Real(-0.5) + fcx(i,j,k,1);  // (j,k) is in lower half of stencil so z < 0
          else
             z =  Real(0.5) + fcx(i,j,k,1);  // (j,k) is in upper half of stencil so z > 0

          if (j < jj && k < kk) // (j,k) is lower left, (j+1,k+1) is upper right
          {
              Real y_ll = Real(-0.5)+y01;
              Real z_ll = Real(-0.5)+z01;
              Real y_hl =  Real(0.5)+y01_jj;
              Real z_hl = Real(-0.5)+z01_jj;
              Real y_hh =  Real(0.5)+y01_jj_kk;
              Real z_hh =  Real(0.5)+z01_jj_kk;
              Real y_lh = Real(-0.5)+y01_kk;
              Real z_lh =  Real(0.5)+z01_kk;
              phi_x(i,j,k,n) = EB_interp_in_quad(y,z, phi01, phi01_jj, phi01_jj_kk, phi01_kk,
                                                 y_ll, z_ll, y_hl, z_hl, y_hh, z_hh, y_lh, z_lh);
          }
          else if (jj < j && k < kk) // (j-1,k) is lower left, (j,k+1) is upper right
          {
              Real y_ll = Real(-0.5)+y01_jj;
              Real z_ll = Real(-0.5)+z01_jj;
              Real y_hl =  Real(0.5)+y01;
              Real z_hl = Real(-0.5)+z01;
              Real y_hh =  Real(0.5)+y01_kk;
              Real z_hh =  Real(0.5)+z01_kk;
              Real y_lh = Real(-0.5)+y01_jj_kk;
              Real z_lh =  Real(0.5)+z01_jj_kk;
              phi_x(i,j,k,n) = EB_interp_in_quad(y,z,phi01_jj,phi01,phi01_kk, phi01_jj_kk,
                                                 y_ll, z_ll, y_hl, z_hl, y_hh, z_hh, y_lh, z_lh);
          }
          else if (j < jj && kk < k) // (j,k-1) is lower left, (j+1,k) is upper right
          {
              Real y_ll = Real(-0.5)+y01_kk;
              Real z_ll = Real(-0.5)+z01_kk;
              Real y_hl =  Real(0.5)+y01_jj_kk;
              Real z_hl = Real(-0.5)+z01_jj_kk;
              Real y_hh =  Real(0.5)+y01_jj;
              Real z_hh =  Real(0.5)+z01_jj;
              Real y_lh = Real(-0.5)+y01;
              Real z_lh =  Real(0.5)+z01;
              phi_x(i,j,k,n) = EB_interp_in_quad(y,z,phi01_kk,phi01_jj_kk,phi01_jj, phi01,
                                                 y_ll, z_ll, y_hl, z_hl, y_hh, z_hh, y_lh, z_lh);
          }
          else if (jj < j && kk < k) // (j-1,k-1) is lower left, (j,k) is upper right
          {
              Real y_ll = Real(-0.5)+y01_jj_kk;
              Real z_ll = Real(-0.5)+z01_jj_kk;
              Real y_hl =  Real(0.5)+y01_kk;
              Real z_hl = Real(-0.5)+z01_kk;
              Real y_hh =  Real(0.5)+y01;
              Real z_hh =  Real(0.5)+z01;
              Real y_lh = Real(-0.5)+y01_jj;
              Real z_lh =  Real(0.5)+z01_jj;
              phi_x(i,j,k,n) = EB_interp_in_quad(y,z,phi01_jj_kk,phi01_kk,phi01, phi01_jj,
                                                 y_ll, z_ll, y_hl, z_hl, y_hh, z_hh, y_lh, z_lh);
          }
          else
          {
              amrex::Abort("Bad option in interpolation from cell centroid to x-face centroid!");
          }
       }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_centroid2facecent_y (Box const& vbx,
                                    Array4<Real const> const& phi,
                                    Array4<Real const> const& apy,
                                    Array4<Real const> const& cvol,
                                    Array4<Real const> const& ccent,
                                    Array4<Real const> const& fcy,
                                    Array4<Real> const& phi_y,
                                    int ncomp,
                                    const Box& domain,
                                    const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  // Note that ccent holds (x,y,z) of the cell centroids  as components (0/1/2)
  //           fcy   holds (x,  z) of the y-face centroid as components (0/ /1)
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == 0)
        {
#ifdef AMREX_USE_FLOAT
            phi_y(i,j,k,n) = Real(1e20);
#else
            phi_y(i,j,k,n) = 1e40;
#endif
        }
        else if ( (j == domlo.y) && (bc[n].lo(1) == BCType::ext_dir) )
        {
            phi_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) && (bc[n].hi(1) == BCType::ext_dir) )
        {
            phi_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else if (apy(i,j,k) == Real(1.) && cvol(i,j,k) == Real(1.) && cvol(i,j-1,k) == Real(1.))
        {
            phi_y(i,j,k,n) = Real(0.5) * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
        else if ( apy(i,j,k) == Real(1.) && (std::abs(ccent(i,j,k,0)-ccent(i,j-1,k,0)) < Real(1.e-8)) &&
                                      (std::abs(ccent(i,j,k,2)-ccent(i,j-1,k,2)) < Real(1.e-8)) )
        {
            Real d0 = Real(0.5) + ccent(i,j  ,k,1);                                // distance in y-dir from centroid(i,j  ,k) to face(i,j,k)
            Real d1 = Real(0.5) - ccent(i,j-1,k,1);                                // distance in y-dir from centroid(i,j-1,k) to face(i,j,k)
            Real a0 = d1 / (d0 + d1);                                        // coefficient of phi(i,j  ,k,n)
            Real a1 = d0 / (d0 + d1);                                        // coefficient of phi(i,j-1,k,n)
            phi_y(i,j,k,n) = a0*phi(i,j,k,n) + a1*phi(i,j-1,k,n);            // interpolation of phi in y-direction
        }
        else
        {

          // We must add additional tests to avoid the case where fcy is very close to zero, but i-1/i+1 or k-1/k+1
          //    might be covered cells -- this can happen when the EB is exactly aligned with the grid planes
          int ii,kk;
          if (std::abs(fcy(i,j,k,0)) > Real(1.e-8))
              ii = (fcy(i,j,k,0) < Real(0.0)) ? i - 1 : i + 1;
          else if (apy(i-1,j,k) > Real(0.))
              ii = i-1;
          else
              ii = i+1;

          if (std::abs(fcy(i,j,k,1)) > Real(1.e-8))
              kk = (fcy(i,j,k,1) < Real(0.0)) ? k - 1 : k + 1;
          else if (apy(i,j,k-1) > Real(0.))
              kk = k-1;
          else
              kk = k+1;

          // If any of these cells has zero volume we don't want to use this stencil
          Real test_zero =  cvol(ii,j-1,k) * cvol(i,j-1,kk) *  cvol(ii,j-1,kk) *
                            cvol(ii,j  ,k) * cvol(i,j  ,kk) *  cvol(ii,j  ,kk);

          if (test_zero ==  Real(0.))
          {
              int ialt = 2*i - ii;
              int kalt = 2*k - kk;
              Real test_zero_ialt  = cvol(ialt,j-1,k) * cvol(i,j-1,kk  ) *  cvol(ialt,j-1,kk  ) *
                                     cvol(ialt,j  ,k) * cvol(i,j  ,kk  ) *  cvol(ialt,j  ,kk  );

              Real test_zero_kalt  = cvol(ii  ,j-1,k) * cvol(i,j-1,kalt) *  cvol(ii  ,j-1,kalt) *
                                     cvol(ii  ,j  ,k) * cvol(i,j  ,kalt) *  cvol(ii  ,j  ,kalt);

              Real test_zero_ikalt = cvol(ialt,j-1,k) * cvol(i,j-1,kalt) *  cvol(ialt,j-1,kalt) *
                                     cvol(ialt,j  ,k) * cvol(i,j  ,kalt) *  cvol(ialt,j  ,kalt);

              if (test_zero_ialt > Real(0.))
                 ii = ialt;
              else if (test_zero_kalt > Real(0.))
                 kk = kalt;
              else if (test_zero_ikalt > Real(0.))
              {
                 ii = ialt;
                 kk = kalt;
              }
          }

          Real d0       = Real(0.5) + ccent( i,j  , k,1); // distance in y-dir from centroid( i,j  , k) to face(i,j,k)
          Real d1       = Real(0.5) - ccent( i,j-1, k,1); // distance in y-dir from centroid( i,j-1, k) to face(i,j,k)

          Real d0_ii    = Real(0.5) + ccent(ii,j  , k,1); // distance in y-dir from centroid(ii,j  , k) to face(ii,j,k)
          Real d1_ii    = Real(0.5) - ccent(ii,j-1, k,1); // distance in y-dir from centroid(ii,j-1, k) to face(ii,j,k)

          Real d0_kk    = Real(0.5) + ccent( i,j  ,kk,1); // distance in y-dir from centroid( i,j  ,kk) to face( i,j,kk)
          Real d1_kk    = Real(0.5) - ccent( i,j-1,kk,1); // distance in y-dir from centroid( i,j-1,kk) to face( i,j,kk)

          Real d0_ii_kk = Real(0.5) + ccent(ii,j  ,kk,1); // distance in y-dir from centroid(ii,j  ,kk) to face(ii,j,kk)
          Real d1_ii_kk = Real(0.5) - ccent(ii,j-1,kk,1); // distance in y-dir from centroid(ii,j-1,kk) to face(ii,j,kk)

          Real a0    = d1 / (d0 + d1);                              // coefficient of phi( i,j  ,k,n)
          Real a1    = d0 / (d0 + d1);                              // coefficient of phi( i,j-1,k,n)

          Real a0_ii    = d1_ii / (d0_ii + d1_ii);                  // coefficient of phi(ii,j  ,k,n)
          Real a1_ii    = d0_ii / (d0_ii + d1_ii);                  // coefficient of phi(ii,j-1,k,n)

          Real a0_kk    = d1_kk / (d0_kk + d1_kk);                  // coefficient of phi( i,j  ,kk,n)
          Real a1_kk    = d0_kk / (d0_kk + d1_kk);                  // coefficient of phi( i,j-1,kk,n)

          Real a0_ii_kk = d1_ii_kk / (d0_ii_kk + d1_ii_kk);         // coefficient of phi(ii,j  ,kk,n)
          Real a1_ii_kk = d0_ii_kk / (d0_ii_kk + d1_ii_kk);         // coefficient of phi(ii,j-1,kk,n)

          Real phi01       = a0       *   phi( i, j, k,n)  + a1       *   phi( i, j-1, k,n); // interpolation in y-dir of phi
          Real x01         = a0       * ccent( i, j, k,0)  + a1       * ccent( i, j-1, k,0); // interpolation in y-dir of x-loc
          Real z01         = a0       * ccent( i, j, k,2)  + a1       * ccent( i, j-1, k,2); // interpolation in y-dir of z-loc

          Real phi01_ii    = a0_ii    *   phi(ii, j, k,n)  + a1_ii    *   phi(ii, j-1, k,n); // interpolation in y-dir of phi
          Real x01_ii      = a0_ii    * ccent(ii, j, k,0)  + a1_ii    * ccent(ii, j-1, k,0); // interpolation in y-dir of x-loc
          Real z01_ii      = a0_ii    * ccent(ii, j, k,2)  + a1_ii    * ccent(ii, j-1, k,2); // interpolation in y-dir of z-loc

          Real phi01_kk    = a0_kk    *   phi( i, j,kk,n)  + a1_kk    *   phi( i, j-1,kk,n); // interpolation in y-dir of phi
          Real x01_kk      = a0_kk    * ccent( i, j,kk,0)  + a1_kk    * ccent( i, j-1,kk,0); // interpolation in y-dir of x-loc
          Real z01_kk      = a0_kk    * ccent( i, j,kk,2)  + a1_kk    * ccent( i, j-1,kk,2); // interpolation in y-dir of z-loc

          Real phi01_ii_kk = a0_ii_kk *   phi(ii, j,kk,n)  + a1_ii_kk *   phi(ii, j-1,kk,n); // interpolation in y-dir of phi
          Real x01_ii_kk   = a0_ii_kk * ccent(ii, j,kk,0)  + a1_ii_kk * ccent(ii, j-1,kk,0); // interpolation in y-dir of x-loc
          Real z01_ii_kk   = a0_ii_kk * ccent(ii, j,kk,2)  + a1_ii_kk * ccent(ii, j-1,kk,2); // interpolation in y-dir of z-loc

          // We order these in order to form a set of four points on the x-face ...
          // (x0,y0) : lower left
          // (x1,y1) : lower right
          // (x2,y2) : upper right
          // (x3,y3) : upper left

          Real x,z;

          // This is the location of the face centroid relative to the central node
          // Recall fcy holds (x,z) of the x-face centroid as components (0/ /1)
          if (i < ii)
             x = Real(-0.5) + fcy(i,j,k,0);  // (i,k) is in lower half of stencil so x < 0
          else
             x =  Real(0.5) + fcy(i,j,k,0);  // (i,k) is in upper half of stencil so x > 0

          if (k < kk)
             z = Real(-0.5) + fcy(i,j,k,1);  // (i,k) is in lower half of stencil so z < 0
          else
             z =  Real(0.5) + fcy(i,j,k,1);  // (i,k) is in upper half of stencil so z > 0

          if (i < ii && k < kk) // (i,k) is lower left, (i+1,k+1) is upper right
          {
              Real x_ll = Real(-0.5)+x01;
              Real z_ll = Real(-0.5)+z01;
              Real x_hl =  Real(0.5)+x01_ii;
              Real z_hl = Real(-0.5)+z01_ii;
              Real x_hh =  Real(0.5)+x01_ii_kk;
              Real z_hh =  Real(0.5)+z01_ii_kk;
              Real x_lh = Real(-0.5)+x01_kk;
              Real z_lh =  Real(0.5)+z01_kk;
              phi_y(i,j,k,n) = EB_interp_in_quad(x,z, phi01, phi01_ii, phi01_ii_kk, phi01_kk,
                                                 x_ll, z_ll, x_hl, z_hl, x_hh, z_hh, x_lh, z_lh);
          }
          else if (ii < i && k < kk) // (i-1,k) is lower left, (i,k+1) is upper right
          {
              Real x_ll = Real(-0.5)+x01_ii;
              Real z_ll = Real(-0.5)+z01_ii;
              Real x_hl =  Real(0.5)+x01;
              Real z_hl = Real(-0.5)+z01;
              Real x_hh =  Real(0.5)+x01_kk;
              Real z_hh =  Real(0.5)+z01_kk;
              Real x_lh = Real(-0.5)+x01_ii_kk;
              Real z_lh =  Real(0.5)+z01_ii_kk;
              phi_y(i,j,k,n) = EB_interp_in_quad(x,z,phi01_ii,phi01,phi01_kk, phi01_ii_kk,
                                                 x_ll, z_ll, x_hl, z_hl, x_hh, z_hh, x_lh, z_lh);
          }
          else if (i < ii && kk < k) // (i,k-1) is lower left, (i+1,k) is upper right
          {
              Real x_ll = Real(-0.5)+x01_kk;
              Real z_ll = Real(-0.5)+z01_kk;
              Real x_hl =  Real(0.5)+x01_ii_kk;
              Real z_hl = Real(-0.5)+z01_ii_kk;
              Real x_hh =  Real(0.5)+x01_ii;
              Real z_hh =  Real(0.5)+z01_ii;
              Real x_lh = Real(-0.5)+x01;
              Real z_lh =  Real(0.5)+z01;
              phi_y(i,j,k,n) = EB_interp_in_quad(x,z,phi01_kk,phi01_ii_kk,phi01_ii, phi01,
                                                 x_ll, z_ll, x_hl, z_hl, x_hh, z_hh, x_lh, z_lh);
          }
          else if (ii < i && kk < k) // (i-1,k-1) is lower left, (i,k) is upper right
          {
              Real x_ll = Real(-0.5)+x01_ii_kk;
              Real z_ll = Real(-0.5)+z01_ii_kk;
              Real x_hl =  Real(0.5)+x01_kk;
              Real z_hl = Real(-0.5)+z01_kk;
              Real x_hh =  Real(0.5)+x01;
              Real z_hh =  Real(0.5)+z01;
              Real x_lh = Real(-0.5)+x01_ii;
              Real z_lh =  Real(0.5)+z01_ii;
              phi_y(i,j,k,n) = EB_interp_in_quad(x,z,phi01_ii_kk,phi01_kk,phi01, phi01_ii,
                                                 x_ll, z_ll, x_hl, z_hl, x_hh, z_hh, x_lh, z_lh);
          }
          else
          {
              amrex::Abort("Bad option in interpolation from cell centroid to y-face centroid!");
          }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_centroid2facecent_z (Box const& wbx,
                                    Array4<Real const> const& phi,
                                    Array4<Real const> const& apz,
                                    Array4<Real const> const& cvol,
                                    Array4<Real const> const& ccent,
                                    Array4<Real const> const& fcz,
                                    Array4<Real> const& phi_z,
                                    int ncomp,
                                    const Box& domain,
                                    const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  // Note that ccent holds (x,y,z) of the cell centroids  as components (0/1/2)
  //           fcz   holds (x,y  ) of the z-face centroid as components (0/1  )
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) == 0)
        {
#ifdef AMREX_USE_FLOAT
            phi_z(i,j,k,n) = Real(1e20);
#else
            phi_z(i,j,k,n) = 1e40;
#endif
        }
        else if ( (k == domlo.z) && (bc[n].lo(2) == BCType::ext_dir) )
        {
            phi_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
        }
        else if ( (k == domhi.z+1) && (bc[n].hi(2) == BCType::ext_dir) )
        {
            phi_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
        }
        else if (apz(i,j,k) == Real(1.) && cvol(i,j,k) == Real(1.) && cvol(i,j,k-1) == Real(1.))
        {
            phi_z(i,j,k,n) = Real(0.5) * (phi(i,j,k,n) + phi(i,j,k-1,n));
        }
        else if ( apz(i,j,k) == Real(1.) && (std::abs(ccent(i,j,k,0)-ccent(i,j,k-1,0)) < Real(1.e-8)) &&
                                      (std::abs(ccent(i,j,k,1)-ccent(i,j,k-1,1)) < Real(1.e-8)) )
        {
            Real d0 = Real(0.5) + ccent(i,j,k  ,2);                                // distance in z-dir from centroid(i,j,k  ) to face(i,j,k)
            Real d1 = Real(0.5) - ccent(i,j,k-1,2);                                // distance in z-dir from centroid(i,j,k-1) to face(i,j,k)
            Real a0 = d1 / (d0 + d1);                                        // coefficient of phi(i,j,k  ,n)
            Real a1 = d0 / (d0 + d1);                                        // coefficient of phi(i,j,k-1,n)
            phi_z(i,j,k,n) = a0*phi(i,j,k,n) + a1*phi(i,j,k-1,n);            // interpolation of phi in z-direction
        }
        else
        {
          // We must add additional tests to avoid the case where fcz is very close to zero, but i-1/i+1 or j-1/j+1
          //    might be covered cells -- this can happen when the EB is exactly aligned with the grid planes
          int ii,jj;
          if (std::abs(fcz(i,j,k,0)) > Real(1.e-8))
              ii = (fcz(i,j,k,0) < Real(0.0)) ? i - 1 : i + 1;
          else if (apz(i-1,j,k) > Real(0.))
              ii = i-1;
          else
              ii = i+1;

          if (std::abs(fcz(i,j,k,1)) > Real(1.e-8))
              jj = (fcz(i,j,k,1) < Real(0.0)) ? j - 1 : j + 1;
          else if (apz(i,j-1,k) > Real(0.))
              jj = j-1;
          else
              jj = j+1;

          // If any of these cells has zero volume we don't want to use this stencil
          Real test_zero =  cvol(ii,j,k-1) * cvol(i,jj,k-1) *  cvol(ii,jj,k-1) *
                            cvol(ii,j,k  ) * cvol(i,jj,k  ) *  cvol(ii,jj,k  );

          if (test_zero ==  Real(0.))
          {
              int ialt = 2*i - ii;
              int jalt = 2*j - jj;
              Real test_zero_ialt  = cvol(ialt,j,k-1) * cvol(i,jj  ,k-1) *  cvol(ialt,jj  ,k-1) *
                                     cvol(ialt,j,k  ) * cvol(i,jj  ,k  ) *  cvol(ialt,jj  ,k  );

              Real test_zero_jalt  = cvol(ii  ,j,k-1) * cvol(i,jalt,k-1) *  cvol(ii  ,jalt,k-1) *
                                     cvol(ii  ,j,k  ) * cvol(i,jalt,k  ) *  cvol(ii  ,jalt,k  );

              Real test_zero_ijalt = cvol(ialt,j,k-1) * cvol(i,jalt,k-1) *  cvol(ialt,jalt,k-1) *
                                     cvol(ialt,j,k  ) * cvol(i,jalt,k  ) *  cvol(ialt,jalt,k  );

              if (test_zero_ialt > Real(0.))
                 ii = ialt;
              else if (test_zero_jalt > Real(0.))
                 jj = jalt;
              else if (test_zero_ijalt > Real(0.))
              {
                 ii = ialt;
                 jj = jalt;
              }
          }

          Real d0       = Real(0.5) + ccent( i, j,k  ,2);  // distance in z-dir from centroid(i,j,k  ) to face(i,j,k)
          Real d1       = Real(0.5) - ccent( i, j,k-1,2);  // distance in z-dir from centroid(i,j,k-1) to face(i,j,k)

          Real d0_ii    = Real(0.5) + ccent(ii, j,k  ,2); // distance in z-dir from centroid(ii,j,k  ) to face(ii,j,k)
          Real d1_ii    = Real(0.5) - ccent(ii, j,k-1,2); // distance in z-dir from centroid(ii,j,k-1) to face(ii,j,k)

          Real d0_jj    = Real(0.5) + ccent( i,jj,k  ,2); // distance in z-dir from centroid(ii,j,k  ) to face(ii,j,k)
          Real d1_jj    = Real(0.5) - ccent( i,jj,k-1,2); // distance in z-dir from centroid(ii,j,k-1) to face(ii,j,k)

          Real d0_ii_jj = Real(0.5) + ccent(ii,jj,k  ,2); // distance in z-dir from centroid(ii,j,k  ) to face(ii,j,k)
          Real d1_ii_jj = Real(0.5) - ccent(ii,jj,k-1,2); // distance in z-dir from centroid(ii,j,k-1) to face(ii,j,k)

          Real a0       = d1 / (d0 + d1);                           // coefficient of phi( i, j,k  ,n)
          Real a1       = d0 / (d0 + d1);                           // coefficient of phi( i, j,k-1,n)

          Real a0_ii    = d1_ii / (d0_ii + d1_ii);                  // coefficient of phi(ii, j,k  ,n)
          Real a1_ii    = d0_ii / (d0_ii + d1_ii);                  // coefficient of phi(ii, j,k-1,n)

          Real a0_jj    = d1_jj / (d0_jj + d1_jj);                  // coefficient of phi( i,jj,k  ,n)
          Real a1_jj    = d0_jj / (d0_jj + d1_jj);                  // coefficient of phi( i,jj,k-1,n)

          Real a0_ii_jj = d1_ii_jj / (d0_ii_jj + d1_ii_jj);         // coefficient of phi(ii,jj,k  ,n)
          Real a1_ii_jj = d0_ii_jj / (d0_ii_jj + d1_ii_jj);         // coefficient of phi(ii,jj,k-1,n)

          Real phi01       = a0       *   phi( i, j,k,n)  + a1       *   phi( i, j,k-1,n); // interpolation of phi   in z-direction
          Real x01         = a0       * ccent( i, j,k,0)  + a1       * ccent( i, j,k-1,0); // interpolation of x-loc in z-direction
          Real y01         = a0       * ccent( i, j,k,1)  + a1       * ccent( i, j,k-1,1); // interpolation of y-loc in z-direction

          Real phi01_ii    = a0_ii    *   phi(ii, j,k,n)  + a1_ii    *   phi(ii, j,k-1,n); // interpolation of phi   in z-direction
          Real x01_ii      = a0_ii    * ccent(ii, j,k,0)  + a1_ii    * ccent(ii, j,k-1,0); // interpolation of x-loc in z-direction
          Real y01_ii      = a0_ii    * ccent(ii, j,k,1)  + a1_ii    * ccent(ii, j,k-1,1); // interpolation of y-loc in z-direction

          Real phi01_jj    = a0_jj    *   phi( i,jj,k,n)  + a1_jj    *   phi( i,jj,k-1,n); // interpolation of phi   in z-direction
          Real x01_jj      = a0_jj    * ccent( i,jj,k,0)  + a1_jj    * ccent( i,jj,k-1,0); // interpolation of x-loc in z-direction
          Real y01_jj      = a0_jj    * ccent( i,jj,k,1)  + a1_jj    * ccent( i,jj,k-1,1); // interpolation of y-loc in z-direction

          Real phi01_ii_jj = a0_ii_jj *   phi(ii,jj,k,n)  + a1_ii_jj *   phi(ii,jj,k-1,n); // interpolation of phi   in z-direction
          Real x01_ii_jj   = a0_ii_jj * ccent(ii,jj,k,0)  + a1_ii_jj * ccent(ii,jj,k-1,0); // interpolation of x-loc in z-direction
          Real y01_ii_jj   = a0_ii_jj * ccent(ii,jj,k,1)  + a1_ii_jj * ccent(ii,jj,k-1,1); // interpolation of y-loc in z-direction

          // We order these in order to form a set of four points on the x-face ...
          // (x0,y0) : lower left
          // (x1,y1) : lower right
          // (x2,y2) : upper right
          // (x3,y3) : upper left

          Real x,y;

          // This is the location of the face centroid relative to the central node
          // Recall fcz holds (x,y) of the x-face centroid as components (0/1/ )
          if (i < ii)
             x = Real(-0.5) + fcz(i,j,k,0);  // (i,k) is in lower half of stencil so x < 0
          else
             x =  Real(0.5) + fcz(i,j,k,0);  // (i,k) is in upper half of stencil so x > 0

          if (j < jj)
             y = Real(-0.5) + fcz(i,j,k,1);  // (i,k) is in lower half of stencil so z < 0
          else
             y =  Real(0.5) + fcz(i,j,k,1);  // (i,k) is in upper half of stencil so z > 0

          if (i < ii && j < jj) // (i,j) is lower left, (i+1,j+1) is upper right
          {
              Real x_ll = Real(-0.5)+x01;
              Real y_ll = Real(-0.5)+y01;
              Real x_hl =  Real(0.5)+x01_ii;
              Real y_hl = Real(-0.5)+y01_ii;
              Real x_hh =  Real(0.5)+x01_ii_jj;
              Real y_hh =  Real(0.5)+y01_ii_jj;
              Real x_lh = Real(-0.5)+x01_jj;
              Real y_lh =  Real(0.5)+y01_jj;
              phi_z(i,j,k,n) = EB_interp_in_quad(x,y, phi01, phi01_ii, phi01_ii_jj, phi01_jj,
                                                 x_ll, y_ll, x_hl, y_hl, x_hh, y_hh, x_lh, y_lh);
          }
          else if (ii < i && j < jj) // (i-1,j) is lower left, (i,j+1) is upper right
          {
              Real x_ll = Real(-0.5)+x01_ii;
              Real y_ll = Real(-0.5)+y01_ii;
              Real x_hl =  Real(0.5)+x01;
              Real y_hl = Real(-0.5)+y01;
              Real x_hh =  Real(0.5)+x01_jj;
              Real y_hh =  Real(0.5)+y01_jj;
              Real x_lh = Real(-0.5)+x01_ii_jj;
              Real y_lh =  Real(0.5)+y01_ii_jj;
              phi_z(i,j,k,n) = EB_interp_in_quad(x,y,phi01_ii,phi01,phi01_jj, phi01_ii_jj,
                                                 x_ll, y_ll, x_hl, y_hl, x_hh, y_hh, x_lh, y_lh);
          }
          else if (i < ii && jj < j) // (i,j-1) is lower left, (i+1,j) is upper right
          {
              Real x_ll = Real(-0.5)+x01_jj;
              Real y_ll = Real(-0.5)+y01_jj;
              Real x_hl =  Real(0.5)+x01_ii_jj;
              Real y_hl = Real(-0.5)+y01_ii_jj;
              Real x_hh =  Real(0.5)+x01_ii;
              Real y_hh =  Real(0.5)+y01_ii;
              Real x_lh = Real(-0.5)+x01;
              Real y_lh =  Real(0.5)+y01;
              phi_z(i,j,k,n) = EB_interp_in_quad(x,y,phi01_jj,phi01_ii_jj,phi01_ii, phi01,
                                                 x_ll, y_ll, x_hl, y_hl, x_hh, y_hh, x_lh, y_lh);
          }
          else if (ii < i && jj < j) // (i-1,j-1) is lower left, (i,j) is upper right
          {
              Real x_ll = Real(-0.5)+x01_ii_jj;
              Real y_ll = Real(-0.5)+y01_ii_jj;
              Real x_hl =  Real(0.5)+x01_jj;
              Real y_hl = Real(-0.5)+y01_jj;
              Real x_hh =  Real(0.5)+x01;
              Real y_hh =  Real(0.5)+y01;
              Real x_lh = Real(-0.5)+x01_ii;
              Real y_lh =  Real(0.5)+y01_ii;
              phi_z(i,j,k,n) = EB_interp_in_quad(x,y,phi01_ii_jj,phi01_jj,phi01, phi01_ii,
                                                 x_ll, y_ll, x_hl, y_hl, x_hh, y_hh, x_lh, y_lh);
          }
          else
          {
              amrex::Abort("Bad option in interpolation from cell centroid to z-face centroid!");
          }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2face_x (Box const& ubx,
                          Array4<Real const> const& phi,
                          Array4<Real> const& edg_x,
                          int ncomp,
                          const Box& domain,
                          const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (i == domlo.x) && (bc[n].lo(0) == BCType::ext_dir) )
    {
      edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
    }
    else if ( (i == domhi.x+1) && (bc[n].hi(0) == BCType::ext_dir) )
    {
      edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
    }
    else
    {
      edg_x(i,j,k,n) = Real(0.5) * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2face_y (Box const& vbx,
                          Array4<Real const> const& phi,
                          Array4<Real> const& edg_y,
                          int ncomp,
                          const Box& domain,
                          const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (j == domlo.y) && (bc[n].lo(1) == BCType::ext_dir) )
    {
      edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
    }
    else if ( (j == domhi.y+1) && (bc[n].hi(1) == BCType::ext_dir) )
    {
      edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
    }
    else
    {
      edg_y(i,j,k,n) = Real(0.5) * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2face_z (Box const& wbx,
                          Array4<Real const> const& phi,
                          Array4<Real> const& edg_z,
                          int ncomp,
                          const Box& domain,
                          const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (k == domlo.z) && (bc[n].lo(2) == BCType::ext_dir) )
    {
      edg_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
    }
    else if ( (k == domhi.z+1) && (bc[n].hi(2) == BCType::ext_dir) )
    {
      edg_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
    }
    else
    {
      edg_z(i,j,k,n) = Real(0.5) * (phi(i,j  ,k,n) + phi(i,j,k-1,n));
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_add_divergence_from_flow (int i, int j, int k, int n, Array4<Real> const& divu,
                            Array4<Real const> const& vel_eb, Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& vfrc, Array4<Real const> const& bnorm,
                            Array4<Real const> const& barea, GpuArray<Real,3> const& dxinv)
{
    if (flag(i,j,k).isSingleValued())
    {
        Real Ueb_dot_n = vel_eb(i,j,k,0)*bnorm(i,j,k,0)
                       + vel_eb(i,j,k,1)*bnorm(i,j,k,1)
                       + vel_eb(i,j,k,2)*bnorm(i,j,k,2);

        divu(i,j,k,n) += Ueb_dot_n * barea(i,j,k) * dxinv[0] / vfrc(i,j,k);
    }
}


}
#endif
