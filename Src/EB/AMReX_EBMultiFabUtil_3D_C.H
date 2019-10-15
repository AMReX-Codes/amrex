#ifndef AMREX_EB_MULTIFAB_UTIL_3D_C_H_
#define AMREX_EB_MULTIFAB_UTIL_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real const * AMREX_RESTRICT v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
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
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = fv(ii,jj,kk)*vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
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
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
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
        Real c = 0.0;
        Real a = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
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
        Real c = 0.0;
        Real a = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
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
        Real c = 0.0;
        Real a = 0.0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
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
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = ba(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = 0.0;
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
        divu(i,j,k,n) = 0.0;
    }
    else if (flag(i,j,k).isRegular())
    {
        divu(i,j,k,n) = dxinv[0] * (u(i+1,j,k,n)-u(i,j,k,n))
            +           dxinv[1] * (v(i,j+1,k,n)-v(i,j,k,n))
            +           dxinv[2] * (w(i,j,k+1,n)-w(i,j,k,n));
    }
    else if (already_on_centroids)
    {
        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) * (
                        dxinv[0] * (apx(i+1,j,k)*u(i+1,j,k,n)-apx(i,j,k)*u(i,j,k,n))
            +           dxinv[1] * (apy(i,j+1,k)*v(i,j+1,k,n)-apy(i,j,k)*v(i,j,k,n))
            +           dxinv[2] * (apz(i,j,k+1)*w(i,j,k+1,n)-apz(i,j,k)*w(i,j,k,n)) );
    }
    else
    {
        Real fxm = u(i,j,k,n);
        if (apx(i,j,k) != 0.0 and apx(i,j,k) != 1.0) {
            int jj = j + static_cast<int>(std::copysign(1.0, fcx(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(1.0, fcx(i,j,k,1)));
            Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : 0.0;
            Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : 0.0;
            fxm = (1.0-fracy)*(1.0-fracz)*fxm
                +      fracy *(1.0-fracz)*u(i,jj,k ,n)
                +      fracz *(1.0-fracy)*u(i,j ,kk,n)
                +      fracy *     fracz *u(i,jj,kk,n);
        }

        Real fxp = u(i+1,j,k,n);
        if (apx(i+1,j,k) != 0.0 and apx(i+1,j,k) != 1.0) {
            int jj = j + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k,1)));
            Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k,0)) : 0.0;
            Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? std::abs(fcx(i+1,j,k,1)) : 0.0;
            fxp = (1.0-fracy)*(1.0-fracz)*fxp
                +      fracy *(1.0-fracz)*u(i+1,jj,k ,n)
                +      fracz *(1.0-fracy)*u(i+1,j ,kk,n)
                +      fracy *     fracz *u(i+1,jj,kk,n);
        }

        Real fym = v(i,j,k,n);
        if (apy(i,j,k) != 0.0 and apy(i,j,k) != 1.0) {
            int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j,k,1)));
            Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : 0.0;
            Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : 0.0;
            fym = (1.0-fracx)*(1.0-fracz)*fym
                +      fracx *(1.0-fracz)*v(ii,j,k ,n)
                +      fracz *(1.0-fracx)*v(i ,j,kk,n)
                +      fracx *     fracz *v(ii,j,kk,n);
        }

        Real fyp = v(i,j+1,k,n);
        if (apy(i,j+1,k) != 0.0 and apy(i,j+1,k) != 1.0) {
            int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k,0)));
            int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k,0)) : 0.0;
            Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? std::abs(fcy(i,j+1,k,1)) : 0.0;
            fyp = (1.0-fracx)*(1.0-fracz)*fyp
                +      fracx *(1.0-fracz)*v(ii,j+1,k ,n)
                +      fracz *(1.0-fracx)*v(i ,j+1,kk,n)
                +      fracx *     fracz *v(ii,j+1,kk,n);
        }

        Real fzm = w(i,j,k,n);
        if (apz(i,j,k) != 0.0 and apz(i,j,k) != 1.0) {
            int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k,0)));
            int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k,1)));
            Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : 0.0;
            Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : 0.0;
            fzm = (1.0-fracx)*(1.0-fracy)*fzm
                +      fracx *(1.0-fracy)*w(ii,j ,k,n)
                +      fracy *(1.0-fracx)*w(i ,jj,k,n)
                +      fracx *     fracy *w(ii,jj,k,n);
        }

        Real fzp = w(i,j,k+1,n);
        if (apz(i,j,k+1) != 0.0 and apz(i,j,k+1) != 1.0) {
            int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k+1,0)));
            int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k+1,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? std::abs(fcz(i,j,k+1,0)) : 0.0;
            Real fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? std::abs(fcz(i,j,k+1,1)) : 0.0;
            fzp = (1.0-fracx)*(1.0-fracy)*fzp
                +      fracx *(1.0-fracy)*w(ii,j ,k+1,n)
                +      fracy *(1.0-fracx)*w(i ,jj,k+1,n)
                +      fracx *     fracy *w(ii,jj,k+1,n);
        }

        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) *
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
        cc(i,j,k,n+0) = 0.0;
        cc(i,j,k,n+1) = 0.0;
        cc(i,j,k,n+2) = 0.0;
    } else {
        if (ax(i,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i+1,j,k);
        } else if (ax(i+1,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i,j,k);
        } else {
            cc(i,j,k,n+0) = 0.5 * (fx(i,j,k) + fx(i+1,j,k));
        }

        if (ay(i,j,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j+1,k);
        } else if (ay(i,j+1,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j,k);
        } else {
            cc(i,j,k,n+1) = 0.5 * (fy(i,j,k) + fy(i,j+1,k));
        }

        if (az(i,j,k) == 0.0) {
            cc(i,j,k,n+2) = fz(i,j,k+1);
        } else if (az(i,j,k+1) == 0.0) {
            cc(i,j,k,n+2) = fz(i,j,k);
        } else {
            cc(i,j,k,n+2) = 0.5 * (fz(i,j,k) + fz(i,j,k+1));
        }
    }
}

}
#endif
