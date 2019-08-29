#ifndef AMREX_EB_FLUXREGISTER_3D_C_H_
#define AMREX_EB_FLUXREGISTER_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_crseadd_va(int i, int j, int k, Array4<Real> const& d,
                            Array4<int const> const& flag, Array4<Real const> const& fx,
                            Array4<Real const> const& fy, Array4<Real const> const& fz,
                            Array4<Real const> const& vfrac, Array4<Real const> const& ax,
                            Array4<Real const> const& ay, Array4<Real const> const& az,
                            Real dtdx, Real dtdy, Real dtdz)
{
    if (flag(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell
        and vfrac(i,j,k) > 1.e-14)
    {
        int ncomp = d.nComp();
        Real volinv = 1.0/vfrac(i,j,k);

        if (flag(i-1,j,k) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdx*ax(i,j,k)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) -= tmp*fx(i,j,k,n);
            }
        } else if (flag(i+1,j,k) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdx*ax(i+1,j,k)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) += tmp*fx(i+1,j,k,n);
            }
        }

        if (flag(i,j-1,k) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdy*ay(i,j,k)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) -= tmp*fy(i,j,k,n);
            }
        } else if (flag(i,j+1,k) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdy*ay(i,j+1,k)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) += tmp*fy(i,j+1,k,n);
            }
        }

        if (flag(i,j,k-1) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdz*az(i,j,k)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) -= tmp*fz(i,j,k,n);
            }
        } else if (flag(i,j,k+1) == amrex_yafluxreg_fine_cell) {
            Real tmp = dtdz*az(i,j,k+1)*volinv;
            for (int n = 0; n < ncomp; ++n) {
                d(i,j,k,n) += tmp*fz(i,j,k+1,n);
            }
        }
    }
}

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real eb_flux_reg_cvol (int i, int j, int k, Array4<Real const> const& vfrac,
                           Dim3 const& ratio, Real small) noexcept
    {
        Real cvol = 0.0;
        for         (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
            for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
                for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
                    cvol += vfrac(ii,jj,kk);
                }
            }
        }
        return (cvol > small) ? 1.0/cvol : 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_xlo (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int ii = (i+1)*ratio.x;
    Real fa = 0.0;
    for     (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= -fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_xhi (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int ii = i*ratio.x;
    Real fa = 0.0;
    for     (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_ylo (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int jj = (j+1)*ratio.y;
    Real fa = 0.0;
    for     (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= -fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_yhi (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int jj = j*ratio.y;
    Real fa = 0.0;
    for     (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_zlo (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int kk = (k+1)*ratio.z;
    Real fa = 0.0;
    for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= -fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_zhi (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    int kk = k*ratio.z;
    Real fa = 0.0;
    for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            if (f.contains(ii,jj,kk)) {
                fa += f(ii,jj,kk,n) * a(ii,jj,kk);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,1.e-14);
    fa *= fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_dm (int i, int j, int k, int n, Box const& dmbx, Array4<Real> const& d,
                             Array4<Real const> const& dm, Array4<Real const> const& vfrac,
                             Dim3 const& ratio, Real threshold)
{
    Real dmtot = 0.0;
    for         (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
                if (dmbx.contains(IntVect(ii,jj,kk))) {
                    dmtot += dm(ii,jj,kk,n);
                }
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,k,vfrac,ratio,threshold);
    dmtot *= cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), dmtot);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_rereflux_from_crse (int i, int j, int k, int n, Box const& bx, Array4<Real> const& d,
                            Array4<Real const> const& s, Array4<int const> const& amrflg,
                            Array4<EBCellFlag const> const& ebflg, Array4<Real const> const& vfrac)
{
    if (amrflg(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell)
    {
        auto flag = ebflg(i,j,k);
        if (flag.isRegular())
        {
            if (bx.contains(IntVect(i,j,k)))
            {
                HostDevice::Atomic::Add(d.ptr(i,j,k,n), s(i,j,k,n));
            }
        }
        else if (flag.isSingleValued())
        {
            Real dm = s(i,j,k,n) * vfrac(i,j,k);
            if (bx.contains(IntVect(i,j,k))) {
                HostDevice::Atomic::Add(d.ptr(i,j,k,n), dm);
            }

            Real wtot = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and flag.isConnected(ii,jj,kk)) {
                    wtot += vfrac(i+ii,j+jj,k+kk);
                }
            }}}

            Real drho = dm * ((1.0-vfrac(i,j,k))/wtot);
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and flag.isConnected(ii,jj,kk)) {
                    if (bx.contains(IntVect(i+ii,j+jj,k+kk))) {
                        HostDevice::Atomic::Add(d.ptr(i+ii,j+jj,k+kk,n), drho);
                    }
                }
            }}}
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_rereflux_to_fine (int i, int j, int k, int n, Array4<Real> const& d,
                          Array4<Real const> const& s, Array4<int const> const& msk, Dim3 ratio)
{
    int ic = amrex::coarsen(i,ratio.x);
    int jc = amrex::coarsen(j,ratio.y);
    int kc = amrex::coarsen(k,ratio.z);
    if (msk(ic,jc,kc) == 1) {
        d(i,j,k,n) += s(ic,jc,kc,n);
    }
}

}

#endif
