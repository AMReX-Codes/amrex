#ifndef AMREX_EB_FLUXREGISTER_2D_C_H_
#define AMREX_EB_FLUXREGISTER_2D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_crseadd_va(int i, int j, int k, int n, Array4<Real> const& d,
                            Array4<int const> const& flag, Array4<Real const> const& fx,
                            Array4<Real const> const& fy, Array4<Real const> const& vfrac,
                            Array4<Real const> const& ax, Array4<Real const> const& ay,
                            Real dtdx, Real dtdy)
{
    if (flag(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell
        and vfrac(i,j,k) > 1.e-14)
    {
        Real volinv = 1.0/vfrac(i,j,k);

        if (flag(i-1,j,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) -= dtdx*fx(i,j,k,n)*(ax(i,j,k)*volinv);
        } else if (flag(i+1,j,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) +=  dtdx*fx(i+1,j,k,n)*(ax(i+1,j,k)*volinv);
        }

        if (flag(i,j-1,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) -= dtdy*fy(i,j,k,n)*(ay(i,j,k)*volinv);
        } else if (flag(i,j+1,k) == amrex_yafluxreg_fine_cell) {
            d(i,j,k,n) += dtdy*fy(i,j+1,k,n)*(ay(i,j+1,k)*volinv);
        }
    }
}

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real eb_flux_reg_cvol (int i, int j, Array4<Real const> const& vfrac,
                           Dim3 const& ratio, Real threshold) noexcept
    {
        Real cvol = 0.0;
        constexpr int kk = 0;
        for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
                cvol += vfrac(ii,jj,kk);
            }
        }
        return (cvol > threshold) ? 1.0/cvol : 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_xlo (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    constexpr int kk = 0;
    int ii = (i+1)*ratio.x;
    Real fa = 0.0;
    for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        if (f.contains(ii,jj,kk)) {
            fa += f(ii,jj,kk,n) * a(ii,jj,kk);
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,vfrac,ratio,1.e-14);
    fa *= -fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_xhi (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    constexpr int kk = 0;
    int ii = i*ratio.x;
    Real fa = 0.0;
    for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        if (f.contains(ii,jj,kk)) {
            fa += f(ii,jj,kk,n) * a(ii,jj,kk);
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,vfrac,ratio,1.e-14);
    fa *= fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_ylo (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    constexpr int kk = 0;
    int jj = (j+1)*ratio.y;
    Real fa = 0.0;
    for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
        if (f.contains(ii,jj,kk)) {
            fa += f(ii,jj,kk,n) * a(ii,jj,kk);
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,vfrac,ratio,1.e-14);
    fa *= -fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_va_yhi (int i, int j, int k, int n, Array4<Real> const& d,
                                 Array4<Real const> const& f, Array4<Real const> const& vfrac,
                                 Array4<Real const> const& a, Real fac, Dim3 const& ratio)
{
    constexpr int kk = 0;
    int jj = j*ratio.y;
    Real fa = 0.0;
    for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
        if (f.contains(ii,jj,kk)) {
            fa += f(ii,jj,kk,n) * a(ii,jj,kk);
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,vfrac,ratio,1.e-14);
    fa *= fac*cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), fa);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_flux_reg_fineadd_dm (int i, int j, int k, int n, Box const& dmbx, Array4<Real> const& d,
                             Array4<Real const> const& dm, Array4<Real const> const& vfrac,
                             Dim3 const& ratio, Real threshold)
{
    Real dmtot = 0.0;
    constexpr int kk = 0;
    for     (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            if (dmbx.contains(IntVect(ii,jj))) {
                dmtot += dm(ii,jj,kk,n);
            }
        }
    }
    Real cvol = eb_flux_reg_cvol(i,j,vfrac,ratio,threshold);
    dmtot *= cvol;
    HostDevice::Atomic::Add(d.ptr(i,j,k,n), dmtot);
}

}

#endif
