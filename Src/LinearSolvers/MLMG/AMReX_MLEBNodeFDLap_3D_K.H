#ifndef AMREX_MLEBNODEFDLAP_3D_K_H_
#define AMREX_MLEBNODEFDLAP_3D_K_H_

namespace amrex {

// See the 2D header for the EB position convention.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_scale_rhs (int i, int j, int k, Array4<Real> const& rhs,
                            Array4<int const> const& dmsk, Array4<Real const> const& ebpx,
                            Array4<Real const> const& ebpy, Array4<Real const> const& ebpz) noexcept
{
    if (!dmsk(i,j,k)) {
        Real const s = amrex::min(mlebndfdlap_hm(ebpx(i-1,j,k)), ebpx(i,j,k),
                                  mlebndfdlap_hm(ebpy(i,j-1,k)), ebpy(i,j,k),
                                  mlebndfdlap_hm(ebpz(i,j,k-1)), ebpz(i,j,k));
        rhs(i,j,k) *= s;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_scaled_hinv (Real scale, Real h) noexcept
{
    return (h == Real(0.0)) ? Real(1.0) : scale / h;
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int i, int j, int k, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<int const> const& dmsk,
                        S const& sig, Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const sigxm = sig.x(i-1,j,k);
        Real const sigxp = sig.x(i  ,j,k);
        Real const sigym = sig.y(i,j-1,k);
        Real const sigyp = sig.y(i,j  ,k);
        Real const sigzm = sig.z(i,j,k-1);
        Real const sigzp = sig.z(i,j,k  );
        y(i,j,k) = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
            +      by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
            +      bz * (sigzm*x(i,j,k-1) + sigzp*x(i,j,k+1))
            - (bx*(sigxm+sigxp) + by*(sigym+sigyp) + bz*(sigzm+sigzp)) * x(i,j,k);
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb (int i, int j, int k, Array4<Real> const& x,
                       Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                       S const& sig, Real bx, Real by, Real bz, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const sigxm = sig.x(i-1,j,k);
            Real const sigxp = sig.x(i  ,j,k);
            Real const sigym = sig.y(i,j-1,k);
            Real const sigyp = sig.y(i,j  ,k);
            Real const sigzm = sig.z(i,j,k-1);
            Real const sigzp = sig.z(i,j,k  );
            Real const gamma = -(bx*(sigxm+sigxp) + by*(sigym+sigyp) + bz*(sigzm+sigzp));
            Real const Ax = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
                +           by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
                +           bz * (sigzm*x(i,j,k-1) + sigzp*x(i,j,k+1))
                + gamma * x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

template <typename S, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                Array4<Real const> const& x, Array4<Real const> const& levset,
                                Array4<int const> const& dmsk,
                                Array4<Real const> const& ebpx, Array4<Real const> const& ebpy,
                                Array4<Real const> const& ebpz,
                                S const& sig, F const& xeb, Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = ebpx(i  ,j  ,k  );
        Real const hmx = mlebndfdlap_hm(ebpx(i-1,j  ,k  ));
        Real const hpy = ebpy(i  ,j  ,k  );
        Real const hmy = mlebndfdlap_hm(ebpy(i  ,j-1,k  ));
        Real const hpz = ebpz(i  ,j  ,k  );
        Real const hmz = mlebndfdlap_hm(ebpz(i  ,j  ,k-1));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

        Real out = Real(0.0);
        Real fac = bx * Real(2.0) / (hpx+hmx);
        Real sigma = sig.x(i,j,k);
        if (levset(i+1,j,k) < Real(0.0) && hpx == Real(1.0)) {
            out += sigma * fac * scale * (x(i+1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx) * (xeb(i+1,j,k) - xc);
        }
        sigma = sig.x(i-1,j,k);
        if (levset(i-1,j,k) < Real(0.0) && hmx == Real(1.0)) {
            out += sigma * fac * scale * (x(i-1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx) * (xeb(i-1,j,k) - xc);
        }

        fac = by * Real(2.0) / (hpy+hmy);
        sigma = sig.y(i,j,k);
        if (levset(i,j+1,k) < Real(0.0) && hpy == Real(1.0)) {
            out += sigma * fac * scale * (x(i,j+1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy) * (xeb(i,j+1,k) - xc);
        }
        sigma = sig.y(i,j-1,k);
        if (levset(i,j-1,k) < Real(0.0) && hmy == Real(1.0)) {
            out += sigma * fac * scale * (x(i,j-1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy) * (xeb(i,j-1,k) - xc);
        }

        fac = bz * Real(2.0) / (hpz+hmz);
        sigma = sig.z(i,j,k);
        if (levset(i,j,k+1) < Real(0.0) && hpz == Real(1.0)) {
            out += sigma * fac * scale * (x(i,j,k+1) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpz) * (xeb(i,j,k+1) - xc);
        }
        sigma = sig.z(i,j,k-1);
        if (levset(i,j,k-1) < Real(0.0) && hmz == Real(1.0)) {
            out += sigma * fac * scale * (x(i,j,k-1) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmz) * (xeb(i,j,k-1) - xc);
        }

        y(i,j,k) = out;
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ebpx, Array4<Real const> const& ebpy,
                           Array4<Real const> const& ebpz,
                           S const& sig, Real xeb, Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ebpx, ebpy, ebpz, sig,
                              [=] (int, int, int) -> Real { return xeb; },
                              bx, by, bz);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ebpx, Array4<Real const> const& ebpy,
                           Array4<Real const> const& ebpz,
                           S const& sig, Array4<Real const> const& xeb,
                           Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ebpx, ebpy, ebpz, sig,
                              [=] (int i1, int i2, int i3) -> Real {
                                  return xeb(i1,i2,i3); },
                              bx, by, bz);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<Real const> const& levset,
                          Array4<int const> const& dmsk,
                          Array4<Real const> const& ebpx, Array4<Real const> const& ebpy,
                          Array4<Real const> const& ebpz,
                          S const& sig, Real bx, Real by, Real bz, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = ebpx(i  ,j  ,k  );
            Real const hmx = mlebndfdlap_hm(ebpx(i-1,j  ,k  ));
            Real const hpy = ebpy(i  ,j  ,k  );
            Real const hmy = mlebndfdlap_hm(ebpy(i  ,j-1,k  ));
            Real const hpz = ebpz(i  ,j  ,k  );
            Real const hmz = mlebndfdlap_hm(ebpz(i  ,j  ,k-1));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

            Real gamma = Real(0.0);
            Real rho = Real(0.0);
            Real fac = bx * Real(2.0) / (hpx+hmx);
            Real sigma = sig.x(i,j,k);
            if (levset(i+1,j,k) < Real(0.0) && hpx == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i+1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx);
            }
            sigma = sig.x(i-1,j,k);
            if (levset(i-1,j,k) < Real(0.0) && hmx == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i-1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx);
            }

            fac = by * Real(2.0) / (hpy+hmy);
            sigma = sig.y(i,j,k);
            if (levset(i,j+1,k) < Real(0.0) && hpy == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j+1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy);
            }
            sigma = sig.y(i,j-1,k);
            if (levset(i,j-1,k) < Real(0.0) && hmy == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j-1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy);
            }

            fac = bz * Real(2.0) / (hpz+hmz);
            sigma = sig.z(i,j,k);
            if (levset(i,j,k+1) < Real(0.0) && hpz == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j,k+1);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpz);
            }
            sigma = sig.z(i,j,k-1);
            if (levset(i,j,k-1) < Real(0.0) && hmz == Real(1.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j,k-1);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmz);
            }

            Real const Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

}

#endif
