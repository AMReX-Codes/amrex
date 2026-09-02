#ifndef AMREX_MLEBNODEFDLAP_3D_K_H_
#define AMREX_MLEBNODEFDLAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_scale_rhs (int i, int j, int k, Array4<Real> const& rhs,
                            Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                            Array4<Real const> const& ecy, Array4<Real const> const& ecz) noexcept
{
    if (!dmsk(i,j,k)) {
        Real hmx = (ecx(i-1,j  ,k  ) == Real(1.)) ? Real(1.) : Real(1.) - Real(2.)*ecx(i-1,j  ,k  );
        Real hpx = (ecx(i  ,j  ,k  ) == Real(1.)) ? Real(1.) : Real(1.) + Real(2.)*ecx(i  ,j  ,k  );
        Real hmy = (ecy(i  ,j-1,k  ) == Real(1.)) ? Real(1.) : Real(1.) - Real(2.)*ecy(i  ,j-1,k  );
        Real hpy = (ecy(i  ,j  ,k  ) == Real(1.)) ? Real(1.) : Real(1.) + Real(2.)*ecy(i  ,j  ,k  );
        Real hmz = (ecz(i  ,j  ,k-1) == Real(1.)) ? Real(1.) : Real(1.) - Real(2.)*ecz(i  ,j  ,k-1);
        Real hpz = (ecz(i  ,j  ,k  ) == Real(1.)) ? Real(1.) : Real(1.) + Real(2.)*ecz(i  ,j  ,k  );
        Real const s = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);
        rhs(i,j,k) *= s;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_scaled_hinv (Real scale, Real h) noexcept
{
    return (h == Real(0.0)) ? Real(1.0) : scale / h;
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                Array4<Real const> const& x, Array4<Real const> const& levset,
                                Array4<int const> const& dmsk,
                                Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                Array4<Real const> const& ecz, F const& xeb,
                                Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
        Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
        Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real const hpz = (ecz(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.0)*ecz(i,j,k));
        Real const hmz = (ecz(i,j,k-1) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecz(i,j,k-1));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

        Real out = Real(0.0);
        Real fac = bx * Real(2.0) / (hpx+hmx);
        if (levset(i+1,j,k) < Real(0.0)) { // regular
            out += fac * scale * (x(i+1,j,k) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hpx) * (xeb(i+1,j,k) - xc);
        }
        if (levset(i-1,j,k) < Real(0.0)) {
            out += fac * scale * (x(i-1,j,k) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hmx) * (xeb(i-1,j,k) - xc);
        }

        fac = by * Real(2.0) / (hpy+hmy);
        if (levset(i,j+1,k) < Real(0.0)) {
            out += fac * scale * (x(i,j+1,k) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hpy) * (xeb(i,j+1,k) - xc);
        }
        if (levset(i,j-1,k) < Real(0.0)) {
            out += fac * scale * (x(i,j-1,k) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hmy) * (xeb(i,j-1,k) - xc);
        }

        fac = bz * Real(2.0) / (hpz+hmz);
        if (levset(i,j,k+1) < Real(0.0)) {
            out += fac * scale * (x(i,j,k+1) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hpz) * (xeb(i,j,k+1) - xc);
        }
        if (levset(i,j,k-1) < Real(0.0)) {
            out += fac * scale * (x(i,j,k-1) - xc);
        } else {
            out += fac * mlebndfdlap_scaled_hinv(scale, hmz) * (xeb(i,j,k-1) - xc);
        }

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Array4<Real const> const& ecz, Real xeb,
                           Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, ecz,
                              [=] (int, int, int) -> Real { return xeb; },
                              bx, by, bz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Array4<Real const> const& ecz, Array4<Real const> const& xeb,
                           Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, ecz,
                              [=] (int i1, int i2, int i3) -> Real {
                                  return xeb(i1,i2,i3); },
                              bx, by, bz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int i, int j, int k, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<int const> const& dmsk,
                        Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        y(i,j,k) = bx * (x(i-1,j,k) + x(i+1,j,k))
            +      by * (x(i,j-1,k) + x(i,j+1,k))
            +      bz * (x(i,j,k-1) + x(i,j,k+1))
            - (Real(2.0)*(bx+by+bz)) * x(i,j,k);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<Real const> const& levset,
                          Array4<int const> const& dmsk,
                          Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                          Array4<Real const> const& ecz, Real bx, Real by, Real bz,
                          int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real const hpz = (ecz(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.0)*ecz(i,j,k));
            Real const hmz = (ecz(i,j,k-1) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecz(i,j,k-1));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

            Real gamma = Real(0.0);
            Real rho = Real(0.0);
            Real fac = bx * Real(2.0) / (hpx+hmx);
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                gamma -= fac * scale;
                rho += fac * scale * x(i+1,j,k);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hpx);
            }
            if (levset(i-1,j,k) < Real(0.0)) {
                gamma -= fac * scale;
                rho += fac * scale * x(i-1,j,k);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hmx);
            }

            fac = by * Real(2.0) / (hpy+hmy);
            if (levset(i,j+1,k) < Real(0.0)) {
                gamma -= fac * scale;
                rho += fac * scale * x(i,j+1,k);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hpy);
            }
            if (levset(i,j-1,k) < Real(0.0)) {
                gamma -= fac * scale;
                rho += fac * scale * x(i,j-1,k);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hmy);
            }

            fac = bz * Real(2.0) / (hpz+hmz);
            if (levset(i,j,k+1) < Real(0.0)) {
                gamma -= fac * scale;
                rho += fac * scale * x(i,j,k+1);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hpz);
            }
            if (levset(i,j,k-1) < Real(0.0)) {
                gamma -= fac * scale;
                rho += fac * scale * x(i,j,k-1);
            } else {
                gamma -= fac * mlebndfdlap_scaled_hinv(scale, hmz);
            }

            Real Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb (int i, int j, int k, Array4<Real> const& x,
                       Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                       Real bx, Real by, Real bz, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real gamma = Real(-2.0)*(bx+by+bz);
            Real Ax = bx * (x(i-1,j,k) + x(i+1,j,k))
                +     by * (x(i,j-1,k) + x(i,j+1,k))
                +     bz * (x(i,j,k-1) + x(i,j,k+1))
                + gamma * x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx (int i, int j, int k, Array4<Real> const& y,
                            Array4<Real const> const& x,
                            Array4<int const> const& dmsk,
                            Array4<Real const> const& sig,
                            Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real sigxm = Real(0.25)*(sig(i-1,j-1,k-1) +
                                 sig(i-1,j  ,k-1) +
                                 sig(i-1,j-1,k  ) +
                                 sig(i-1,j  ,k  ));
        Real sigxp = Real(0.25)*(sig(i  ,j-1,k-1) +
                                 sig(i  ,j  ,k-1) +
                                 sig(i  ,j-1,k  ) +
                                 sig(i  ,j  ,k  ));
        Real sigym = Real(0.25)*(sig(i-1,j-1,k-1) +
                                 sig(i  ,j-1,k-1) +
                                 sig(i-1,j-1,k  ) +
                                 sig(i  ,j-1,k  ));
        Real sigyp = Real(0.25)*(sig(i-1,j  ,k-1) +
                                 sig(i  ,j  ,k-1) +
                                 sig(i-1,j  ,k  ) +
                                 sig(i  ,j  ,k  ));
        Real sigzm = Real(0.25)*(sig(i-1,j-1,k-1) +
                                 sig(i  ,j-1,k-1) +
                                 sig(i-1,j  ,k-1) +
                                 sig(i  ,j  ,k-1));
        Real sigzp = Real(0.25)*(sig(i-1,j-1,k  ) +
                                 sig(i  ,j-1,k  ) +
                                 sig(i-1,j  ,k  ) +
                                 sig(i  ,j  ,k  ));
        y(i,j,k) = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
            +      by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
            +      bz * (sigzm*x(i,j,k-1) + sigzp*x(i,j,k+1))
            - (bx*(sigxm+sigxp) + by*(sigym+sigyp) + bz*(sigzm+sigzp)) * x(i,j,k);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb (int i, int j, int k, Array4<Real> const& x,
                           Array4<Real const> const& rhs,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& sig,
                           Real bx, Real by, Real bz, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real sigxm = Real(0.25)*(sig(i-1,j-1,k-1) +
                                     sig(i-1,j  ,k-1) +
                                     sig(i-1,j-1,k  ) +
                                     sig(i-1,j  ,k  ));
            Real sigxp = Real(0.25)*(sig(i  ,j-1,k-1) +
                                     sig(i  ,j  ,k-1) +
                                     sig(i  ,j-1,k  ) +
                                     sig(i  ,j  ,k  ));
            Real sigym = Real(0.25)*(sig(i-1,j-1,k-1) +
                                     sig(i  ,j-1,k-1) +
                                     sig(i-1,j-1,k  ) +
                                     sig(i  ,j-1,k  ));
            Real sigyp = Real(0.25)*(sig(i-1,j  ,k-1) +
                                     sig(i  ,j  ,k-1) +
                                     sig(i-1,j  ,k  ) +
                                     sig(i  ,j  ,k  ));
            Real sigzm = Real(0.25)*(sig(i-1,j-1,k-1) +
                                     sig(i  ,j-1,k-1) +
                                     sig(i-1,j  ,k-1) +
                                     sig(i  ,j  ,k-1));
            Real sigzp = Real(0.25)*(sig(i-1,j-1,k  ) +
                                     sig(i  ,j-1,k  ) +
                                     sig(i-1,j  ,k  ) +
                                     sig(i  ,j  ,k  ));
            Real gamma = -(bx*(sigxm+sigxp) + by*(sigym+sigyp) + bz*(sigzm+sigzp));
            Real Ax = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
                +     by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
                +     bz * (sigzm*x(i,j,k-1) + sigzp*x(i,j,k+1))
                + gamma * x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                    Array4<Real const> const& x,
                                    Array4<Real const> const& levset,
                                    Array4<int const> const& dmsk,
                                    Array4<Real const> const& ecx,
                                    Array4<Real const> const& ecy,
                                    Array4<Real const> const& ecz,
                                    Array4<Real const> const& sig,
                                    Array4<Real const> const& vfrc,
                                    F const& xeb,
                                    Real bx, Real by, Real bz) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
        Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
        Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real const hpz = (ecz(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.0)*ecz(i,j,k));
        Real const hmz = (ecz(i,j,k-1) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecz(i,j,k-1));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

        Real out = Real(0.0);
        Real sigma = (sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                      sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1) +
                      sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ) +
                      sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
            /        (                   vfrc(i  ,j-1,k-1) +
                                         vfrc(i  ,j  ,k-1) +
                                         vfrc(i  ,j-1,k  ) +
                                         vfrc(i  ,j  ,k  ));
        Real fac = bx * Real(2.0) / (hpx+hmx);
        if (levset(i+1,j,k) < Real(0.0)) { // regular
            out += sigma * fac * scale * (x(i+1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx) * (xeb(i+1,j,k) - xc);
        }

        sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                 sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                 sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                 sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ))
            /   (                   vfrc(i-1,j-1,k-1) +
                                    vfrc(i-1,j  ,k-1) +
                                    vfrc(i-1,j-1,k  ) +
                                    vfrc(i-1,j  ,k  ));
        if (levset(i-1,j,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i-1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx) * (xeb(i-1,j,k) - xc);
        }

        sigma = (sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                 sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1) +
                 sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ) +
                 sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
            /   (                   vfrc(i-1,j  ,k-1) +
                                    vfrc(i  ,j  ,k-1) +
                                    vfrc(i-1,j  ,k  ) +
                                    vfrc(i  ,j  ,k  ));
        fac = by * Real(2.0) / (hpy+hmy);
        if (levset(i,j+1,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j+1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy) * (xeb(i,j+1,k) - xc);
        }

        sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                 sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                 sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                 sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ))
            /   (                   vfrc(i-1,j-1,k-1) +
                                    vfrc(i  ,j-1,k-1) +
                                    vfrc(i-1,j-1,k  ) +
                                    vfrc(i  ,j-1,k  ));
        if (levset(i,j-1,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j-1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy) * (xeb(i,j-1,k) - xc);
        }

        sigma = (sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                 sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ) +
                 sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ) +
                 sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
            /   (                   vfrc(i-1,j-1,k  ) +
                                    vfrc(i  ,j-1,k  ) +
                                    vfrc(i-1,j  ,k  ) +
                                    vfrc(i  ,j  ,k  ));
        fac = bz * Real(2.0) / (hpz+hmz);
        if (levset(i,j,k+1) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j,k+1) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpz) * (xeb(i,j,k+1) - xc);
        }

        sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                 sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                 sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                 sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1))
            /   (                   vfrc(i-1,j-1,k-1) +
                                    vfrc(i  ,j-1,k-1) +
                                    vfrc(i-1,j  ,k-1) +
                                    vfrc(i  ,j  ,k-1));
        if (levset(i,j,k-1) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j,k-1) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmz) * (xeb(i,j,k-1) - xc);
        }

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                               Array4<Real const> const& x,
                               Array4<Real const> const& levset,
                               Array4<int const> const& dmsk,
                               Array4<Real const> const& ecx,
                               Array4<Real const> const& ecy,
                               Array4<Real const> const& ecz,
                               Array4<Real const> const& sig,
                               Array4<Real const> const& vfrc,
                               Real xeb, Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_sig_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, ecz, sig, vfrc,
                                  [=] (int, int, int) -> Real { return xeb; },
                                  bx, by, bz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                               Array4<Real const> const& x,
                               Array4<Real const> const& levset,
                               Array4<int const> const& dmsk,
                               Array4<Real const> const& ecx,
                               Array4<Real const> const& ecy,
                               Array4<Real const> const& ecz,
                               Array4<Real const> const& sig,
                               Array4<Real const> const& vfrc,
                               Array4<Real const> const& xeb,
                               Real bx, Real by, Real bz) noexcept
{
    mlebndfdlap_sig_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, ecz, sig, vfrc,
                                  [=] (int i1, int i2, int i3) -> Real {
                                      return xeb(i1,i2,i3); },
                                  bx, by, bz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                              Array4<Real const> const& rhs,
                              Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx,
                              Array4<Real const> const& ecy,
                              Array4<Real const> const& ecz,
                              Array4<Real const> const& sig,
                              Array4<Real const> const& vfrc,
                              Real bx, Real by, Real bz, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real const hpz = (ecz(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.0)*ecz(i,j,k));
            Real const hmz = (ecz(i,j,k-1) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecz(i,j,k-1));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);

            Real gamma = Real(0.0);
            Real rho = Real(0.0);
            Real sigma = (sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                          sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1) +
                          sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ) +
                          sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
                /        (                   vfrc(i  ,j-1,k-1) +
                                             vfrc(i  ,j  ,k-1) +
                                             vfrc(i  ,j-1,k  ) +
                                             vfrc(i  ,j  ,k  ));
            Real fac = bx * Real(2.0) / (hpx+hmx);
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i+1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx);
            }

            sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                     sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                     sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                     sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ))
                /   (                   vfrc(i-1,j-1,k-1) +
                                        vfrc(i-1,j  ,k-1) +
                                        vfrc(i-1,j-1,k  ) +
                                        vfrc(i-1,j  ,k  ));
            if (levset(i-1,j,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i-1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx);
            }

            sigma = (sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                     sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1) +
                     sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ) +
                     sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
                /   (                   vfrc(i-1,j  ,k-1) +
                                        vfrc(i  ,j  ,k-1) +
                                        vfrc(i-1,j  ,k  ) +
                                        vfrc(i  ,j  ,k  ));
            fac = by * Real(2.0) / (hpy+hmy);
            if (levset(i,j+1,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j+1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy);
            }

            sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                     sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                     sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                     sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ))
                /   (                   vfrc(i-1,j-1,k-1) +
                                        vfrc(i  ,j-1,k-1) +
                                        vfrc(i-1,j-1,k  ) +
                                        vfrc(i  ,j-1,k  ));
            if (levset(i,j-1,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j-1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy);
            }

            sigma = (sig(i-1,j-1,k  ) * vfrc(i-1,j-1,k  ) +
                     sig(i  ,j-1,k  ) * vfrc(i  ,j-1,k  ) +
                     sig(i-1,j  ,k  ) * vfrc(i-1,j  ,k  ) +
                     sig(i  ,j  ,k  ) * vfrc(i  ,j  ,k  ))
                /   (                   vfrc(i-1,j-1,k  ) +
                                        vfrc(i  ,j-1,k  ) +
                                        vfrc(i-1,j  ,k  ) +
                                        vfrc(i  ,j  ,k  ));
            fac = bz * Real(2.0) / (hpz+hmz);
            if (levset(i,j,k+1) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j,k+1);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpz);
            }

            sigma = (sig(i-1,j-1,k-1) * vfrc(i-1,j-1,k-1) +
                     sig(i  ,j-1,k-1) * vfrc(i  ,j-1,k-1) +
                     sig(i-1,j  ,k-1) * vfrc(i-1,j  ,k-1) +
                     sig(i  ,j  ,k-1) * vfrc(i  ,j  ,k-1))
                /   (                   vfrc(i-1,j-1,k-1) +
                                        vfrc(i  ,j-1,k-1) +
                                        vfrc(i-1,j  ,k-1) +
                                        vfrc(i  ,j  ,k-1));
            if (levset(i,j,k-1) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j,k-1);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmz);
            }

            Real Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

}

#endif
