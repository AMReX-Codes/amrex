#ifndef AMREX_MLEBNODEFDLAP_2D_K_H_
#define AMREX_MLEBNODEFDLAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_scale_rhs (int i, int j, int, Array4<Real> const& rhs,
                            Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                            Array4<Real const> const& ecy) noexcept
{
    if (!dmsk(i,j,0)) {
        Real hmx = (ecx(i-1,j  ,0) == Real(1.)) ? Real(1.) : Real(1.) - Real(2.)*ecx(i-1,j  ,0);
        Real hpx = (ecx(i  ,j  ,0) == Real(1.)) ? Real(1.) : Real(1.) + Real(2.)*ecx(i  ,j  ,0);
        Real hmy = (ecy(i  ,j-1,0) == Real(1.)) ? Real(1.) : Real(1.) - Real(2.)*ecy(i  ,j-1,0);
        Real hpy = (ecy(i  ,j  ,0) == Real(1.)) ? Real(1.) : Real(1.) + Real(2.)*ecy(i  ,j  ,0);
        Real const s = amrex::min(hmx,hpx,hmy,hpy);
        rhs(i,j,0) *= s;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_scaled_hinv (Real scale, Real h) noexcept
{
    return (h == Real(0.0)) ? Real(1.0) : scale / h;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlebndfdlap_scaled_h2inv (Real scale, Real h) noexcept
{
    return (h == Real(0.0)) ? Real(1.0) : scale / (h*h);
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                Array4<Real const> const& x, Array4<Real const> const& levset,
                                Array4<int const> const& dmsk,
                                Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                F const& xeb, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
        Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
        Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy);

        Real out = Real(0.0);
        Real fac = bx * Real(2.0) / (hpx+hmx);
        if (levset(i+1,j,k) < Real(0.0)) {
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

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Real xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy,
                              [=] (int, int, int) -> Real { return xeb; },
                              bx, by);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Array4<Real const> const& xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy,
                              [=] (int i1, int i2, int i3) -> Real {
                                  return xeb(i1,i2,i3); },
                              bx, by);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int i, int j, int k, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<int const> const& dmsk,
                        Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        y(i,j,k) = bx * (x(i-1,j,k) + x(i+1,j,k))
            +      by * (x(i,j-1,k) + x(i,j+1,k))
            - (Real(2.0)*(bx+by)) * x(i,j,k);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<Real const> const& levset,
                          Array4<int const> const& dmsk,
                          Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                          Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy);

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

            Real Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb (int i, int j, int k, Array4<Real> const& x,
                       Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                       Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real gamma = Real(-2.0)*(bx+by);
            Real Ax = bx * (x(i-1,j,k) + x(i+1,j,k))
                +     by * (x(i,j-1,k) + x(i,j+1,k))
                + gamma * x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

// RZ

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                   Array4<Real const> const& x, Array4<Real const> const& levset,
                                   Array4<int const> const& dmsk,
                                   Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                   F const& xeb, Real sigr, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i) * dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real tmp;
        Real hp, hm, out;

        Real const hpz = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmz = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real scale;

        if (r == Real(0.0)) {
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                scale = amrex::min(Real(1.0), hmz, hpz);
                out = scale * Real(4.0) * sigr * (x(i+1,j,k)-x(i,j,k)) / (dr*dr);
            } else {
                hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                scale = amrex::min(hp, hmz, hpz);
                out = Real(4.0) * sigr * mlebndfdlap_scaled_h2inv(scale, hp)
                    * (xeb(i+1,j,k)-x(i,j,k)) / (dr*dr);
            }
        } else {
            hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            hm = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            scale = amrex::min(hm, hp, hmz, hpz);
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                tmp = scale * (x(i+1,j,k) - x(i,j,k)) * (r + Real(0.5) * dr);
            } else {
                tmp = (xeb(i+1,j,k) - x(i,j,k)) * mlebndfdlap_scaled_hinv(scale, hp)
                    * (r + Real(0.5) * hp * dr);
            }

            if (levset(i-1,j,k) < Real(0.0)) {
                tmp += scale * (x(i-1,j,k) - x(i,j,k)) * (r - Real(0.5) * dr);
            } else {
                tmp += (xeb(i-1,j,k) - x(i,j,k)) * mlebndfdlap_scaled_hinv(scale, hm)
                    * (r - Real(0.5) * hm * dr);
            }

            out = tmp * Real(2.0) * sigr / ((hp+hm) * r * dr * dr);
        }

        if (levset(i,j+1,k) < Real(0.0)) {
            tmp = scale * (x(i,j+1,k) - x(i,j,k));
        } else {
            tmp = (xeb(i,j+1,k) - x(i,j,k)) * mlebndfdlap_scaled_hinv(scale, hpz);
        }


        if (levset(i,j-1,k) < Real(0.0)) {
            tmp += scale * (x(i,j-1,k) - x(i,j,k));
        } else {
            tmp += (xeb(i,j-1,k) - x(i,j,k)) * mlebndfdlap_scaled_hinv(scale, hmz);
        }

        out += tmp * Real(2.0) / ((hpz+hmz) * dz * dz);

        if (r != Real(0.0)) {
            out -= scale * alpha/(r*r) * x(i,j,k);
        }

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                              Array4<Real const> const& x, Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                              Real xeb, Real sigr, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    mlebndfdlap_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy,
                                 [=] (int, int, int) -> Real { return xeb; },
                                 sigr, dr, dz, rlo, alpha);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                              Array4<Real const> const& x, Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                              Array4<Real const> const& xeb, Real sigr, Real dr, Real dz, Real rlo,
                              Real alpha) noexcept
{
    mlebndfdlap_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy,
                                 [=] (int i1, int i2, int i3) -> Real {
                                     return xeb(i1,i2,i3); },
                                 sigr, dr, dz, rlo, alpha);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<int const> const& dmsk,
                           Real sigr, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i)*dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real Ax = (x(i,j-1,k) - Real(2.0)*x(i,j,k) + x(i,j+1,k)) / (dz*dz);
        if (r == Real(0.0)) {
            Ax += Real(4.0) * sigr * (x(i+1,j,k)-x(i,j,k)) / (dr*dr);
        } else {
            Real const rp = r + Real(0.5)*dr;
            Real const rm = r - Real(0.5)*dr;
            Ax += sigr * (rp*x(i+1,j,k) - (rp+rm)*x(i,j,k) + rm*x(i-1,j,k)) / (r*dr*dr);
            Ax -= alpha/(r*r) * x(i,j,k);
        }
        y(i,j,k) = Ax;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_rz_eb (int i, int j, int k, Array4<Real> const& x,
                             Array4<Real const> const& rhs, Array4<Real const> const& levset,
                             Array4<int const> const& dmsk,
                             Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                             Real sigr, Real dr, Real dz, Real rlo, int redblack, Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i) * dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real tmp, tmp0;
            Real hp, hm, Ax, gamma;

            Real const hpz = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmz = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real scale;

            if (r == Real(0.0)) {
                if (levset(i+1,j,k) < Real(0.0)) { // regular
                    scale = amrex::min(Real(1.0), hmz, hpz);
                    Ax = (scale * Real(4.0) * sigr / (dr*dr)) * (x(i+1,j,k)-x(i,j,k));
                    gamma = -(scale * Real(4.0) * sigr / (dr*dr));
                } else {
                    hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                    scale = amrex::min(hp, hmz, hpz);
                    gamma = -(Real(4.0) * sigr * mlebndfdlap_scaled_h2inv(scale, hp) / (dr*dr));
                    Ax = gamma * x(i,j,k);
                }
            } else {
                hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                hm = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
                scale = amrex::min(hm, hp, hmz, hpz);
                if (levset(i+1,j,k) < Real(0.0)) { // regular
                    tmp = scale * (x(i+1,j,k) - x(i,j,k)) * (r + Real(0.5) * dr);
                    tmp0 = -scale * (r + Real(0.5) * dr);
                } else {
                    tmp0 = -mlebndfdlap_scaled_hinv(scale, hp)
                        * (r + Real(0.5) * hp * dr);
                    tmp = tmp0 * x(i,j,k);
                }

                if (levset(i-1,j,k) < Real(0.0)) {
                    tmp += scale * (x(i-1,j,k) - x(i,j,k)) * (r - Real(0.5) * dr);
                    tmp0 += -scale * (r - Real(0.5) * dr);
                } else {
                    tmp += -x(i,j,k) * mlebndfdlap_scaled_hinv(scale, hm)
                        * (r - Real(0.5) * hm * dr);
                    tmp0 += -mlebndfdlap_scaled_hinv(scale, hm)
                        * (r - Real(0.5) * hm * dr);
                }

                Ax = tmp * Real(2.0) * sigr / ((hp+hm) * r * dr * dr);
                gamma = tmp0 * Real(2.0) * sigr / ((hp+hm) * r * dr * dr);
            }

            if (levset(i,j+1,k) < Real(0.0)) {
                tmp = scale * (x(i,j+1,k) - x(i,j,k));
                tmp0 = -scale;
            } else {
                tmp0 = -mlebndfdlap_scaled_hinv(scale, hpz);
                tmp = tmp0 * x(i,j,k);
            }

            if (levset(i,j-1,k) < Real(0.0)) {
                tmp += scale * (x(i,j-1,k) - x(i,j,k));
                tmp0 += -scale;
            } else {
                tmp += -x(i,j,k) * mlebndfdlap_scaled_hinv(scale, hmz);
                tmp0 += -mlebndfdlap_scaled_hinv(scale, hmz);
            }

            Ax += tmp * Real(2.0) / ((hpz+hmz) * dz * dz);
            gamma += tmp0 * Real(2.0) / ((hpz+hmz) * dz * dz);

            if (r != Real(0.0)) {
                Ax -= scale * alpha/(r*r) * x(i,j,k);
                gamma -= scale * alpha/(r*r);
            }

            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_rz (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                          Real sigr, Real dr, Real dz, Real rlo, int redblack, Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i)*dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real Ax = (x(i,j-1,k) - Real(2.0)*x(i,j,k) + x(i,j+1,k)) / (dz*dz);
            Real gamma = -Real(2.0) / (dz*dz);
            if (r == Real(0.0)) {
                Ax += (Real(4.0)*sigr/(dr*dr)) * (x(i+1,j,k)-x(i,j,k));
                gamma += -(Real(4.0)*sigr/(dr*dr));
            } else {
                Real const rp = r + Real(0.5)*dr;
                Real const rm = r - Real(0.5)*dr;
                Ax += sigr*(rp*x(i+1,j,k) - (rp+rm)*x(i,j,k) + rm*x(i-1,j,k)) / (r*dr*dr);
                gamma += -sigr*(rp+rm) / (r*dr*dr);
                Ax -= alpha/(r*r) * x(i,j,k);
                gamma -= alpha/(r*r);
            }
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
                            Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real sigxm = Real(0.5)*(sig(i-1,j-1,k)+sig(i-1,j  ,k));
        Real sigxp = Real(0.5)*(sig(i  ,j-1,k)+sig(i  ,j  ,k));
        Real sigym = Real(0.5)*(sig(i-1,j-1,k)+sig(i  ,j-1,k));
        Real sigyp = Real(0.5)*(sig(i-1,j  ,k)+sig(i  ,j  ,k));
        y(i,j,k) = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
            +      by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
            - (bx*(sigxm+sigxp) + by*(sigym+sigyp)) * x(i,j,k);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb (int i, int j, int k, Array4<Real> const& x,
                           Array4<Real const> const& rhs,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& sig,
                           Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real sigxm = Real(0.5)*(sig(i-1,j-1,k)+sig(i-1,j  ,k));
            Real sigxp = Real(0.5)*(sig(i  ,j-1,k)+sig(i  ,j  ,k));
            Real sigym = Real(0.5)*(sig(i-1,j-1,k)+sig(i  ,j-1,k));
            Real sigyp = Real(0.5)*(sig(i-1,j  ,k)+sig(i  ,j  ,k));
            Real gamma = -(bx*(sigxm+sigxp) + by*(sigym+sigyp));
            Real Ax = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
                +     by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
                + gamma * x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_rz (int i, int j, int k, Array4<Real> const& y,
                                Array4<Real const> const& x,
                                Array4<int const> const& dmsk,
                                Array4<Real const> const& sig,
                                Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i)*dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const rp = r + Real(0.5)*dr;
        Real const rm = (r >= Real(0.5)*dr) ? (r - Real(0.5)*dr) : (Real(0.5)*dr - r);
        Real const sigyp = (sig(i-1,j  ,k)*rm + sig(i,j  ,k)*rp) / (rm+rp);
        Real const sigym = (sig(i-1,j-1,k)*rm + sig(i,j-1,k)*rp) / (rm+rp);
        Real Ax = (sigyp*(x(i,j+1,k)-x(i,j,k)) + sigym*(x(i,j-1,k)-x(i,j,k))) / (dz*dz);
        if (r == Real(0.0)) {
            Real const sigxp = Real(0.5)*(sig(i,j-1,k)+sig(i,j,k));
            Ax += Real(4.0) * sigxp * (x(i+1,j,k)-x(i,j,k)) / (dr*dr);
        } else {
            Real const sigxp = Real(0.5)*(sig(i,  j-1,k)+sig(i,  j,k));
            Real const sigxm = Real(0.5)*(sig(i-1,j-1,k)+sig(i-1,j,k));
            Ax += (sigxp*rp*(x(i+1,j,k)-x(i,j,k)) + sigxm*rm*(x(i-1,j,k)-x(i,j,k))) / (r*dr*dr);
            Ax -= alpha/(r*r) * x(i,j,k);
        }
        y(i,j,k) = Ax;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb_rz (int i, int j, int k, Array4<Real> const& x,
                               Array4<Real const> const& rhs,
                               Array4<int const> const& dmsk,
                               Array4<Real const> const& sig,
                               Real dr, Real dz, Real rlo, int redblack, Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i)*dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real const rp = r + Real(0.5)*dr;
            Real const rm = (r >= Real(0.5)*dr) ? (r - Real(0.5)*dr) : (Real(0.5)*dr - r);
            Real const sigyp = (sig(i-1,j  ,k)*rm + sig(i,j  ,k)*rp) / (rm+rp);
            Real const sigym = (sig(i-1,j-1,k)*rm + sig(i,j-1,k)*rp) / (rm+rp);
            Real Ax = (sigyp*(x(i,j+1,k)-x(i,j,k)) + sigym*(x(i,j-1,k)-x(i,j,k))) / (dz*dz);
            Real gamma = -(sigyp+sigym) / (dz*dz);
            if (r == Real(0.0)) {
                Real const sigxp = Real(0.5)*(sig(i,j-1,k)+sig(i,j,k));
                Ax += (Real(4.0)*sigxp/(dr*dr)) * (x(i+1,j,k)-x(i,j,k));
                gamma += -(Real(4.0)*sigxp/(dr*dr));
            } else {
                Real const sigxp = Real(0.5)*(sig(i,  j-1,k)+sig(i,  j,k));
                Real const sigxm = Real(0.5)*(sig(i-1,j-1,k)+sig(i-1,j,k));
                Ax += (sigxp*rp*(x(i+1,j,k)-x(i,j,k)) + sigxm*rm*(x(i-1,j,k)-x(i,j,k))) / (r*dr*dr);
                gamma += -(sigxp*rp+sigxm*rm) / (r*dr*dr);
                Ax -= alpha/(r*r) * x(i,j,k);
                gamma -= alpha/(r*r);
            }
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                    Array4<Real const> const& x, Array4<Real const> const& levset,
                                    Array4<int const> const& dmsk,
                                    Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                    Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                                    F const& xeb, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
        Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
        Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy);

        Real out = Real(0.0);
        Real sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
            /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
        Real fac = bx * Real(2.0) / (hpx+hmx);
        if (levset(i+1,j,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i+1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx) * (xeb(i+1,j,k) - xc);
        }

        sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i-1,j,k)*vfrc(i-1,j,k))
            /  (vfrc(i-1,j-1,k)                + vfrc(i-1,j,k));
        if (levset(i-1,j,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i-1,j,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx) * (xeb(i-1,j,k) - xc);
        }

        sigma = (sig(i-1,j,k)*vfrc(i-1,j,k) + sig(i,j,k)*vfrc(i,j,k))
            /  (vfrc(i-1,j,k)              + vfrc(i,j,k));
        fac = by * Real(2.0) / (hpy+hmy);
        if (levset(i,j+1,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j+1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy) * (xeb(i,j+1,k) - xc);
        }

        sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i,j-1,k)*vfrc(i,j-1,k))
            /  (vfrc(i-1,j-1,k)                + vfrc(i,j-1,k));
        if (levset(i,j-1,k) < Real(0.0)) {
            out += sigma * fac * scale * (x(i,j-1,k) - xc);
        } else {
            out += sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy) * (xeb(i,j-1,k) - xc);
        }

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                               Array4<Real const> const& x, Array4<Real const> const& levset,
                               Array4<int const> const& dmsk,
                               Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                               Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                               Real xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_sig_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, sig, vfrc,
                                  [=] (int, int, int) -> Real { return xeb; },
                                  bx, by);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                               Array4<Real const> const& x, Array4<Real const> const& levset,
                               Array4<int const> const& dmsk,
                               Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                               Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                               Array4<Real const> const& xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_sig_adotx_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, sig, vfrc,
                                  [=] (int i1, int i2, int i3) -> Real {
                                      return xeb(i1,i2,i3); },
                                  bx, by);
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_rz_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                        Array4<Real const> const& x, Array4<Real const> const& levset,
                                        Array4<int const> const& dmsk,
                                        Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                        Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                                        F const& xeb, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i) * dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real tmp, sigma;
        Real hp, hm, out;

        Real const hpz = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
        Real const hmz = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
        Real scale;

        if (r == Real(0.0)) {
            sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
                /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                scale = amrex::min(Real(1.0), hmz, hpz);
                out = scale * Real(4.0) * sigma * (x(i+1,j,k)-x(i,j,k)) / (dr*dr);
            } else {
                hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                scale = amrex::min(hp, hmz, hpz);
                out = Real(4.0) * sigma * mlebndfdlap_scaled_h2inv(scale, hp)
                    * (xeb(i+1,j,k)-x(i,j,k)) / (dr*dr);
            }
        } else {
            sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
                /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
            hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            hm = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            scale = amrex::min(hm, hp, hmz, hpz);
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                tmp = scale * sigma * (x(i+1,j,k) - x(i,j,k)) * (r + Real(0.5) * dr);
            } else {
                tmp = sigma * (xeb(i+1,j,k) - x(i,j,k))
                    * mlebndfdlap_scaled_hinv(scale, hp) * (r + Real(0.5) * hp * dr);
            }

            sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i-1,j,k)*vfrc(i-1,j,k))
                /  (vfrc(i-1,j-1,k)                + vfrc(i-1,j,k));
            if (levset(i-1,j,k) < Real(0.0)) {
                tmp += scale * sigma * (x(i-1,j,k) - x(i,j,k)) * (r - Real(0.5) * dr);
            } else {
                tmp += sigma * (xeb(i-1,j,k) - x(i,j,k))
                    * mlebndfdlap_scaled_hinv(scale, hm) * (r - Real(0.5) * hm * dr);
            }

            out = tmp * Real(2.0) / ((hp+hm) * r * dr * dr);
        }

        // Cell-centre radii for cells at i and i-1; abs handles the axis (r=0) where
        // the cell at i-1 is the reflected image at the same physical radius as i.
        Real const rp = r + Real(0.5) * dr;
        Real const rm = (r >= Real(0.5) * dr) ? (r - Real(0.5) * dr) : (Real(0.5) * dr - r);

        sigma = (sig(i-1,j,k)*vfrc(i-1,j,k)*rm + sig(i,j,k)*vfrc(i,j,k)*rp)
            /  (vfrc(i-1,j,k)*rm              + vfrc(i,j,k)*rp);
        if (levset(i,j+1,k) < Real(0.0)) {
            tmp = scale * sigma * (x(i,j+1,k) - x(i,j,k));
        } else {
            tmp = sigma * (xeb(i,j+1,k) - x(i,j,k))
                * mlebndfdlap_scaled_hinv(scale, hpz);
        }

        sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k)*rm + sig(i,j-1,k)*vfrc(i,j-1,k)*rp)
            /  (vfrc(i-1,j-1,k)*rm                + vfrc(i,j-1,k)*rp);
        if (levset(i,j-1,k) < Real(0.0)) {
            tmp += scale * sigma * (x(i,j-1,k) - x(i,j,k));
        } else {
            tmp += sigma * (xeb(i,j-1,k) - x(i,j,k))
                * mlebndfdlap_scaled_hinv(scale, hmz);
        }

        out += tmp * Real(2.0) / ((hpz+hmz) * dz * dz);

        if (r != Real(0.0)) {
            out -= scale * alpha/(r*r) * x(i,j,k);
        }

        y(i,j,k) = out;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                                   Array4<Real const> const& x, Array4<Real const> const& levset,
                                   Array4<int const> const& dmsk,
                                   Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                   Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                                   Real xeb, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    mlebndfdlap_sig_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, sig, vfrc,
                                      [=] (int, int, int) -> Real { return xeb; },
                                      dr, dz, rlo, alpha);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                                   Array4<Real const> const& x, Array4<Real const> const& levset,
                                   Array4<int const> const& dmsk,
                                   Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                   Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                                   Array4<Real const> const& xeb, Real dr, Real dz, Real rlo,
                                   Real alpha) noexcept
{
    mlebndfdlap_sig_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, ecx, ecy, sig, vfrc,
                                      [=] (int i1, int i2, int i3) -> Real {
                                          return xeb(i1,i2,i3); },
                                      dr, dz, rlo, alpha);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                              Array4<Real const> const& rhs,
                              Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx,
                              Array4<Real const> const& ecy,
                              Array4<Real const> const& sig,
                              Array4<Real const> const& vfrc,
                              Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            Real const hmx = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            Real const hpy = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmy = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy);

            Real gamma = Real(0.0);
            Real rho = Real(0.0);
            Real sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
                /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
            Real fac = bx * Real(2.0) / (hpx+hmx);
            if (levset(i+1,j,k) < Real(0.0)) { // regular
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i+1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpx);
            }

            sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i-1,j,k)*vfrc(i-1,j,k))
                /  (vfrc(i-1,j-1,k)                + vfrc(i-1,j,k));
            if (levset(i-1,j,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i-1,j,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmx);
            }

            sigma = (sig(i-1,j,k)*vfrc(i-1,j,k) + sig(i,j,k)*vfrc(i,j,k))
                /  (vfrc(i-1,j,k)              + vfrc(i,j,k));
            fac = by * Real(2.0) / (hpy+hmy);
            if (levset(i,j+1,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j+1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hpy);
            }

            sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i,j-1,k)*vfrc(i,j-1,k))
                /  (vfrc(i-1,j-1,k)                + vfrc(i,j-1,k));
            if (levset(i,j-1,k) < Real(0.0)) {
                gamma -= sigma * fac * scale;
                rho += sigma * fac * scale * x(i,j-1,k);
            } else {
                gamma -= sigma * fac * mlebndfdlap_scaled_hinv(scale, hmy);
            }

            Real Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb_rz_eb (int i, int j, int k, Array4<Real> const& x,
                                  Array4<Real const> const& rhs, Array4<Real const> const& levset,
                                  Array4<int const> const& dmsk,
                                  Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                  Array4<Real const> const& sig, Array4<Real const> const& vfrc,
                                  Real dr, Real dz, Real rlo, int redblack, Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i) * dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real tmp0, tmp1, sigma;
            Real hp, hm, Ax, gamma;

            Real const hpz = (ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            Real const hmz = (ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            Real scale;

            if (r == Real(0.0)) {
                sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
                    /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
                if (levset(i+1,j,k) < Real(0.0)) { // regular
                    scale = amrex::min(Real(1.0), hmz, hpz);
                    Ax = (scale * Real(4.0) * sigma / (dr*dr)) * (x(i+1,j,k)-x(i,j,k));
                    gamma = -(scale * Real(4.0) * sigma / (dr*dr));
                } else {
                    hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                    scale = amrex::min(hp, hmz, hpz);
                    gamma = -(Real(4.0) * sigma * mlebndfdlap_scaled_h2inv(scale, hp) / (dr*dr));
                    Ax = gamma * x(i,j,k);
                }
            } else {
                sigma = (sig(i,j-1,k)*vfrc(i,j-1,k) + sig(i,j,k)*vfrc(i,j,k))
                    /  (vfrc(i,j-1,k)              + vfrc(i,j,k));
                hp = (ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
                hm = (ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
                scale = amrex::min(hm, hp, hmz, hpz);
                if (levset(i+1,j,k) < Real(0.0)) { // regular
                    tmp1 = scale * sigma * x(i+1,j,k) * (r + Real(0.5) * dr);
                    tmp0 = scale * sigma * (-(r + Real(0.5) * dr));
                } else {
                    tmp0 = -sigma * mlebndfdlap_scaled_hinv(scale, hp)
                        * (r + Real(0.5) * hp * dr);
                    tmp1 = Real(0.0);
                }

                sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k) + sig(i-1,j,k)*vfrc(i-1,j,k))
                    /  (vfrc(i-1,j-1,k)                + vfrc(i-1,j,k));
                if (levset(i-1,j,k) < Real(0.0)) {
                    tmp1 += scale * sigma * x(i-1,j,k) * (r - Real(0.5) * dr);
                    tmp0 += scale * sigma * (-(r - Real(0.5) * dr));
                } else {
                    tmp0 += -sigma * mlebndfdlap_scaled_hinv(scale, hm)
                        * (r - Real(0.5) * hm * dr);
                }

                Ax = (tmp1 + tmp0*x(i,j,k)) * Real(2.0) / ((hp+hm) * r * dr * dr);
                gamma = tmp0 * Real(2.0) / ((hp+hm) * r * dr * dr);
            }

            Real const rp = r + Real(0.5) * dr;
            Real const rm = (r >= Real(0.5) * dr) ? (r - Real(0.5) * dr) : (Real(0.5) * dr - r);

            sigma = (sig(i-1,j,k)*vfrc(i-1,j,k)*rm + sig(i,j,k)*vfrc(i,j,k)*rp)
                /  (vfrc(i-1,j,k)*rm              + vfrc(i,j,k)*rp);
            if (levset(i,j+1,k) < Real(0.0)) {
                tmp0 = -scale * sigma;
                tmp1 = scale * sigma * x(i,j+1,k);
            } else {
                tmp0 = -sigma * mlebndfdlap_scaled_hinv(scale, hpz);
                tmp1 = Real(0.0);
            }

            sigma = (sig(i-1,j-1,k)*vfrc(i-1,j-1,k)*rm + sig(i,j-1,k)*vfrc(i,j-1,k)*rp)
                /  (vfrc(i-1,j-1,k)*rm                + vfrc(i,j-1,k)*rp);
            if (levset(i,j-1,k) < Real(0.0)) {
                tmp0 += -scale * sigma;
                tmp1 += scale * sigma * x(i,j-1,k);
            } else {
                tmp0 += -sigma * mlebndfdlap_scaled_hinv(scale, hmz);
            }

            Ax += (tmp1 + tmp0*x(i,j,k)) * Real(2.0) / ((hpz+hmz) * dz * dz);
            gamma += tmp0 * Real(2.0) / ((hpz+hmz) * dz * dz);

            if (r != Real(0.0)) {
                Ax -= scale * alpha/(r*r) * x(i,j,k);
                gamma -= scale * alpha/(r*r);
            }

            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

}

#endif
