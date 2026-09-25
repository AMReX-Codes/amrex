#ifndef AMREX_MLEBNODEFDLAP_2D_K_H_
#define AMREX_MLEBNODEFDLAP_2D_K_H_

namespace amrex {

// Edge data convention used by these kernels: elx(i,j,k) is the open length
// of the x-edge between nodes i and i+1 measured from node i, and
// mlebndfdlap_hm(elx(i,j,k)) is the open length measured from node i+1.
// An edge is open all the way to a neighbor only if the length is 1 and the
// neighbor node is not covered (levset < 0).

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_scale_rhs (int i, int j, int, Array4<Real> const& rhs,
                            Array4<int const> const& dmsk, Array4<Real const> const& elx,
                            Array4<Real const> const& ely) noexcept
{
    if (!dmsk(i,j,0)) {
        Real const s = amrex::min(mlebndfdlap_hm(elx(i-1,j,0)), elx(i,j,0),
                                  mlebndfdlap_hm(ely(i,j-1,0)), ely(i,j,0));
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

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int i, int j, int k, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<int const> const& dmsk,
                        S const& sig, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const sigxm = sig.x(i-1,j,k);
        Real const sigxp = sig.x(i  ,j,k);
        Real const sigym = sig.y(i,j-1,k);
        Real const sigyp = sig.y(i,j  ,k);
        y(i,j,k) = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
            +      by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
            - (bx*(sigxm+sigxp) + by*(sigym+sigyp)) * x(i,j,k);
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb (int i, int j, int k, Array4<Real> const& x,
                       Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                       S const& sig, Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const sigxm = sig.x(i-1,j,k);
            Real const sigxp = sig.x(i  ,j,k);
            Real const sigym = sig.y(i,j-1,k);
            Real const sigyp = sig.y(i,j  ,k);
            Real const gamma = -(bx*(sigxm+sigxp) + by*(sigym+sigyp));
            Real const Ax = bx * (sigxm*x(i-1,j,k) + sigxp*x(i+1,j,k))
                +           by * (sigym*x(i,j-1,k) + sigyp*x(i,j+1,k))
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
                                Array4<Real const> const& elx, Array4<Real const> const& ely,
                                S const& sig, F const& xeb, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hpx = elx(i  ,j  ,k);
        Real const hmx = mlebndfdlap_hm(elx(i-1,j  ,k));
        Real const hpy = ely(i  ,j  ,k);
        Real const hmy = mlebndfdlap_hm(ely(i  ,j-1,k));
        Real const scale = amrex::min(hmx,hpx,hmy,hpy);

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

        y(i,j,k) = out;
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& elx, Array4<Real const> const& ely,
                           S const& sig, Real xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, elx, ely, sig,
                              [=] (int, int, int) -> Real { return xeb; },
                              bx, by);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<Real const> const& levset,
                           Array4<int const> const& dmsk,
                           Array4<Real const> const& elx, Array4<Real const> const& ely,
                           S const& sig, Array4<Real const> const& xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, levset, dmsk, elx, ely, sig,
                              [=] (int i1, int i2, int i3) -> Real {
                                  return xeb(i1,i2,i3); },
                              bx, by);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_eb (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<Real const> const& levset,
                          Array4<int const> const& dmsk,
                          Array4<Real const> const& elx, Array4<Real const> const& ely,
                          S const& sig, Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real const hpx = elx(i  ,j  ,k);
            Real const hmx = mlebndfdlap_hm(elx(i-1,j  ,k));
            Real const hpy = ely(i  ,j  ,k);
            Real const hmy = mlebndfdlap_hm(ely(i  ,j-1,k));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy);

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

            Real const Ax = rho + gamma*x(i,j,k);
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

// RZ

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<int const> const& dmsk,
                           S const& sig, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i)*dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const sigyp = sig.y(i,j  ,k);
        Real const sigym = sig.y(i,j-1,k);
        Real Ax = (sigyp*(x(i,j+1,k)-xc) + sigym*(x(i,j-1,k)-xc)) / (dz*dz);
        if (r == Real(0.0)) {
            Real const sigxp = sig.x(i,j,k);
            Ax += Real(4.0) * sigxp * (x(i+1,j,k)-xc) / (dr*dr);
        } else {
            Real const rp = r + Real(0.5)*dr;
            Real const rm = r - Real(0.5)*dr;
            Real const sigxp = sig.x(i  ,j,k);
            Real const sigxm = sig.x(i-1,j,k);
            Ax += (sigxp*rp*(x(i+1,j,k)-xc) + sigxm*rm*(x(i-1,j,k)-xc)) / (r*dr*dr);
            Ax -= alpha/(r*r) * xc;
        }
        y(i,j,k) = Ax;
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_rz (int i, int j, int k, Array4<Real> const& x,
                          Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                          S const& sig, Real dr, Real dz, Real rlo, int redblack,
                          Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i)*dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real const xc = x(i,j,k);
            Real const sigyp = sig.y(i,j  ,k);
            Real const sigym = sig.y(i,j-1,k);
            Real Ax = (sigyp*(x(i,j+1,k)-xc) + sigym*(x(i,j-1,k)-xc)) / (dz*dz);
            Real gamma = -(sigyp+sigym) / (dz*dz);
            if (r == Real(0.0)) {
                Real const sigxp = sig.x(i,j,k);
                Ax += (Real(4.0)*sigxp/(dr*dr)) * (x(i+1,j,k)-xc);
                gamma += -(Real(4.0)*sigxp/(dr*dr));
            } else {
                Real const rp = r + Real(0.5)*dr;
                Real const rm = r - Real(0.5)*dr;
                Real const sigxp = sig.x(i  ,j,k);
                Real const sigxm = sig.x(i-1,j,k);
                Ax += (sigxp*rp*(x(i+1,j,k)-xc) + sigxm*rm*(x(i-1,j,k)-xc)) / (r*dr*dr);
                gamma += -(sigxp*rp+sigxm*rm) / (r*dr*dr);
                Ax -= alpha/(r*r) * xc;
                gamma -= alpha/(r*r);
            }
            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

template <typename S, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                   Array4<Real const> const& x, Array4<Real const> const& levset,
                                   Array4<int const> const& dmsk,
                                   Array4<Real const> const& elx, Array4<Real const> const& ely,
                                   S const& sig, F const& xeb,
                                   Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    Real const r = rlo + Real(i) * dr;
    if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
        y(i,j,k) = Real(0.0);
    } else {
        Real const xc = x(i,j,k);
        Real const hp  = elx(i  ,j  ,k);
        Real const hm  = mlebndfdlap_hm(elx(i-1,j  ,k));
        Real const hpz = ely(i  ,j  ,k);
        Real const hmz = mlebndfdlap_hm(ely(i  ,j-1,k));
        Real scale, out, tmp;

        if (r == Real(0.0)) {
            Real const sigp = sig.x(i,j,k);
            if (levset(i+1,j,k) < Real(0.0) && hp == Real(1.0)) { // regular
                scale = amrex::min(Real(1.0), hmz, hpz);
                out = scale * Real(4.0) * sigp * (x(i+1,j,k)-xc) / (dr*dr);
            } else {
                scale = amrex::min(hp, hmz, hpz);
                out = Real(4.0) * sigp * mlebndfdlap_scaled_h2inv(scale, hp)
                    * (xeb(i+1,j,k)-xc) / (dr*dr);
            }
        } else {
            scale = amrex::min(hm, hp, hmz, hpz);
            Real const sigp = sig.x(i  ,j,k);
            Real const sigm = sig.x(i-1,j,k);
            if (levset(i+1,j,k) < Real(0.0) && hp == Real(1.0)) { // regular
                tmp = scale * sigp * (x(i+1,j,k) - xc) * (r + Real(0.5) * dr);
            } else {
                tmp = sigp * (xeb(i+1,j,k) - xc) * mlebndfdlap_scaled_hinv(scale, hp)
                    * (r + Real(0.5) * hp * dr);
            }
            if (levset(i-1,j,k) < Real(0.0) && hm == Real(1.0)) {
                tmp += scale * sigm * (x(i-1,j,k) - xc) * (r - Real(0.5) * dr);
            } else {
                tmp += sigm * (xeb(i-1,j,k) - xc) * mlebndfdlap_scaled_hinv(scale, hm)
                    * (r - Real(0.5) * hm * dr);
            }
            out = tmp * Real(2.0) / ((hp+hm) * r * dr * dr);
        }

        Real const sigp = sig.y(i,j  ,k);
        Real const sigm = sig.y(i,j-1,k);
        if (levset(i,j+1,k) < Real(0.0) && hpz == Real(1.0)) {
            tmp = scale * sigp * (x(i,j+1,k) - xc);
        } else {
            tmp = sigp * (xeb(i,j+1,k) - xc) * mlebndfdlap_scaled_hinv(scale, hpz);
        }
        if (levset(i,j-1,k) < Real(0.0) && hmz == Real(1.0)) {
            tmp += scale * sigm * (x(i,j-1,k) - xc);
        } else {
            tmp += sigm * (xeb(i,j-1,k) - xc) * mlebndfdlap_scaled_hinv(scale, hmz);
        }
        out += tmp * Real(2.0) / ((hpz+hmz) * dz * dz);

        if (r != Real(0.0)) {
            out -= scale * alpha/(r*r) * xc;
        }

        y(i,j,k) = out;
    }
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                              Array4<Real const> const& x, Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& elx, Array4<Real const> const& ely,
                              S const& sig, Real xeb, Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    mlebndfdlap_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, elx, ely, sig,
                                 [=] (int, int, int) -> Real { return xeb; },
                                 dr, dz, rlo, alpha);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_rz_eb (int i, int j, int k, Array4<Real> const& y,
                              Array4<Real const> const& x, Array4<Real const> const& levset,
                              Array4<int const> const& dmsk,
                              Array4<Real const> const& elx, Array4<Real const> const& ely,
                              S const& sig, Array4<Real const> const& xeb,
                              Real dr, Real dz, Real rlo, Real alpha) noexcept
{
    mlebndfdlap_adotx_rz_eb_doit(i, j, k, y, x, levset, dmsk, elx, ely, sig,
                                 [=] (int i1, int i2, int i3) -> Real {
                                     return xeb(i1,i2,i3); },
                                 dr, dz, rlo, alpha);
}

template <typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_rz_eb (int i, int j, int k, Array4<Real> const& x,
                             Array4<Real const> const& rhs, Array4<Real const> const& levset,
                             Array4<int const> const& dmsk,
                             Array4<Real const> const& elx, Array4<Real const> const& ely,
                             S const& sig, Real dr, Real dz, Real rlo, int redblack,
                             Real alpha) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        Real const r = rlo + Real(i) * dr;
        if (dmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
            x(i,j,k) = Real(0.);
        } else {
            Real const xc = x(i,j,k);
            Real const hp  = elx(i  ,j  ,k);
            Real const hm  = mlebndfdlap_hm(elx(i-1,j  ,k));
            Real const hpz = ely(i  ,j  ,k);
            Real const hmz = mlebndfdlap_hm(ely(i  ,j-1,k));
            // Ax = tmp1 + tmp0*xc for each direction, gamma is the diagonal.
            Real scale, Ax, gamma, tmp0, tmp1;

            if (r == Real(0.0)) {
                Real const sigp = sig.x(i,j,k);
                if (levset(i+1,j,k) < Real(0.0) && hp == Real(1.0)) { // regular
                    scale = amrex::min(Real(1.0), hmz, hpz);
                    gamma = -(scale * Real(4.0) * sigp / (dr*dr));
                    Ax = gamma * (xc - x(i+1,j,k));
                } else {
                    scale = amrex::min(hp, hmz, hpz);
                    gamma = -(Real(4.0) * sigp * mlebndfdlap_scaled_h2inv(scale, hp) / (dr*dr));
                    Ax = gamma * xc;
                }
            } else {
                scale = amrex::min(hm, hp, hmz, hpz);
                Real const sigp = sig.x(i  ,j,k);
                Real const sigm = sig.x(i-1,j,k);
                if (levset(i+1,j,k) < Real(0.0) && hp == Real(1.0)) { // regular
                    tmp0 = -scale * sigp * (r + Real(0.5) * dr);
                    tmp1 = -tmp0 * x(i+1,j,k);
                } else {
                    tmp0 = -sigp * mlebndfdlap_scaled_hinv(scale, hp)
                        * (r + Real(0.5) * hp * dr);
                    tmp1 = Real(0.0);
                }
                if (levset(i-1,j,k) < Real(0.0) && hm == Real(1.0)) {
                    Real const c = -scale * sigm * (r - Real(0.5) * dr);
                    tmp0 += c;
                    tmp1 -= c * x(i-1,j,k);
                } else {
                    tmp0 += -sigm * mlebndfdlap_scaled_hinv(scale, hm)
                        * (r - Real(0.5) * hm * dr);
                }
                Real const f = Real(2.0) / ((hp+hm) * r * dr * dr);
                Ax = (tmp1 + tmp0*xc) * f;
                gamma = tmp0 * f;
            }

            Real const sigp = sig.y(i,j  ,k);
            Real const sigm = sig.y(i,j-1,k);
            if (levset(i,j+1,k) < Real(0.0) && hpz == Real(1.0)) {
                tmp0 = -scale * sigp;
                tmp1 = scale * sigp * x(i,j+1,k);
            } else {
                tmp0 = -sigp * mlebndfdlap_scaled_hinv(scale, hpz);
                tmp1 = Real(0.0);
            }
            if (levset(i,j-1,k) < Real(0.0) && hmz == Real(1.0)) {
                tmp0 += -scale * sigm;
                tmp1 += scale * sigm * x(i,j-1,k);
            } else {
                tmp0 += -sigm * mlebndfdlap_scaled_hinv(scale, hmz);
            }
            Real const f = Real(2.0) / ((hpz+hmz) * dz * dz);
            Ax += (tmp1 + tmp0*xc) * f;
            gamma += tmp0 * f;

            if (r != Real(0.0)) {
                Ax -= scale * alpha/(r*r) * xc;
                gamma -= scale * alpha/(r*r);
            }

            constexpr Real omega = Real(1.25);
            x(i,j,k) += (rhs(i,j,k) - Ax) * (omega / gamma);
        }
    }
}

}

#endif
