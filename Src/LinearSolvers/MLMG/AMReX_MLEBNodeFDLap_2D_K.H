#ifndef AMREX_MLEBNODEFDLAP_2D_K_H_
#define AMREX_MLEBNODEFDLAP_2D_K_H_

namespace amrex {

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb_doit (int i, int j, int k, Array4<Real> const& y,
                                Array4<Real const> const& x, Array4<int const> const& dmsk,
                                Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                                F && xeb, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else if (!ecx) {
        y(i,j,k) = bx * (x(i-1,j,k) + x(i+1,j,k))
            +      by * (x(i,j-1,k) + x(i,j+1,k))
            - Real(2.0)*(bx+by) * x(i,j,k);
    } else {
        Real tmp;
        Real hp, hm;
        if (ecx(i,j,k) == Real(1.0)) { // regular
            hp = Real(1.0);
            tmp = x(i+1,j,k) - x(i,j,k);
        } else {
            hp = Real(1.0) + Real(2.) * ecx(i,j,k);
            tmp = (xeb(i+1,j,k) - x(i,j,k)) / hp;
        }

        if (ecx(i-1,j,k) == Real(1.0)) {
            hm = Real(1.0);
            tmp += x(i-1,j,k) - x(i,j,k);
        } else {
            hm = Real(1.0) - Real(2.) * ecx(i-1,j,k);
            tmp += (xeb(i-1,j,k) - x(i,j,k)) / hm;
        }

        y(i,j,k) = tmp * bx * Real(2.0) / (hp+hm);

        if (ecy(i,j,k) == Real(1.0)) {
            hp = Real(1.0);
            tmp = x(i,j+1,k) - x(i,j,k);
        } else {
            hp = Real(1.0) + Real(2.) * ecy(i,j,k);
            tmp = (xeb(i,j+1,k) - x(i,j,k)) / hp;
        }

        if (ecy(i,j-1,k) == Real(1.0)) {
            hm = Real(1.0);
            tmp += x(i,j-1,k) - x(i,j,k);
        } else {
            hm = Real(1.0) - Real(2.) * ecy(i,j-1,k);
            tmp += (xeb(i,j-1,k) - x(i,j,k)) / hm;
        }

        y(i,j,k) += tmp * by * Real(2.0) / (hp+hm);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Real xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, dmsk, ecx, ecy,
                              [=] (int, int, int) -> Real { return xeb; },
                              bx, by);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx_eb (int i, int j, int k, Array4<Real> const& y,
                           Array4<Real const> const& x, Array4<int const> const& dmsk,
                           Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                           Array4<Real const> const& xeb, Real bx, Real by) noexcept
{
    mlebndfdlap_adotx_eb_doit(i, j, k, y, x, dmsk, ecx, ecy,
                              [=] (int i1, int i2, int i3) -> Real {
                                  return xeb(i1,i2,i3); },
                              bx, by);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int i, int j, int k, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<int const> const& dmsk,
                        Array4<Real const> const& a, Real bx, Real by) noexcept
{
    if (dmsk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        y(i,j,k) = bx * (x(i-1,j,k) + x(i+1,j,k))
            +      by * (x(i,j-1,k) + x(i,j+1,k))
            - (Real(2.0)*(bx+by) + a(i,j,k)) * x(i,j,k);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb_eb(int i, int j, int k, Array4<Real> const& x,
                         Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                         Array4<Real const> const& ecx, Array4<Real const> const& ecy,
                         Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else if (!ecx) {
            Real gamma = Real(-2.0)*(bx+by);
            Real rho = bx * (x(i-1,j,k) + x(i+1,j,k))
                +      by * (x(i,j-1,k) + x(i,j+1,k));
            x(i,j,k) = (rhs(i,j,k) - rho) / gamma;
        } else {
            Real tmp0, tmp1;
            Real hp, hm;
            if (ecx(i,j,k) == Real(1.0)) { // regular
                hp = Real(1.0);
                tmp0 = Real(-1.0);
                tmp1 = x(i+1,j,k);
            } else {
                hp = Real(1.0) + Real(2.) * ecx(i,j,k);
                tmp0 = Real(-1.0) / hp;
                tmp1 = Real(0.0);
            }

            if (ecx(i-1,j,k) == Real(1.0)) {
                hm = Real(1.0);
                tmp0 += Real(-1.0);
                tmp1 += x(i-1,j,k);
            } else {
                hm = Real(1.0) - Real(2.) * ecx(i-1,j,k);
                tmp0 += Real(-1.0) / hm;
            }

            Real gamma = tmp0 * (bx * Real(2.0) / (hp+hm));
            Real rho   = tmp1 * (bx * Real(2.0) / (hp+hm));

            if (ecy(i,j,k) == Real(1.0)) {
                hp = Real(1.0);
                tmp0 = Real(-1.0);
                tmp1 = x(i,j+1,k);
            } else {
                hp = Real(1.0) + Real(2.) * ecy(i,j,k);
                tmp0 = Real(-1.0) / hp;
                tmp1 = Real(0.0);
            }

            if (ecy(i,j-1,k) == Real(1.0)) {
                hm = Real(1.0);
                tmp0 += Real(-1.0);
                tmp1 += x(i,j-1,k);
            } else {
                hm = Real(1.0) - Real(2.) * ecy(i,j-1,k);
                tmp0 += Real(-1.0) / hm;
            }

            gamma += tmp0 * (by * Real(2.0) / (hp+hm));
            rho   += tmp1 * (by * Real(2.0) / (hp+hm));

            x(i,j,k) = (rhs(i,j,k) - rho) / gamma;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb(int i, int j, int k, Array4<Real> const& x,
                      Array4<Real const> const& rhs, Array4<int const> const& dmsk,
                      Array4<Real const> const& a, Real bx, Real by, int redblack) noexcept
{
    if ((i+j+k+redblack)%2 == 0) {
        if (dmsk(i,j,k)) {
            x(i,j,k) = Real(0.);
        } else {
            Real gamma = Real(-2.0)*(bx+by) - a(i,j,k);
            Real rho = bx * (x(i-1,j,k) + x(i+1,j,k))
                +      by * (x(i,j-1,k) + x(i,j+1,k));
            x(i,j,k) = (rhs(i,j,k) - rho) / gamma;
        }
    }
}

}

#endif
