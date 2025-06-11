#ifndef AMREX_ML_NODE_LINOP_3D_K_H_
#define AMREX_ML_NODE_LINOP_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
inline void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                             GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                             GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    Box gdomain = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (! bflo[idim]) { gdomain.growLo(idim,1); }
        if (! bfhi[idim]) { gdomain.growHi(idim,1); }
    }

    if (gdomain.strictly_contains(vbx)) { return; }

    const int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Box const& sbox = amrex::grow(vbx,1);
    AMREX_HOST_DEVICE_FOR_3D(sbox, i, j, k,
    {
        if (! gdomain.contains(IntVect(i,j,k))) {
            // xlo & ylo & zlo
            if (i == dlo.x-1 && j == dlo.y-1 && k == dlo.z-1 && (bflo[0] || bflo[1] || bflo[2]))
            {
                if (bflo[0] && bflo[1] && bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k+1+offset);
                }
                else if (bflo[0] && bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0] && bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bflo[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & ylo & zlo
            else if (i == dhi.x+1 && j == dlo.y-1 && k == dlo.z-1 && (bfhi[0] || bflo[1] || bflo[2]))
            {
                if (bfhi[0] && bflo[1] && bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k+1+offset);
                }
                else if (bfhi[0] && bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0] && bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bflo[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & yhi & zlo
            else if (i == dlo.x-1 && j == dhi.y+1 && k == dlo.z-1 && (bflo[0] || bfhi[1] || bflo[2]))
            {
                if (bflo[0] && bfhi[1] && bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k+1+offset);
                }
                else if (bflo[0] && bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0] && bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bfhi[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & yhi & zlo
            else if (i == dhi.x+1 && j == dhi.y+1 && k == dlo.z-1 && (bfhi[0] || bfhi[1] || bflo[2]))
            {
                if (bfhi[0] && bfhi[1] && bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k+1+offset);
                }
                else if (bfhi[0] && bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0] && bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bfhi[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & ylo & zhi
            else if (i == dlo.x-1 && j == dlo.y-1 && k == dhi.z+1 && (bflo[0] || bflo[1] || bfhi[2]))
            {
                if (bflo[0] && bflo[1] && bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k-1-offset);
                }
                else if (bflo[0] && bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0] && bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bflo[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & ylo & zhi
            else if (i == dhi.x+1 && j == dlo.y-1 && k == dhi.z+1 && (bfhi[0] || bflo[1] || bfhi[2]))
            {
                if (bfhi[0] && bflo[1] && bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k-1-offset);
                }
                else if (bfhi[0] && bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0] && bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bflo[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xlo & yhi & zhi
            else if (i == dlo.x-1 && j == dhi.y+1 && k == dhi.z+1 && (bflo[0] || bfhi[1] || bfhi[2]))
            {
                if (bflo[0] && bfhi[1] && bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k-1-offset);
                }
                else if (bflo[0] && bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0] && bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bfhi[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & yhi & zhi
            else if (i == dhi.x+1 && j == dhi.y+1 && k == dhi.z+1 && (bfhi[0] || bfhi[1] || bfhi[2]))
            {
                if (bfhi[0] && bfhi[1] && bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k-1-offset);
                }
                else if (bfhi[0] && bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0] && bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bfhi[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xlo & ylo
            else if (i == dlo.x-1 && j == dlo.y-1 && (bflo[0] || bflo[1]))
            {
                if (bflo[0] && bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xhi & ylo
            else if (i == dhi.x+1 && j == dlo.y-1 && (bfhi[0] || bflo[1]))
            {
                if (bfhi[0] && bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xlo & yhi
            else if (i == dlo.x-1 && j == dhi.y+1 && (bflo[0] || bfhi[1]))
            {
                if (bflo[0] && bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            // xhi & yhi
            else if (i == dhi.x+1 && j == dhi.y+1 && (bfhi[0] || bfhi[1]))
            {
                if (bfhi[0] && bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            // xlo & zlo
            else if (i == dlo.x-1 && k == dlo.z-1 && (bflo[0] || bflo[2]))
            {
                if (bflo[0] && bflo[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k+1+offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xhi & zlo
            else if (i == dhi.x+1 && k == dlo.z-1 && (bfhi[0] || bflo[2]))
            {
                if (bfhi[0] && bflo[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k+1+offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // xlo & zhi
            else if (i == dlo.x-1 && k == dhi.z+1 && (bflo[0] || bfhi[2]))
            {
                if (bflo[0] && bfhi[2])
                {
                    a(i,j,k) = a(i+1+offset, j, k-1-offset);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // xhi & zhi
            else if (i == dhi.x+1 && k == dhi.z+1 && (bfhi[0] || bfhi[2]))
            {
                if (bfhi[0] && bfhi[2])
                {
                    a(i,j,k) = a(i-1-offset, j, k-1-offset);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // ylo & zlo
            else if (j == dlo.y-1 && k == dlo.z-1 && (bflo[1] || bflo[2]))
            {
                if (bflo[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k+1+offset);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // yhi & zlo
            else if (j == dhi.y+1 && k == dlo.z-1 && (bfhi[1] || bflo[2]))
            {
                if (bfhi[1] && bflo[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k+1+offset);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bflo[2])
                {
                    a(i,j,k) = a(i, j, k+1+offset);
                }
            }
            // ylo & zhi
            else if (j == dlo.y-1 && k == dhi.z+1 && (bflo[1] || bfhi[2]))
            {
                if (bflo[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j+1+offset, k-1-offset);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            // yhi & zhi
            else if (j == dhi.y+1 && k == dhi.z+1 && (bfhi[1] || bfhi[2]))
            {
                if (bfhi[1] && bfhi[2])
                {
                    a(i,j,k) = a(i, j-1-offset, k-1-offset);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
                else if (bfhi[2])
                {
                    a(i,j,k) = a(i, j, k-1-offset);
                }
            }
            else if (i == dlo.x-1 && bflo[0])
            {
                a(i,j,k) = a(i+1+offset, j, k);
            }
            else if (i == dhi.x+1 && bfhi[0])
            {
                a(i,j,k) = a(i-1-offset, j, k);
            }
            else if (j == dlo.y-1 && bflo[1])
            {
                a(i,j,k) = a(i, j+1+offset, k);
            }
            else if (j == dhi.y+1 && bfhi[1])
            {
                a(i,j,k) = a(i, j-1-offset, k);
            }
            else if (k == dlo.z-1 && bflo[2])
            {
                a(i,j,k) = a(i, j, k+1+offset);
            }
            else if (k == dhi.z+1 && bfhi[2])
            {
                a(i,j,k) = a(i, j, k-1-offset);
            }
        }
    });
}

//
// restriction
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = Real(0.0);
    } else {
        crse(i,j,k) = Real(1./64.)*(fine(ii-1,jj-1,kk-1)+fine(ii+1,jj-1,kk-1)
                                   +fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1)
                                   +fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1)
                                   +fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1))
                    + Real(1./32.)*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1)
                                   +fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)
                                   +fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1)
                                   +fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)
                                   +fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  )
                                   +fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ))
                    + Real(1./16.)*(fine(ii-1,jj,kk)+fine(ii+1,jj,kk)
                                   +fine(ii,jj-1,kk)+fine(ii,jj+1,kk)
                                   +fine(ii,jj,kk-1)+fine(ii,jj,kk+1))
                      + Real(1./8.)*fine(ii,jj,kk);
    }
}

template <int rr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk,
                          Box const& fdom,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const int ii = i*rr;
    const int jj = j*rr;
    const int kk = k*rr;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = Real(0.0);
    } else {
        const auto ndlo = amrex::lbound(fdom);
        const auto ndhi = amrex::ubound(fdom);
        Real tmp = Real(0.0);
        for (int koff = -rr+1; koff <= rr-1; ++koff) {
            Real wz = rr - std::abs(koff);
            for (int joff = -rr+1; joff <= rr-1; ++joff) {
                Real wy = rr - std::abs(joff);
                for (int ioff = -rr+1; ioff <= rr-1; ++ioff) {
                    Real wx = rr - std::abs(ioff);
                    int itmp = ii + ioff;
                    int jtmp = jj + joff;
                    int ktmp = kk + koff;
                    if ((itmp < ndlo.x && (bclo[0] == LinOpBCType::Neumann ||
                                           bclo[0] == LinOpBCType::inflow)) ||
                        (itmp > ndhi.x && (bchi[0] == LinOpBCType::Neumann ||
                                           bchi[0] == LinOpBCType::inflow))) {
                        itmp = ii - ioff;
                    }
                    if ((jtmp < ndlo.y && (bclo[1] == LinOpBCType::Neumann ||
                                           bclo[1] == LinOpBCType::inflow)) ||
                        (jtmp > ndhi.y && (bchi[1] == LinOpBCType::Neumann ||
                                           bchi[1] == LinOpBCType::inflow))) {
                        jtmp = jj - joff;
                    }
                    if ((ktmp < ndlo.z && (bclo[2] == LinOpBCType::Neumann ||
                                           bclo[2] == LinOpBCType::inflow)) ||
                        (ktmp > ndhi.z && (bchi[2] == LinOpBCType::Neumann ||
                                           bchi[2] == LinOpBCType::inflow))) {
                        ktmp = kk - koff;
                    }
                    tmp += wx*wy*wz*fine(itmp,jtmp,ktmp);
                }
            }
        }
        crse(i,j,k) = tmp/Real(rr*rr*rr*rr*rr*rr);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int i, int j, int k, Array4<Real> const& crse,
                               Array4<Real const> const& fine, Array4<int const> const& msk, int idir) noexcept
{
    if (idir == 2)
    {
        int ii = i*2;
        int jj = j*2;
        int kk = k;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = Real(0.0);
        } else { // use 2-D formula
            crse(i,j,k) = Real(1./16.)*(fine(ii-1,jj-1,kk) + Real(2.)*fine(ii  ,jj-1,kk) +          fine(ii+1,jj-1,kk)
                             + Real(2.)*fine(ii-1,jj  ,kk) + Real(4.)*fine(ii  ,jj  ,kk) + Real(2.)*fine(ii+1,jj  ,kk)
                                      + fine(ii-1,jj+1,kk) + Real(2.)*fine(ii  ,jj+1,kk) +          fine(ii+1,jj+1,kk));
        }
    }
    else if (idir == 1)
    {
        int ii = i*2;
        int jj = j;
        int kk = k*2;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = Real(0.0);
        } else { // use 2-D formula
            crse(i,j,k) = Real(1./16.)*(fine(ii-1,jj,kk-1) + Real(2.)*fine(ii  ,jj,kk-1) +          fine(ii+1,jj,kk-1)
                             + Real(2.)*fine(ii-1,jj  ,kk) + Real(4.)*fine(ii  ,jj,kk  ) + Real(2.)*fine(ii+1,jj,kk  )
                                      + fine(ii-1,jj,kk+1) + Real(2.)*fine(ii  ,jj,kk+1) +          fine(ii+1,jj,kk+1));
        }
    }
    else if (idir == 0)
    {
        int ii = i;
        int jj = j*2;
        int kk = k*2;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = Real(0.0);
        } else { // use 2-D formula
            crse(i,j,k) = Real(1./16.)*(fine(ii,jj-1,kk-1) + Real(2.)*fine(ii  ,jj,kk-1) +          fine(ii,jj+1,kk-1)
                             + Real(2.)*fine(ii,jj-1  ,kk) + Real(4.)*fine(ii  ,jj,kk  ) + Real(2.)*fine(ii,jj+1,kk  )
                                      + fine(ii,jj-1,kk+1) + Real(2.)*fine(ii  ,jj,kk+1) +          fine(ii,jj+1,kk+1));
        }
    }
    else
    {
        amrex::Abort("mlndlap_semi_restriction semi direction wrong semi-direction. ");
    }
}

//
// masks
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int i, int j, int k, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    using namespace nodelap_detail;

    int s = cmsk(i-1,j-1,k-1) + cmsk(i  ,j-1,k-1)
        +   cmsk(i-1,j  ,k-1) + cmsk(i  ,j  ,k-1)
        +   cmsk(i-1,j-1,k  ) + cmsk(i  ,j-1,k  )
        +   cmsk(i-1,j  ,k  ) + cmsk(i  ,j  ,k  );
    if (s == 8*crse_cell) {
        nmsk(i,j,k) = crse_node;
    }
    else if (s == 8*fine_cell) {
        nmsk(i,j,k) = fine_node;
    } else {
        nmsk(i,j,k) = crse_fine_node;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if (!dmsk(i,j,k)) {
            dmsk(i,j,k) = (omsk(i-1,j-1,k-1) == 1 || omsk(i,j-1,k-1) == 1 ||
                           omsk(i-1,j  ,k-1) == 1 || omsk(i,j  ,k-1) == 1 ||
                           omsk(i-1,j-1,k  ) == 1 || omsk(i,j-1,k  ) == 1 ||
                           omsk(i-1,j  ,k  ) == 1 || omsk(i,j  ,k  ) == 1);
        }
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if (bclo[0] == LinOpBCType::Dirichlet && lo.x == domlo.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) = 1;
        }}
    }

    if (bchi[0] == LinOpBCType::Dirichlet && hi.x == domhi.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) = 1;
        }}
    }

    if (bclo[1] == LinOpBCType::Dirichlet && lo.y == domlo.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) = 1;
        }}
    }

    if (bchi[1] == LinOpBCType::Dirichlet && hi.y == domhi.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) = 1;
        }}
    }

    if (bclo[2] == LinOpBCType::Dirichlet && lo.z == domlo.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) = 1;
        }}
    }

    if (bchi[2] == LinOpBCType::Dirichlet && hi.z == domhi.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) = 1;
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,k) = static_cast<Real>(omsk(i,j,k));
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
        && lo.x == domlo.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) *= Real(0.5);
        }}
    }

    if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
        && hi.x == domhi.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) *= Real(0.5);
        }}
    }

    if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
        && lo.y == domlo.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) *= Real(0.5);
        }}
    }

    if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
        && hi.y == domhi.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) *= Real(0.5);
        }}
    }

    if ((bclo[2] == LinOpBCType::Neumann || bclo[2] == LinOpBCType::inflow)
        && lo.z == domlo.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) *= Real(0.5);
        }}
    }

    if ((bchi[2] == LinOpBCType::Neumann || bchi[2] == LinOpBCType::inflow)
        && hi.z == domhi.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) *= Real(0.5);
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk,
                           Array4<int const> const& fmsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if (fmsk(i,j,k) == 0) {
            dmsk(i,j,k) = static_cast<Real>(omsk(i,j,k));
        } else {
            dmsk(i,j,k) = Real(0);
        }
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
        && lo.x == domlo.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) *= Real(0.5);
        }}
    }

    if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
        && hi.x == domhi.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) *= Real(0.5);
        }}
    }

    if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
        && lo.y == domlo.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) *= Real(0.5);
        }}
    }

    if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
        && hi.y == domhi.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) *= Real(0.5);
        }}
    }

    if ((bclo[2] == LinOpBCType::Neumann || bclo[2] == LinOpBCType::inflow)
        && lo.z == domlo.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) *= Real(0.5);
        }}
    }

    if ((bchi[2] == LinOpBCType::Neumann || bchi[2] == LinOpBCType::inflow)
        && hi.z == domhi.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) *= Real(0.5);
        }}
    }
}

}

#endif
