#ifndef AMREX_MLNODELAP_3D_K_H_
#define AMREX_MLNODELAP_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int i, int j, int k, Array4<T> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,k-1) == fine_flag &&
        msk(i  ,j-1,k-1) == fine_flag &&
        msk(i-1,j  ,k-1) == fine_flag &&
        msk(i  ,j  ,k-1) == fine_flag &&
        msk(i-1,j-1,k  ) == fine_flag &&
        msk(i  ,j-1,k  ) == fine_flag &&
        msk(i-1,j  ,k  ) == fine_flag &&
        msk(i  ,j  ,k  ) == fine_flag)
    {
        phi(i,j,k) = T(0.0);
    }
}

//
// coeffs
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T cl = fine(2*i  ,2*j,2*k  )+fine(2*i  ,2*j+1,2*k  )+
              fine(2*i  ,2*j,2*k+1)+fine(2*i  ,2*j+1,2*k+1);
    T cr = fine(2*i+1,2*j,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i+1,2*j,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = T(0.5)*cl*cr/(cl+cr);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_y (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1);
    T cr = fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = T(0.5)*cl*cr/(cl+cr);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_z (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine) noexcept
{
    T cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  );
    T cr = fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1)+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = T(0.5)*cl*cr/(cl+cr);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine, int idir) noexcept
{
    if (idir == 2) {
        T cl = fine(2*i  ,2*j,k) + fine(2*i  ,2*j+1,k);
        T cr = fine(2*i+1,2*j,k) + fine(2*i+1,2*j+1,k);
        crse(i,j,k) = cl*cr/(cl+cr);
    } else if (idir == 1) {
        T cl = fine(2*i  ,j,2*k) + fine(2*i  ,j,2*k+1);
        T cr = fine(2*i+1,j,2*k) + fine(2*i+1,j,2*k+1);
        crse(i,j,k) = cl*cr/(cl+cr);
    } else {
        T cl = fine(i,2*j  ,2*k) + fine(i,2*j  ,2*k+1);
        T cr = fine(i,2*j+1,2*k) + fine(i,2*j+1,2*k+1);
        crse(i,j,k) = cl*cr/(cl+cr);
    }
}

//
// operator
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_ha (int i, int j, int k, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                       Array4<std::type_identity_t<T> const> const& sz, Array4<int const> const& msk,
                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1./36.)*dxinv[0]*dxinv[0];
        T facy = T(1./36.)*dxinv[1]*dxinv[1];
        T facz = T(1./36.)*dxinv[2]*dxinv[2];
        T y   = x(i,j,k)*T(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                             +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                       +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                             +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                       +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                             +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        y        += x(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                   +facy*sy(i-1,j-1,k-1)
                                   +facz*sz(i-1,j-1,k-1))
                  + x(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                   +facy*sy(i  ,j-1,k-1)
                                   +facz*sz(i  ,j-1,k-1))
                  + x(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                   +facy*sy(i-1,j  ,k-1)
                                   +facz*sz(i-1,j  ,k-1))
                  + x(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                   +facy*sy(i  ,j  ,k-1)
                                   +facz*sz(i  ,j  ,k-1))
                  + x(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                   +facy*sy(i-1,j-1,k  )
                                   +facz*sz(i-1,j-1,k  ))
                  + x(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                   +facy*sy(i  ,j-1,k  )
                                   +facz*sz(i  ,j-1,k  ))
                  + x(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                   +facy*sy(i-1,j  ,k  )
                                   +facz*sz(i-1,j  ,k  ))
                  + x(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                   +facy*sy(i  ,j  ,k  )
                                   +facz*sz(i  ,j  ,k  ));
        y        += x(i  ,j-1,k-1)*(          -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                    +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                    +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                  + x(i  ,j+1,k-1)*(          -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                    +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                    +T(2.0)*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                  + x(i  ,j-1,k+1)*(          -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                    +T(2.0)*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                    +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                  + x(i  ,j+1,k+1)*(          -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                    +T(2.0)*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                    +T(2.0)*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                  + x(i-1,j  ,k-1)*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                              -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                    +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                  + x(i+1,j  ,k-1)*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                              -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                    +T(2.0)*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                  + x(i-1,j  ,k+1)*( T(2.0)*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                              -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                    +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                  + x(i+1,j  ,k+1)*( T(2.0)*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                              -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                    +T(2.0)*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                  + x(i-1,j-1,k  )*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                    +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                              -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                  + x(i+1,j-1,k  )*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                    +T(2.0)*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                              -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                  + x(i-1,j+1,k  )*( T(2.0)*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                    +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                              -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                  + x(i+1,j+1,k  )*( T(2.0)*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                    +T(2.0)*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                              -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)));
        y            += T(2.0)*x(i-1,j,k)*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                                        -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                                        -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                      + T(2.0)*x(i+1,j,k)*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                                        -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                                        -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                      + T(2.0)*x(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                        +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                                  -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                      + T(2.0)*x(i,j+1,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                        +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                                  -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                      + T(2.0)*x(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                                  -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                        +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                      + T(2.0)*x(i,j,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                                  -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        return y;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_aa (int i, int j, int k, Array4<T const> const& x,
                       Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& msk,
                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
        T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
        T fxyz = facx + facy + facz;
        T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
        T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
        T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
        T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
        T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
        T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;
        return x(i,j,k)*T(-4.0)*fxyz*
            (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
            + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + x(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + x(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + x(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + x(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_c (int i, int j, int k, Array4<T const> const& x,
                      std::type_identity_t<T> sig, Array4<int const> const& msk,
                      GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
        T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
        T fxyz = facx + facy + facz;
        T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
        T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
        T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
        T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
        T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
        T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;
        return sig * (x(i,j,k)*T(-4.0)*fxyz*T(8.)
            + fxyz*(x(i-1,j-1,k-1)
                  + x(i+1,j-1,k-1)
                  + x(i-1,j+1,k-1)
                  + x(i+1,j+1,k-1)
                  + x(i-1,j-1,k+1)
                  + x(i+1,j-1,k+1)
                  + x(i-1,j+1,k+1)
                  + x(i+1,j+1,k+1))
            + fmx2y2z*(x(i  ,j-1,k-1)*T(2.)
                     + x(i  ,j+1,k-1)*T(2.)
                     + x(i  ,j-1,k+1)*T(2.)
                     + x(i  ,j+1,k+1)*T(2.))
            + f2xmy2z*(x(i-1,j  ,k-1)*T(2.)
                     + x(i+1,j  ,k-1)*T(2.)
                     + x(i-1,j  ,k+1)*T(2.)
                     + x(i+1,j  ,k+1)*T(2.))
            + f2x2ymz*(x(i-1,j-1,k  )*T(2.)
                     + x(i+1,j-1,k  )*T(2.)
                     + x(i-1,j+1,k  )*T(2.)
                     + x(i+1,j+1,k  )*T(2.))
            + f4xm2ym2z*(x(i-1,j,k)*T(4.)
                       + x(i+1,j,k)*T(4.))
            + fm2x4ym2z*(x(i,j-1,k)*T(4.)
                       + x(i,j+1,k)*T(4.))
            + fm2xm2y4z*(x(i,j,k-1)*T(4.)
                       + x(i,j,k+1)*T(4.)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (int i, int j, int k, Array4<T> const& x, Array4<std::type_identity_t<T> const> const& sx,
                           Array4<std::type_identity_t<T> const> const& sy, Array4<std::type_identity_t<T> const> const& sz,
                           Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];

    if (!msk(i,j,k)) {
        x(i,j,k) = x(i,j,k)/(T(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                                +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                         +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                                +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                         +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                                +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  ))));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (int i, int j, int k, Array4<T> const& x, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
    T fxyz = facx + facy + facz;

    if (!msk(i,j,k)) {
        x(i,j,k) = x(i,j,k) /
           (T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                             +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy, Array4<std::type_identity_t<T> const> const& sz,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(-4.0 / 36.0)*dxinv[0]*dxinv[0];
    T facy = T(-4.0 / 36.0)*dxinv[1]*dxinv[1];
    T facz = T(-4.0 / 36.0)*dxinv[2]*dxinv[2];

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                    +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
              +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                    +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
              +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                    +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
    }
}

template <typename T>
inline
void mlndlap_jacobi_ha (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy, Array4<std::type_identity_t<T> const> const& sz,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(-4.0 / 36.0)*dxinv[0]*dxinv[0];
    T facy = T(-4.0 / 36.0)*dxinv[1]*dxinv[1];
    T facz = T(-4.0 / 36.0)*dxinv[2]*dxinv[2];

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                        +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                  +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                        +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                  +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                        +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        }
    });
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fxyz = T(-4.0 / 36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                    +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_c (int i, int j, int k, Array4<T> const& sol, std::type_identity_t<T> Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fxyz = T(-4.0 / 36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (fxyz*T(8.)*sig);
    }
}

template <typename T>
inline
void mlndlap_jacobi_aa (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                        Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                        Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fxyz = T(-4.0 / 36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                        +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

template <typename T>
inline
void mlndlap_jacobi_c (Box const& bx, Array4<T> const& sol, Array4<std::type_identity_t<T> const> const& Ax,
                       Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                       Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T fxyz = T(-4.0 / 36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            sol(i,j,k) += T(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fxyz*T(8.)*sig);
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                              Array4<std::type_identity_t<T> const> const& sy, Array4<std::type_identity_t<T> const> const& sz,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = T(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                       +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                 +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                       +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                 +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                       +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
            T Ax = sol(i,j,k)*s0
                     + sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                        +facy*sy(i-1,j-1,k-1)
                                        +facz*sz(i-1,j-1,k-1))
                     + sol(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                        +facy*sy(i  ,j-1,k-1)
                                        +facz*sz(i  ,j-1,k-1))
                     + sol(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                        +facy*sy(i-1,j  ,k-1)
                                        +facz*sz(i-1,j  ,k-1))
                     + sol(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                        +facy*sy(i  ,j  ,k-1)
                                        +facz*sz(i  ,j  ,k-1))
                     + sol(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                        +facy*sy(i-1,j-1,k  )
                                        +facz*sz(i-1,j-1,k  ))
                     + sol(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                        +facy*sy(i  ,j-1,k  )
                                        +facz*sz(i  ,j-1,k  ))
                     + sol(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                        +facy*sy(i-1,j  ,k  )
                                        +facz*sz(i-1,j  ,k  ))
                     + sol(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                        +facy*sy(i  ,j  ,k  )
                                        +facz*sz(i  ,j  ,k  ))
                     +sol(i  ,j-1,k-1)*(          -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                        +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                        +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                     +sol(i  ,j+1,k-1)*(          -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                        +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                        +T(2.0)*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                     +sol(i  ,j-1,k+1)*(          -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                        +T(2.0)*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                        +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                     +sol(i  ,j+1,k+1)*(          -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                        +T(2.0)*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                        +T(2.0)*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                     +sol(i-1,j  ,k-1)*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                                  -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                        +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                     +sol(i+1,j  ,k-1)*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                                  -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                        +T(2.0)*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                     +sol(i-1,j  ,k+1)*( T(2.0)*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                                  -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                        +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                     +sol(i+1,j  ,k+1)*( T(2.0)*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                                  -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                        +T(2.0)*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                     +sol(i-1,j-1,k  )*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                        +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                                  -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                     +sol(i+1,j-1,k  )*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                        +T(2.0)*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                                  -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                     +sol(i-1,j+1,k  )*( T(2.0)*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                        +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                                  -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                     +sol(i+1,j+1,k  )*( T(2.0)*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                        +T(2.0)*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                                  -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                     + T(2.0)*sol(i-1,j,k)*(T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                                        -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                                        -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                     + T(2.0)*sol(i+1,j,k)*(T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                                        -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                                        -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                     + T(2.0)*sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                        +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                                  -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                     + T(2.0)*sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                        +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                                  -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                     + T(2.0)*sol(i,j,k-1)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                                  -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                        +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                     + T(2.0)*sol(i,j,k+1)*(   -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                                  -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));

                sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
    T fxyz = facx + facy + facz;
    T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
    T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
    T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
    T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
    T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
    T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                      +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
            T Ax = sol(i,j,k)*s0
                + fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                         + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                         + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                         + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                         + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                         + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                         + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                         + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                         + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                         + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                           + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                           + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                           + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

template <typename T>
inline
void mlndlap_gauss_seidel_c (Box const& bx, Array4<T> const& sol,
                             Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                             Array4<int const> const& msk,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
    T fxyz = facx + facy + facz;
    T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
    T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
    T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
    T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
    T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
    T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else {
            T s0 = T(-4.0)*fxyz*T(8.);
            T Ax = sol(i,j,k)*s0
                + fxyz*(sol(i-1,j-1,k-1)
                      + sol(i+1,j-1,k-1)
                      + sol(i-1,j+1,k-1)
                      + sol(i+1,j+1,k-1)
                      + sol(i-1,j-1,k+1)
                      + sol(i+1,j-1,k+1)
                      + sol(i-1,j+1,k+1)
                      + sol(i+1,j+1,k+1))
                + fmx2y2z*(sol(i  ,j-1,k-1)*T(2.)
                         + sol(i  ,j+1,k-1)*T(2.)
                         + sol(i  ,j-1,k+1)*T(2.)
                         + sol(i  ,j+1,k+1)*T(2.))
                + f2xmy2z*(sol(i-1,j  ,k-1)*T(2.)
                         + sol(i+1,j  ,k-1)*T(2.)
                         + sol(i-1,j  ,k+1)*T(2.)
                         + sol(i+1,j  ,k+1)*T(2.))
                + f2x2ymz*(sol(i-1,j-1,k  )*T(2.)
                         + sol(i+1,j-1,k  )*T(2.)
                         + sol(i-1,j+1,k  )*T(2.)
                         + sol(i+1,j+1,k  )*T(2.))
                + f4xm2ym2z*(sol(i-1,j,k)*T(4.)
                           + sol(i+1,j,k)*T(4.))
                + fm2x4ym2z*(sol(i,j-1,k)*T(4.)
                           + sol(i,j+1,k)*T(4.))
                + fm2xm2y4z*(sol(i,j,k-1)*T(4.)
                           + sol(i,j,k+1)*T(4.));

            sol(i,j,k) += (rhs(i,j,k) - Ax*sig) / (s0*sig);
        }
    });
}

template <typename T>
AMREX_FORCE_INLINE
void tridiagonal_solve (Array1D<T,0,31>& a_ls, Array1D<T,0,31>& b_ls, Array1D<T,0,31>& c_ls,
                        Array1D<T,0,31>& r_ls, Array1D<T,0,31>& u_ls, Array1D<T,0,31>& gam,
                        int ilen )
{
    T bet = b_ls(0);
    u_ls(0) = r_ls(0) / bet;

    for (int i = 1; i <= ilen - 1; i++) {
        gam(i) = c_ls(i-1) / bet;
        bet = b_ls(i) - a_ls(i)*gam(i);
        if (bet == 0) { amrex::Abort(">>>TRIDIAG FAILED"); }
        u_ls(i) = (r_ls(i)-a_ls(i)*u_ls(i-1)) / bet;
    }
    for (int i = ilen-2; i >= 0; i--) {
        u_ls(i) = u_ls(i) - gam(i+1)*u_ls(i+1);
    }
}

template <typename T>
inline
void mlndlap_gauss_seidel_with_line_solve_aa (Box const& bx, Array4<T> const& sol,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv)
{
    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
    T fxyz = facx + facy + facz;
    T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
    T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
    T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
    T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
    T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
    T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    int idir = -1;
    int ilen = 33;

    if ( (dxinv[0] <= dxinv[2]) && (dxinv[1] <= dxinv[2]) ) {
        idir = 2;
        ilen = hi.z - lo.z + 1;
    }
    if ( (dxinv[0] <= dxinv[1]) && (dxinv[2] <= dxinv[1]) ) {
        idir = 1;
        ilen = hi.y - lo.y + 1;
    }
    if ( (dxinv[1] <= dxinv[0]) && (dxinv[2] <= dxinv[0]) ) {
        idir = 0;
        ilen = hi.x - lo.x + 1;
    }

    if (ilen > 32) {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is hard-wired to be no longer than 32");
    }

    Array1D<T,0,31> a_ls,b_ls,c_ls,u_ls,r_ls,gam;


    if ( idir == 2 )
    {
        for (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                for (int k = lo.z; k <= hi.z; ++k)
                {
                    if (msk(i,j,k))
                    {
                        a_ls(k-lo.z) = T(0.);
                        b_ls(k-lo.z) = T(1.);
                        c_ls(k-lo.z) = T(0.);
                        u_ls(k-lo.z) = T(0.);
                        r_ls(k-lo.z) = T(0.);
                    }
                    else
                    {
                        T s0 = T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                                 + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        T Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                            + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                                     + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                                     + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                                     + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                            + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                                     + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                                     + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                                     + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                            + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                                     + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                                     + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                                     + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                            + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                                       + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                            + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                                       + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)));

                        a_ls(k-lo.z) = fm2xm2y4z*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1));
                        b_ls(k-lo.z) = s0;
                        c_ls(k-lo.z) = fm2xm2y4z*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
                        u_ls(k-lo.z) = T(0.);
                        r_ls(k-lo.z) = rhs(i,j,k) - Ax;
                    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int k = lo.z; k <= hi.z; ++k)
                {
                    sol(i,j,k) = u_ls(k-lo.z);
                }
            }
        }
    }
    else if (idir == 1)
    {
        for (int k = lo.z; k <= hi.z; ++k)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                for (int j = lo.y; j <= hi.y; ++j)
                {
                    if (msk(i,j,k)) {
                        a_ls(j-lo.y) = T(0.);
                        b_ls(j-lo.y) = T(1.);
                        c_ls(j-lo.y) = T(0.);
                        u_ls(j-lo.y) = T(0.);
                        r_ls(j-lo.y) = T(0.);
                    }
                    else
                    {
                        T s0 = T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                                 + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        T Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                            + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                                     + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                                     + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                                     + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                            + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                                     + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                                     + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                                     + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                            + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                                     + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                                     + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                                     + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                            + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                                       + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                            + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                                       + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

                        a_ls(j-lo.y) = fm2x4ym2z*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k));
                        b_ls(j-lo.y) = s0;
                        c_ls(j-lo.y) = fm2x4ym2z*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k));
                        u_ls(j-lo.y) = T(0.);
                        r_ls(j-lo.y) = rhs(i,j,k) - Ax;

                    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int j = lo.y; j <= hi.y; ++j)
                {
                    sol(i,j,k) = u_ls(j-lo.y);
                }
            }
        }
    }
    else if (idir == 0)
    {
        for (int j = lo.y; j <= hi.y; ++j)
        {
            for (int k = lo.z; k <= hi.z; ++k)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if (msk(i,j,k))
                    {
                        a_ls(i-lo.x) = T(0.);
                        b_ls(i-lo.x) = T(1.);
                        c_ls(i-lo.x) = T(0.);
                        u_ls(i-lo.x) = T(0.);
                        r_ls(i-lo.x) = T(0.);
                    }
                    else
                    {
                        T s0 = T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                                 + sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));

                        T Ax = fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                            + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                                     + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                                     + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                                     + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                            + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                                     + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                                     + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                                     + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                            + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                                     + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                                     + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                                     + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                            + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                                     + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                            + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                                     + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

                        a_ls(i-lo.x) = f4xm2ym2z*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k));
                        b_ls(i-lo.x) = s0;
                        c_ls(i-lo.x) = f4xm2ym2z*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k));
                        u_ls(i-lo.x) = T(0.);
                        r_ls(i-lo.x) = rhs(i,j,k) - Ax;
                    }
                }
                tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

                for (int i = lo.x; i <= hi.x; ++i)
                {
                    sol(i,j,k) = u_ls(i-lo.x);
                }
            }
        }
    }
    else
    {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is wrong direction.");
    }
}

//
// interpolation
//

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_line_x (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        T w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic+1,jc,kc))/(w1+w2);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_line_y (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        T w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc+1,kc))/(w1+w2);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_line_z (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        T w2 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc,kc+1))/(w1+w2);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_face_xy (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        T w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        T w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        T w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_face_xz (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        T w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        T w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        T w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T aa_interp_face_yz (Array4<T const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        T w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        T w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        T w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sig,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sig,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_c (int i, int j, int k, Array4<T> const& fine,
                          Array4<std::type_identity_t<T> const> const& crse,
                          Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd && j_is_odd && k_is_odd) {
            // Fine node at center of cell
            fine(i,j,k) += T(0.125) *
                (crse(ic,jc  ,kc  ) + crse(ic+1,jc  ,kc  ) +
                 crse(ic,jc+1,kc  ) + crse(ic+1,jc+1,kc  ) +
                 crse(ic,jc  ,kc+1) + crse(ic+1,jc  ,kc+1) +
                 crse(ic,jc+1,kc+1) + crse(ic+1,jc+1,kc+1));
        } else if (j_is_odd && k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += T(0.25) *
                (crse(ic,jc,kc  ) + crse(ic,jc+1,kc  ) +
                 crse(ic,jc,kc+1) + crse(ic,jc+1,kc+1));
        } else if (i_is_odd && k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += T(0.25) *
                (crse(ic,jc,kc  ) + crse(ic+1,jc,kc  ) +
                 crse(ic,jc,kc+1) + crse(ic+1,jc,kc+1));
        } else if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += T(0.25) *
                (crse(ic,jc  ,kc) + crse(ic+1,jc  ,kc) +
                 crse(ic,jc+1,kc) + crse(ic+1,jc+1,kc));
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += T(0.5)*(crse(ic,jc,kc)+crse(ic+1,jc,kc));
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += T(0.5)*(crse(ic,jc,kc)+crse(ic,jc+1,kc));
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += T(0.5)*(crse(ic,jc,kc)+crse(ic,jc,kc+1));
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int j, int k, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd && j_is_odd && k_is_odd) {
            // Fine node at center of cell
            T w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
            T w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
            T w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
            T w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
            T w5 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
            T w6 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
            fine(i,j,k) += (w1 * aa_interp_face_yz(crse,sig,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * aa_interp_face_yz(crse,sig,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * aa_interp_face_xz(crse,sig,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * aa_interp_face_xz(crse,sig,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * aa_interp_face_xy(crse,sig,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * aa_interp_face_xy(crse,sig,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd && k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd && k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int i, int j, int k, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sig,
                           Array4<int const> const& msk, int idir)
{
    if (idir == 2 )
    {
        if (!msk(i,j,k)) {
            int ic = amrex::coarsen(i,2);
            int jc = amrex::coarsen(j,2);
            int kc = k;
            bool i_is_odd = (ic*2 != i);
            bool j_is_odd = (jc*2 != j);

            if (i_is_odd && j_is_odd) {
                // Node on a X-Y face
                fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,kc);
            }
        }
    } else if (idir ==1 ){
        if (!msk(i,j,k)) {
            int ic = amrex::coarsen(i,2);
            int jc = j;
            int kc = amrex::coarsen(k,2);
            bool i_is_odd = (ic*2 != i);
            bool k_is_odd = (kc*2 != k);

            if (i_is_odd && k_is_odd) {
                // Node on a X-Z face
                fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
            } else if (i_is_odd) {
                // Node on X line
                fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,kc);
            }
        }
    } else if (idir == 0 ) {
        if (!msk(i,j,k)) {
            int ic = i;
            int jc = amrex::coarsen(j,2);
            int kc = amrex::coarsen(k,2);
            bool j_is_odd = (jc*2 != j);
            bool k_is_odd = (kc*2 != k);

            if (j_is_odd && k_is_odd) {
                // Node on a Y-Z face
                fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
            } else if (j_is_odd) {
                // Node on Y line
                fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
            } else if (k_is_odd) {
                // Node on Z line
                fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
            } else {
                // Node coincident with coarse node
                fine(i,j,k) += crse(ic,jc,kc);
            }
        }
    } else {
        amrex::Abort("mlndlap_semi_interpolation semi direction wrong semi-direction. ");
    }
}


    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T ha_interp_face_xy (Array4<T const> const& crse,
                            Array4<std::type_identity_t<T> const> const& sigx, Array4<std::type_identity_t<T> const> const& sigy,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        T w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        T w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        T w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T ha_interp_face_xz (Array4<T const> const& crse,
                            Array4<std::type_identity_t<T> const> const& sigx, Array4<std::type_identity_t<T> const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        T w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        T w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        T w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T ha_interp_face_yz (Array4<T const> const& crse,
                            Array4<std::type_identity_t<T> const> const& sigy, Array4<std::type_identity_t<T> const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        T w1 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        T w2 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        T w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        T w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sigy,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sigy,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int i, int j, int k, Array4<T> const& fine,
                           Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sigx,
                           Array4<std::type_identity_t<T> const> const& sigy, Array4<std::type_identity_t<T> const> const& sigz,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd && j_is_odd && k_is_odd) {
            // Fine node at center of cell
            T w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
            T w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
            T w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
            T w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
            T w5 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
            T w6 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
            fine(i,j,k) += (w1 * ha_interp_face_yz(crse,sigy,sigz,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * ha_interp_face_yz(crse,sigy,sigz,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * ha_interp_face_xz(crse,sigx,sigz,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * ha_interp_face_xz(crse,sigx,sigz,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd && k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += ha_interp_face_yz(crse,sigy,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd && k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += ha_interp_face_xz(crse,sigx,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += ha_interp_face_xy(crse,sigx,sigy,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sigx,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sigy,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sigz,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

//
// rhs & u
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                   Box const& nodal_domain,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    T facx = T(0.25)*dxinv[0];
    T facy = T(0.25)*dxinv[1];
    T facz = T(0.25)*dxinv[2];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (msk(i,j,k)) {
        rhs(i,j,k) = T(0.0);
    } else {

        T zero_ilo = T(1.0);
        T zero_ihi = T(1.0);
        T zero_jlo = T(1.0);
        T zero_jhi = T(1.0);
        T zero_klo = T(1.0);
        T zero_khi = T(1.0);

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
            && i == domlo.x)
        {
            zero_ilo = T(0.0);
        }
        if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
            && i == domhi.x)
        {
            zero_ihi = T(0.0);
        }
        if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
            && j == domlo.y)
        {
            zero_jlo = T(0.0);
        }
        if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
            && j == domhi.y)
        {
            zero_jhi = T(0.0);
        }
        if ((bclo[2] == LinOpBCType::Neumann || bclo[2] == LinOpBCType::inflow)
            && k == domlo.z)
        {
            zero_klo = T(0.0);
        }
        if ((bchi[2] == LinOpBCType::Neumann || bchi[2] == LinOpBCType::inflow)
            && k == domhi.z)
        {
            zero_khi = T(0.0);
        }

        rhs(i,j,k) = facx*(-vel(i-1,j-1,k-1,0)*zero_jlo*zero_klo+vel(i,j-1,k-1,0)*zero_jlo*zero_klo
                           -vel(i-1,j  ,k-1,0)*zero_jhi*zero_klo+vel(i,j  ,k-1,0)*zero_jhi*zero_klo
                           -vel(i-1,j-1,k  ,0)*zero_jlo*zero_khi+vel(i,j-1,k  ,0)*zero_jlo*zero_khi
                           -vel(i-1,j  ,k  ,0)*zero_jhi*zero_khi+vel(i,j  ,k  ,0)*zero_jhi*zero_khi)

                   + facy*(-vel(i-1,j-1,k-1,1)*zero_ilo*zero_klo-vel(i,j-1,k-1,1)*zero_ihi*zero_klo
                           +vel(i-1,j  ,k-1,1)*zero_ilo*zero_klo+vel(i,j  ,k-1,1)*zero_ihi*zero_klo
                           -vel(i-1,j-1,k  ,1)*zero_ilo*zero_khi-vel(i,j-1,k  ,1)*zero_ihi*zero_khi
                           +vel(i-1,j  ,k  ,1)*zero_ilo*zero_khi+vel(i,j  ,k  ,1)*zero_ihi*zero_khi)

                   + facz*(-vel(i-1,j-1,k-1,2)*zero_ilo*zero_jlo-vel(i,j-1,k-1,2)*zero_ihi*zero_jlo
                           -vel(i-1,j  ,k-1,2)*zero_ilo*zero_jhi-vel(i,j  ,k-1,2)*zero_ihi*zero_jhi
                           +vel(i-1,j-1,k  ,2)*zero_ilo*zero_jlo+vel(i,j-1,k  ,2)*zero_ihi*zero_jlo
                           +vel(i-1,j  ,k  ,2)*zero_ilo*zero_jhi+vel(i,j  ,k  ,2)*zero_ihi*zero_jhi);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_rhcc (int i, int j, int k, Array4<T const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    T r;
    if (msk(i,j,k)) {
        r = T(0.0);
    } else {
        r = T(0.125) * (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) +
                           rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ));
    }
    return r;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                     Array4<std::type_identity_t<T> const> const& sig, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(0.25)*dxinv[0];
    T facy = T(0.25)*dxinv[1];
    T facz = T(0.25)*dxinv[2];
    u(i,j,k,0) -= sig(i,j,k)*facx
        * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,1) -= sig(i,j,k)*facy
        * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,2) -= sig(i,j,k)*facz
        * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
           +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
}

// Harmonic-average (anisotropic sigma) variant.  Used when the node-
// centered Poisson operator was set up with a multi-component sigma
// (see MLNodeLaplacian::setSigma and setMapped) so each velocity
// component is corrected with its own sigma direction.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_ha (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                        Array4<std::type_identity_t<T> const> const& sx,
                        Array4<std::type_identity_t<T> const> const& sy,
                        Array4<std::type_identity_t<T> const> const& sz,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(0.25)*dxinv[0];
    T facy = T(0.25)*dxinv[1];
    T facz = T(0.25)*dxinv[2];
    u(i,j,k,0) -= sx(i,j,k)*facx
        * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,1) -= sy(i,j,k)*facy
        * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,2) -= sz(i,j,k)*facz
        * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
           +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_c (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                       std::type_identity_t<T> sig, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T facx = T(0.25)*dxinv[0];
    T facy = T(0.25)*dxinv[1];
    T facz = T(0.25)*dxinv[2];
    u(i,j,k,0) -= sig*facx
        * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,1) -= sig*facy
        * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,2) -= sig*facz
        * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
           +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Df (int ii, int jj, int kk, std::type_identity_t<T> facx, std::type_identity_t<T> facy, std::type_identity_t<T> facz,
                         Array4<T const> const& vel, Box const& velbx) noexcept
    {
        T Df = T(0.0);
        if (velbx.contains(ii-1,jj-1,kk-1)) {
            Df += -facx*vel(ii-1,jj-1,kk-1,0) - facy*vel(ii-1,jj-1,kk-1,1) - facz*vel(ii-1,jj-1,kk-1,2);
        }
        if (velbx.contains(ii,jj-1,kk-1)) {
            Df += facx*vel(ii,jj-1,kk-1,0) - facy*vel(ii,jj-1,kk-1,1) - facz*vel(ii,jj-1,kk-1,2);
        }
        if (velbx.contains(ii-1,jj,kk-1)) {
            Df += -facx*vel(ii-1,jj,kk-1,0) + facy*vel(ii-1,jj,kk-1,1) - facz*vel(ii-1,jj,kk-1,2);
        }
        if (velbx.contains(ii,jj,kk-1)) {
            Df += facx*vel(ii,jj,kk-1,0) + facy*vel(ii,jj,kk-1,1) - facz*vel(ii,jj,kk-1,2);
        }
        if (velbx.contains(ii-1,jj-1,kk)) {
            Df += -facx*vel(ii-1,jj-1,kk,0) - facy*vel(ii-1,jj-1,kk,1) + facz*vel(ii-1,jj-1,kk,2);
        }
        if (velbx.contains(ii,jj-1,kk)) {
            Df += facx*vel(ii,jj-1,kk,0) - facy*vel(ii,jj-1,kk,1) + facz*vel(ii,jj-1,kk,2);
        }
        if (velbx.contains(ii-1,jj,kk)) {
            Df += -facx*vel(ii-1,jj,kk,0) + facy*vel(ii-1,jj,kk,1) + facz*vel(ii-1,jj,kk,2);
        }
        if (velbx.contains(ii,jj,kk)) {
            Df += facx*vel(ii,jj,kk,0) + facy*vel(ii,jj,kk,1) + facz*vel(ii,jj,kk,2);
        }
        return Df;
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_fine_contrib (int i, int j, int k, Box const& fvbx, Box const& velbx,
                                Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                                Array4<std::type_identity_t<T> const> const& frhs, Array4<int const> const& msk,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const int ii = rr*i;
    const int jj = rr*j;
    const int kk = rr*k;
    if (msk(ii,jj,kk)) {
        const T facx = T(0.25)*dxinv[0];
        const T facy = T(0.25)*dxinv[1];
        const T facz = T(0.25)*dxinv[2];

        T Df = T(0.0);

        const int ilo = amrex::max(ii-rr+1, fvbx.smallEnd(0));
        const int ihi = amrex::min(ii+rr-1, fvbx.bigEnd  (0));
        const int jlo = amrex::max(jj-rr+1, fvbx.smallEnd(1));
        const int jhi = amrex::min(jj+rr-1, fvbx.bigEnd  (1));
        const int klo = amrex::max(kk-rr+1, fvbx.smallEnd(2));
        const int khi = amrex::min(kk+rr-1, fvbx.bigEnd  (2));

        for (int koff = klo; koff <= khi; ++koff) {
        for (int joff = jlo; joff <= jhi; ++joff) {
        for (int ioff = ilo; ioff <= ihi; ++ioff) {
            T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                           (rr-std::abs(jj-joff)) *
                                           (rr-std::abs(kk-koff)));
            if (fvbx.strictly_contains(ioff,joff,koff)) {
                Df += scale * frhs(ioff,joff,koff);
            } else {
                Df += scale * mlndlap_sum_Df(ioff, joff, koff, facx, facy, facz, vel, velbx);
            }
        }}}

        rhs(i,j,k) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr*rr*rr));
    } else {
        rhs(i,j,k) = T(0.0);
    }
}

template<int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int i, int j, int k, Box const& ccbx,
                                Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    const int ii = rr*i;
    const int jj = rr*j;
    const int kk = rr*k;
    if (msk(ii,jj,kk)) {
        T tmp = T(0.0);

        const int ilo = amrex::max(ii-rr  , ccbx.smallEnd(0));
        const int ihi = amrex::min(ii+rr-1, ccbx.bigEnd  (0));
        const int jlo = amrex::max(jj-rr  , ccbx.smallEnd(1));
        const int jhi = amrex::min(jj+rr-1, ccbx.bigEnd  (1));
        const int klo = amrex::max(kk-rr  , ccbx.smallEnd(2));
        const int khi = amrex::min(kk+rr-1, ccbx.bigEnd (2));

        for (int koff = klo; koff <= khi; ++koff) {
        for (int joff = jlo; joff <= jhi; ++joff) {
        for (int ioff = ilo; ioff <= ihi; ++ioff) {
            T scale = (static_cast<T>(rr)-std::abs(static_cast<T>(ioff-ii)+T(0.5)))
                *        (static_cast<T>(rr)-std::abs(static_cast<T>(joff-jj)+T(0.5)))
                *        (static_cast<T>(rr)-std::abs(static_cast<T>(koff-kk)+T(0.5)));
            tmp += cc(ioff,joff,koff) * scale;
        }}}

        rhs(i,j,k) += tmp * (T(1.0)/T(rr*rr*rr*rr*rr*rr));
    }
}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T neumann_scale (int i, int j, int k, Box const& nddom,
                        GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                        GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
    {
        T val = T(1.0);

        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);

        if ((i == ndlo.x && ( bclo[0] == LinOpBCType::Neumann ||
                              bclo[0] == LinOpBCType::inflow)) ||
            (i == ndhi.x && ( bchi[0] == LinOpBCType::Neumann ||
                              bchi[0] == LinOpBCType::inflow))) {
            val *= T(2.);
        }

        if ((j == ndlo.y && ( bclo[1] == LinOpBCType::Neumann ||
                              bclo[1] == LinOpBCType::inflow)) ||
            (j == ndhi.y && ( bchi[1] == LinOpBCType::Neumann ||
                              bchi[1] == LinOpBCType::inflow))) {
            val *= T(2.);
        }

        if ((k == ndlo.z && ( bclo[2] == LinOpBCType::Neumann ||
                              bclo[2] == LinOpBCType::inflow)) ||
            (k == ndhi.z && ( bchi[2] == LinOpBCType::Neumann ||
                              bchi[2] == LinOpBCType::inflow))) {
            val *= T(2.);
        }

        return val;
    }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int i, int j, int k, Array4<T> const& rhs,
                              Array4<std::type_identity_t<T> const> const& vel, Array4<std::type_identity_t<T> const> const& fc,
                              Array4<std::type_identity_t<T> const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                              Box const& ccdom_p, Box const& veldom, Box const& nddom,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,k) && ndmsk(i,j,k) == crse_fine_node) {
        T facx = T(0.25) * dxinv[0];
        T facy = T(0.25) * dxinv[1];
        T facz = T(0.25) * dxinv[2];
        T tmp = fc(i,j,k);

        // Where there is inflow, veldom there is bigger than ccdom_p by one cell.
        // ccdom_p is cc domain grown at periodic boundaries.

        if (ccmsk(i-1,j-1,k-1) == crse_cell && veldom.contains(i-1,j-1,k-1)) {
            tmp += -facx*vel(i-1,j-1,k-1,0) - facy*vel(i-1,j-1,k-1,1) - facz*vel(i-1,j-1,k-1,2);
            if (rhcc && ccdom_p.contains(i-1,j-1,k-1)) {
                tmp += T(0.125) * rhcc(i-1,j-1,k-1);
            }
        }

        if (ccmsk(i,j-1,k-1) == crse_cell && veldom.contains(i,j-1,k-1)) {
            tmp += facx*vel(i,j-1,k-1,0) - facy*vel(i,j-1,k-1,1) - facz*vel(i,j-1,k-1,2);
            if (rhcc && ccdom_p.contains(i,j-1,k-1)) {
                tmp += T(0.125) * rhcc(i,j-1,k-1);
            }
        }

        if (ccmsk(i-1,j,k-1) == crse_cell && veldom.contains(i-1,j,k-1)) {
            tmp += -facx*vel(i-1,j,k-1,0) + facy*vel(i-1,j,k-1,1) - facz*vel(i-1,j,k-1,2);
            if (rhcc && ccdom_p.contains(i-1,j,k-1)) {
                tmp += T(0.125) * rhcc(i-1,j,k-1);
            }
        }

        if (ccmsk(i,j,k-1) == crse_cell && veldom.contains(i,j,k-1)) {
            tmp += facx*vel(i,j,k-1,0) + facy*vel(i,j,k-1,1) - facz*vel(i,j,k-1,2);
            if (rhcc && ccdom_p.contains(i,j,k-1)) {
                tmp += T(0.125) * rhcc(i,j,k-1);
            }
        }

        if (ccmsk(i-1,j-1,k) == crse_cell && veldom.contains(i-1,j-1,k)) {
            tmp += -facx*vel(i-1,j-1,k,0) - facy*vel(i-1,j-1,k,1) + facz*vel(i-1,j-1,k,2);
            if (rhcc && ccdom_p.contains(i-1,j-1,k)) {
                tmp += T(0.125) * rhcc(i-1,j-1,k);
            }
        }

        if (ccmsk(i,j-1,k) == crse_cell && veldom.contains(i,j-1,k)) {
            tmp += facx*vel(i,j-1,k,0) - facy*vel(i,j-1,k,1) + facz*vel(i,j-1,k,2);
            if (rhcc && ccdom_p.contains(i,j-1,k)) {
                tmp += T(0.125) * rhcc(i,j-1,k);
            }
        }

        if (ccmsk(i-1,j,k) == crse_cell && veldom.contains(i-1,j,k)) {
            tmp += -facx*vel(i-1,j,k,0) + facy*vel(i-1,j,k,1) + facz*vel(i-1,j,k,2);
            if (rhcc && ccdom_p.contains(i-1,j,k)) {
                tmp += T(0.125) * rhcc(i-1,j,k);
            }
        }

        if (ccmsk(i,j,k) == crse_cell && veldom.contains(i,j,k)) {
            tmp += facx*vel(i,j,k,0) + facy*vel(i,j,k,1) + facz*vel(i,j,k,2);
            if (rhcc && ccdom_p.contains(i,j,k)) {
                tmp += T(0.125) * rhcc(i,j,k);
            }
        }

        rhs(i,j,k) = tmp * neumann_scale<T>(i, j, k, nddom, bclo, bchi);
    }
}

//
// residual
//

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int i, int j, int k, Array4<T> const& resid,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                         bool neumann_doubling) noexcept
{
    if ( msk(i-1,j-1,k-1) == 0 ||
         msk(i  ,j-1,k-1) == 0 ||
         msk(i-1,j  ,k-1) == 0 ||
         msk(i  ,j  ,k-1) == 0 ||
         msk(i-1,j-1,k  ) == 0 ||
         msk(i  ,j-1,k  ) == 0 ||
         msk(i-1,j  ,k  ) == 0 ||
         msk(i  ,j  ,k  ) == 0 )
    {
        T fac = T(1.0);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if ((i == ndlo.x && ( bclo[0] == LinOpBCType::Neumann ||
                                  bclo[0] == LinOpBCType::inflow)) ||
                (i == ndhi.x && ( bchi[0] == LinOpBCType::Neumann ||
                                  bchi[0] == LinOpBCType::inflow))) {
                fac *= T(2.);
            }
            if ((j == ndlo.y && ( bclo[1] == LinOpBCType::Neumann ||
                                  bclo[1] == LinOpBCType::inflow)) ||
                (j == ndhi.y && ( bchi[1] == LinOpBCType::Neumann ||
                                  bchi[1] == LinOpBCType::inflow))) {
                fac *= T(2.);
            }
            if ((k == ndlo.z && ( bclo[2] == LinOpBCType::Neumann ||
                                  bclo[2] == LinOpBCType::inflow)) ||
                (k == ndhi.z && ( bchi[2] == LinOpBCType::Neumann ||
                                  bchi[2] == LinOpBCType::inflow))) {
                fac *= T(2.);
            }
        }
        resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac;
    } else {
        resid(i,j,k) = T(0.);
    }
}

//
// sync residual
//

    template <typename P, typename S, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Ax (P const& pred, S const& sig,
                         int i, int j, int k, std::type_identity_t<T> facx, std::type_identity_t<T> facy, std::type_identity_t<T> facz,
                         Array4<T const> const& phi) noexcept
    {
        T Ax = T(0.0);
        if (pred(i-1,j-1,k-1)) {
            Ax += sig(i-1,j-1,k-1)*(facx*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                         +T(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                         +T(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                         +         (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                  + facy*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                         +T(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                         +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                         +         (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1)))
                                  + facz*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                         +T(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                         +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                         +         (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  ))));
        }
        if (pred(i,j-1,k-1)) {
            Ax += sig(i,j-1,k-1)*(facx*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +T(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                + facy*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1)))
                                + facz*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                       +         (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  ))));
        }
        if (pred(i-1,j,k-1)) {
            Ax += sig(i-1,j,k-1)*(facx*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                       +T(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1)))
                                 + facy*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                        +T(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                        +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                        +         (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1)))
                                 + facz*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                        +T(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                        +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                        +         (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  ))));
        }
        if (pred(i,j,k-1)) {
            Ax += sig(i,j,k-1)*(facx*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +T(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1)))
                              + facy*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1)))
                              + facz*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                     +         (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  ))));
        }
        if (pred(i-1,j-1,k)) {
            Ax += sig(i-1,j-1,k)*(facx*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +T(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                       +         (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1)))
                                + facy*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                       +         (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1)))
                                + facz*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                       +         (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  ))));
        }
        if (pred(i,j-1,k)) {
            Ax += sig(i,j-1,k)*(facx*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                     +T(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1)))
                              + facy*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1)))
                              + facz*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                     +         (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  ))));
        }
        if (pred(i-1,j,k)) {
            Ax += sig(i-1,j,k)*(facx*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +T(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1)))
                              + facy*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1)))
                              + facz*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                     +         (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  ))));
        }
        if (pred(i,j,k)) {
            Ax += sig(i,j,k)*(facx*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                   +T(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1)))
                            + facy*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1)))
                            + facz*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                   +         (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  ))));
        }
        return Ax;
    }

    template <int rr, typename S, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mlndlap_Ax_fine_contrib_doit (S const& sig,
                                       int i, int j, int k, Box const& ndbx, Box const& ccbx,
                                       Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                       Array4<std::type_identity_t<T> const> const& rhs,
                                       Array4<std::type_identity_t<T> const> const& phi,
                                       Array4<int const> const& msk,
                                       GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
    {
        const int ii = rr*i;
        const int jj = rr*j;
        const int kk = rr*k;
        if (msk(ii,jj,kk)) {
            T facx = T(1./36.)*dxinv[0]*dxinv[0];
            T facy = T(1./36.)*dxinv[1]*dxinv[1];
            T facz = T(1./36.)*dxinv[2]*dxinv[2];

            auto is_fine = [&ccbx] (int ix, int iy, int iz) -> bool {
                return ccbx.contains(ix,iy,iz);
            };

            T Df = T(0.0);

            const int ilo = amrex::max(ii-rr+1, ndbx.smallEnd(0));
            const int ihi = amrex::min(ii+rr-1, ndbx.bigEnd  (0));
            const int jlo = amrex::max(jj-rr+1, ndbx.smallEnd(1));
            const int jhi = amrex::min(jj+rr-1, ndbx.bigEnd  (1));
            const int klo = amrex::max(kk-rr+1, ndbx.smallEnd(2));
            const int khi = amrex::min(kk+rr-1, ndbx.bigEnd  (2));

            for (int koff = klo; koff <= khi; ++koff) {
            for (int joff = jlo; joff <= jhi; ++joff) {
            for (int ioff = ilo; ioff <= ihi; ++ioff) {
                T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                               (rr-std::abs(jj-joff)) *
                                               (rr-std::abs(kk-koff)));
                if (ndbx.strictly_contains(ioff,joff,koff)) {
                    Df += scale * (rhs(ioff,joff,koff)-res(ioff,joff,koff));
                } else {
                    Df += scale * mlndlap_sum_Ax
                        (is_fine, sig, ioff, joff, koff, facx, facy, facz, phi);
                }
            }}}

            f(i,j,k) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr*rr*rr));
        } else {
            f(i,j,k) = T(0.0);
        }
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib (int i, int j, int k, Box const& ndbx, Box const& ccbx,
                              Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                              Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                              Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& msk,
                              GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit<rr>
        ([&sig] (int ix, int iy, int iz) -> T const& { return sig(ix,iy,iz); },
         i,j,k,ndbx,ccbx,f,res,rhs,phi,msk,dxinv);
}

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib_cs (int i, int j, int k, Box const& ndbx, Box const& ccbx,
                                 Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                 Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                                 std::type_identity_t<T> sig, Array4<int const> const& msk,
                                 GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit<rr>
        ([=] (int, int, int) -> T { return sig; },
         i,j,k,ndbx,ccbx,f,res,rhs,phi,msk,dxinv);
}

// Anisotropic-sigma 3D variants of mlndlap_sum_Ax and
// mlndlap_Ax_fine_contrib used by the AMR-sync paths when m_use_mapped.
// Per-axis sigma: facx-terms use sx, facy-terms use sy, facz-terms use sz.
// Recover the scalar kernels when sx == sy == sz.
    template <typename P, typename SX, typename SY, typename SZ, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T mlndlap_sum_Ax_ha (P const& pred, SX const& sx, SY const& sy, SZ const& sz,
                            int i, int j, int k, std::type_identity_t<T> facx, std::type_identity_t<T> facy, std::type_identity_t<T> facz,
                            Array4<T const> const& phi) noexcept
    {
        T Ax = T(0.0);
        if (pred(i-1,j-1,k-1)) {
            Ax += facx*sx(i-1,j-1,k-1)*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +T(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1)))
                + facy*sy(i-1,j-1,k-1)*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                       +         (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1)))
                + facz*sz(i-1,j-1,k-1)*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                       +T(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                       +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                       +         (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  )));
        }
        if (pred(i,j-1,k-1)) {
            Ax += facx*sx(i,j-1,k-1)*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                     +T(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1)))
                + facy*sy(i,j-1,k-1)*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1)))
                + facz*sz(i,j-1,k-1)*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                     +         (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  )));
        }
        if (pred(i-1,j,k-1)) {
            Ax += facx*sx(i-1,j,k-1)*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +T(2.)*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1)))
                + facy*sy(i-1,j,k-1)*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                     +         (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1)))
                + facz*sz(i-1,j,k-1)*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                     +         (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  )));
        }
        if (pred(i,j,k-1)) {
            Ax += facx*sx(i,j,k-1)*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                   +T(2.)*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                   +         (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1)))
                + facy*sy(i,j,k-1)*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                   +         (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1)))
                + facz*sz(i,j,k-1)*(T(4.)*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                   +         (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  )));
        }
        if (pred(i-1,j-1,k)) {
            Ax += facx*sx(i-1,j-1,k)*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                     +T(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1)))
                + facy*sy(i-1,j-1,k)*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                     +         (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1)))
                + facz*sz(i-1,j-1,k)*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +T(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                     +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                     +         (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  )));
        }
        if (pred(i,j-1,k)) {
            Ax += facx*sx(i,j-1,k)*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                   +T(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1)))
                + facy*sy(i,j-1,k)*(T(4.)*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1)))
                + facz*sz(i,j-1,k)*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                   +         (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  )));
        }
        if (pred(i-1,j,k)) {
            Ax += facx*sx(i-1,j,k)*(T(4.)*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                   +T(2.)*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1)))
                + facy*sy(i-1,j,k)*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                   +         (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1)))
                + facz*sz(i-1,j,k)*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                   +T(2.)*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                   +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                   +         (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  )));
        }
        if (pred(i,j,k)) {
            Ax += facx*sx(i,j,k)*(T(4.)*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                 +T(2.)*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                 +T(2.)*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                 +         (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1)))
                + facy*sy(i,j,k)*(T(4.)*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                 +T(2.)*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                 +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                 +         (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1)))
                + facz*sz(i,j,k)*(T(4.)*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                 +T(2.)*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                 +T(2.)*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                 +         (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  )));
        }
        return Ax;
    }

    template <int rr, typename SX, typename SY, typename SZ, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mlndlap_Ax_fine_contrib_doit_ha (SX const& sx, SY const& sy, SZ const& sz,
                                          int i, int j, int k, Box const& ndbx, Box const& ccbx,
                                          Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                          Array4<std::type_identity_t<T> const> const& rhs,
                                          Array4<std::type_identity_t<T> const> const& phi,
                                          Array4<int const> const& msk,
                                          GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
    {
        const int ii = rr*i;
        const int jj = rr*j;
        const int kk = rr*k;
        if (msk(ii,jj,kk)) {
            T facx = T(1./36.)*dxinv[0]*dxinv[0];
            T facy = T(1./36.)*dxinv[1]*dxinv[1];
            T facz = T(1./36.)*dxinv[2]*dxinv[2];

            auto is_fine = [&ccbx] (int ix, int iy, int iz) -> bool {
                return ccbx.contains(ix,iy,iz);
            };

            T Df = T(0.0);

            const int ilo = amrex::max(ii-rr+1, ndbx.smallEnd(0));
            const int ihi = amrex::min(ii+rr-1, ndbx.bigEnd  (0));
            const int jlo = amrex::max(jj-rr+1, ndbx.smallEnd(1));
            const int jhi = amrex::min(jj+rr-1, ndbx.bigEnd  (1));
            const int klo = amrex::max(kk-rr+1, ndbx.smallEnd(2));
            const int khi = amrex::min(kk+rr-1, ndbx.bigEnd  (2));

            for (int koff = klo; koff <= khi; ++koff) {
            for (int joff = jlo; joff <= jhi; ++joff) {
            for (int ioff = ilo; ioff <= ihi; ++ioff) {
                T scale = static_cast<T>((rr-std::abs(ii-ioff)) *
                                               (rr-std::abs(jj-joff)) *
                                               (rr-std::abs(kk-koff)));
                if (ndbx.strictly_contains(ioff,joff,koff)) {
                    Df += scale * (rhs(ioff,joff,koff)-res(ioff,joff,koff));
                } else {
                    Df += scale * mlndlap_sum_Ax_ha
                        (is_fine, sx, sy, sz, ioff, joff, koff, facx, facy, facz, phi);
                }
            }}}

            f(i,j,k) = Df * (T(1.0)/static_cast<T>(rr*rr*rr*rr*rr*rr));
        } else {
            f(i,j,k) = T(0.0);
        }
    }

template <int rr, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_Ax_fine_contrib_ha (int i, int j, int k, Box const& ndbx, Box const& ccbx,
                                 Array4<T> const& f, Array4<std::type_identity_t<T> const> const& res,
                                 Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& phi,
                                 Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                                 Array4<std::type_identity_t<T> const> const& sz, Array4<int const> const& msk,
                                 GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    mlndlap_Ax_fine_contrib_doit_ha<rr>
        ([&sx] (int ix, int iy, int iz) -> T const& { return sx(ix,iy,iz); },
         [&sy] (int ix, int iy, int iz) -> T const& { return sy(ix,iy,iz); },
         [&sz] (int ix, int iy, int iz) -> T const& { return sz(ix,iy,iz); },
         i,j,k,ndbx,ccbx,f,res,rhs,phi,msk,dxinv);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int i, int j, int k, Array4<T> const& res,
                             Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                             Array4<std::type_identity_t<T> const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<std::type_identity_t<T> const> const& fc,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                             Box const& ccdom_p, Box const& nddom,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                             bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,k) && ndmsk(i,j,k) == crse_fine_node) {
        T facx = T(1./36.)*dxinv[0]*dxinv[0];
        T facy = T(1./36.)*dxinv[1]*dxinv[1];
        T facz = T(1./36.)*dxinv[2]*dxinv[2];

        T Ax = mlndlap_sum_Ax([&ccmsk, &ccdom_p] (int ix, int iy, int iz) -> bool
                                 {
                                     return ccdom_p.contains(ix,iy,iz)
                                         && (ccmsk(ix,iy,iz) == crse_cell);
                                 },
                                 [&sig] (int ix, int iy, int iz) -> T const&
                                 {
                                     return sig(ix,iy,iz);
                                 },
                                 i, j, k, facx, facy, facz, phi);
        Ax += fc(i,j,k);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,k,nddom,bclo,bchi) : T(1.0);
        res(i,j,k) = rhs(i,j,k) - Ax*ns;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_ha (int i, int j, int k, Array4<T> const& res,
                                Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sx, Array4<std::type_identity_t<T> const> const& sy,
                                Array4<std::type_identity_t<T> const> const& sz,
                                Array4<int const> const& dmsk,
                                Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                                Array4<std::type_identity_t<T> const> const& fc,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                                Box const& ccdom_p, Box const& nddom,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                                bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,k) && ndmsk(i,j,k) == crse_fine_node) {
        T facx = T(1./36.)*dxinv[0]*dxinv[0];
        T facy = T(1./36.)*dxinv[1]*dxinv[1];
        T facz = T(1./36.)*dxinv[2]*dxinv[2];

        T Ax = mlndlap_sum_Ax_ha([&ccmsk, &ccdom_p] (int ix, int iy, int iz) -> bool
                                    {
                                        return ccdom_p.contains(ix,iy,iz)
                                            && (ccmsk(ix,iy,iz) == crse_cell);
                                    },
                                    [&sx] (int ix, int iy, int iz) -> T const&
                                    {
                                        return sx(ix,iy,iz);
                                    },
                                    [&sy] (int ix, int iy, int iz) -> T const&
                                    {
                                        return sy(ix,iy,iz);
                                    },
                                    [&sz] (int ix, int iy, int iz) -> T const&
                                    {
                                        return sz(ix,iy,iz);
                                    },
                                    i, j, k, facx, facy, facz, phi);
        Ax += fc(i,j,k);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,k,nddom,bclo,bchi) : T(1.0);
        res(i,j,k) = rhs(i,j,k) - Ax*ns;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_cs (int i, int j, int k, Array4<T> const& res,
                                Array4<std::type_identity_t<T> const> const& phi, Array4<std::type_identity_t<T> const> const& rhs,
                                std::type_identity_t<T> sig, Array4<int const> const& dmsk,
                                Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                                Array4<std::type_identity_t<T> const> const& fc,
                                GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                                Box const& ccdom_p, Box const& nddom,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                                bool neumann_doubling) noexcept
{
    using namespace nodelap_detail;

    if (!dmsk(i,j,k) && ndmsk(i,j,k) == crse_fine_node) {
        T facx = T(1./36.)*dxinv[0]*dxinv[0];
        T facy = T(1./36.)*dxinv[1]*dxinv[1];
        T facz = T(1./36.)*dxinv[2]*dxinv[2];

        T Ax = mlndlap_sum_Ax([&ccmsk, &ccdom_p] (int ix, int iy, int iz) -> bool
                                 {
                                     return ccdom_p.contains(ix,iy,iz)
                                         && (ccmsk(ix,iy,iz) == crse_cell);
                                 },
                                 [=] (int, int, int) -> T
                                 {
                                     return sig;
                                 },
                                 i, j, k, facx, facy, facz, phi);
        Ax += fc(i,j,k);
        T const ns = (neumann_doubling) ? neumann_scale<T>(i,j,k,nddom,bclo,bchi) : T(1.0);
        res(i,j,k) = rhs(i,j,k) - Ax*ns;
    }
}

//
// RAP
//

/// \cond DOXYGEN_IGNORE
namespace nodelap_detail {

    constexpr int ist_000 = 0;
    constexpr int ist_p00 = 1;
    constexpr int ist_0p0 = 2;
    constexpr int ist_00p = 3;
    constexpr int ist_pp0 = 4;
    constexpr int ist_p0p = 5;
    constexpr int ist_0pp = 6;
    constexpr int ist_ppp = 7;
    constexpr int ist_inv = 8;
    constexpr int n_sten  = 9;
}
/// \endcond

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const& bx, Array4<T> const& sten,
                          Array4<std::type_identity_t<T> const> const& sig,
                          GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    using namespace nodelap_detail;

    T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
    T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
    T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
    T fxyz = facx + facy + facz;
    T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
    T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
    T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
    T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
    T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
    T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        // i+1,j,k
        sten(i,j,k,ist_p00) = f4xm2ym2z * (sig(i,j-1,k-1)+sig(i,j,k-1)+sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k: sten(i-1,j,k,ist_p00)

        // i,j+1,k
        sten(i,j,k,ist_0p0) = fm2x4ym2z * (sig(i-1,j,k-1)+sig(i,j,k-1)+sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k: sten(i,j-1,k,ist_0p0)

        // i,j,k+1
        sten(i,j,k,ist_00p) = fm2xm2y4z * (sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        // i,j,k-1: sten(i,j,k-1,ist_00p)

        // i+1,j+1,k
        sten(i,j,k,ist_pp0) = f2x2ymz * (sig(i,j,k-1)+sig(i,j,k));
        // i-1,j-1,k: sten(i-1,j-1,k,ist_pp0)
        // i+1,j-1,k: sten(i  ,j-1,k,ist_pp0)
        // i-1,j+1,k: sten(i-1,j  ,k,ist_pp0)

        // i+1,j,k+1
        sten(i,j,k,ist_p0p) = f2xmy2z * (sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k-1: sten(i-1,j,k-1,ist_p0p)
        // i+1,j,k-1: sten(i  ,j,k-1,ist_p0p)
        // i-1,j,k+1: sten(i-1,j,k  ,ist_p0p)

        // i,j+1,k+1
        sten(i,j,k,ist_0pp) = fmx2y2z * (sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k-1: sten(i,j-1,k-1,ist_0pp)
        // i,j+1,k-1: sten(i,j  ,k-1,ist_0pp)
        // i,j-1,k+1: sten(i,j-1,k  ,ist_0pp)

        // i+1,j+1,k+1
        sten(i,j,k,ist_ppp) = fxyz * sig(i,j,k);
        // i-1,j-1,k-1: sten(i-1,j-1,k-1,ist_ppp)
        // i+1,j-1,k-1: sten(i  ,j-1,k-1,ist_ppp)
        // i-1,j+1,k-1: sten(i-1,j  ,k-1,ist_ppp)
        // i+1,j+1,k-1: sten(i  ,j  ,k-1,ist_ppp)
        // i-1,j-1,k+1: sten(i-1,j-1,k  ,ist_ppp)
        // i+1,j-1,k+1: sten(i  ,j-1,k  ,ist_ppp)
        // i-1,j+1,k+1: sten(i-1,j  ,k  ,ist_ppp)
    });
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int i, int j, int k, Array4<T> const& sten) noexcept
{
    using namespace nodelap_detail;

    sten(i,j,k,ist_000) = -(sten(i-1,j,k,ist_p00) + sten(i,j,k,ist_p00)
                          + sten(i,j-1,k,ist_0p0) + sten(i,j,k,ist_0p0)
                          + sten(i,j,k-1,ist_00p) + sten(i,j,k,ist_00p)
                          + sten(i-1,j-1,k,ist_pp0) + sten(i,j-1,k,ist_pp0)
                          + sten(i-1,j,k,ist_pp0) + sten(i,j,k,ist_pp0)
                          + sten(i-1,j,k-1,ist_p0p) + sten(i,j,k-1,ist_p0p)
                          + sten(i-1,j,k,ist_p0p) + sten(i,j,k,ist_p0p)
                          + sten(i,j-1,k-1,ist_0pp) + sten(i,j,k-1,ist_0pp)
                          + sten(i,j-1,k,ist_0pp) + sten(i,j,k,ist_0pp)
                          + sten(i-1,j-1,k-1,ist_ppp) + sten(i,j-1,k-1,ist_ppp)
                          + sten(i-1,j,k-1,ist_ppp) + sten(i,j,k-1,ist_ppp)
                          + sten(i-1,j-1,k,ist_ppp) + sten(i,j-1,k,ist_ppp)
                          + sten(i-1,j,k,ist_ppp) + sten(i,j,k,ist_ppp));
    sten(i,j,k,ist_inv) = T(1.0) /
        (  std::abs(sten(i-1,j,k,ist_p00)) + std::abs(sten(i,j,k,ist_p00))
         + std::abs(sten(i,j-1,k,ist_0p0)) + std::abs(sten(i,j,k,ist_0p0))
         + std::abs(sten(i,j,k-1,ist_00p)) + std::abs(sten(i,j,k,ist_00p))
         + std::abs(sten(i-1,j-1,k,ist_pp0)) + std::abs(sten(i,j-1,k,ist_pp0))
         + std::abs(sten(i-1,j,k,ist_pp0)) + std::abs(sten(i,j,k,ist_pp0))
         + std::abs(sten(i-1,j,k-1,ist_p0p)) + std::abs(sten(i,j,k-1,ist_p0p))
         + std::abs(sten(i-1,j,k,ist_p0p)) + std::abs(sten(i,j,k,ist_p0p))
         + std::abs(sten(i,j-1,k-1,ist_0pp)) + std::abs(sten(i,j,k-1,ist_0pp))
         + std::abs(sten(i,j-1,k,ist_0pp)) + std::abs(sten(i,j,k,ist_0pp))
         + std::abs(sten(i-1,j-1,k-1,ist_ppp)) + std::abs(sten(i,j-1,k-1,ist_ppp))
         + std::abs(sten(i-1,j,k-1,ist_ppp)) + std::abs(sten(i,j,k-1,ist_ppp))
         + std::abs(sten(i-1,j-1,k,ist_ppp)) + std::abs(sten(i,j-1,k,ist_ppp))
         + std::abs(sten(i-1,j,k,ist_ppp)) + std::abs(sten(i,j,k,ist_ppp)) + eps_v<T>);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int i, int j, int k, Array4<T> const& csten,
                          Array4<std::type_identity_t<T> const> const& fsten) noexcept
{
    using namespace nodelap_detail;

    auto interp_from_mmm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };
    amrex::ignore_unused(interp_from_mmm_to);

    auto interp_from_pmm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mpm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mmp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_-1,k_,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_pmp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_  ,j_  ,k_,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) * fsten(i_,j_  ,k_,ist_inv);
        return p;
    };

    auto interp_from_mpp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T p = T(1.);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps_v<T>);
        p *= std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_0mm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_0mp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_0pm_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_0pp_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_m0m_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_p0m_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_m0p_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_p0p_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_mm0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps_v<T>);
      T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                            +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_mp0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_pm0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_pp0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto interp_from_00m_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_00p_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_0m0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_0p0_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_m00_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_p00_to = [&fsten] (int i_, int j_, int k_) -> T {
        T w1 = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        T w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto Ammm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_-1,k_-1,ist_ppp);
    };
    amrex::ignore_unused(Ammm);

    auto A0mm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_-1,ist_0pp);
    };

    auto Apmm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_-1,ist_ppp);
    };

    auto Am0m = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_-1,ist_p0p);
    };

    auto A00m = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_-1,ist_00p);
    };

    auto Ap0m = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_-1,ist_p0p);
    };

    auto Ampm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_-1,ist_ppp);
    };

    auto A0pm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_-1,ist_0pp);
    };

    auto Appm = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_-1,ist_ppp);
    };

    auto Amm0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_-1,k_  ,ist_pp0);
    };

    auto A0m0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_  ,ist_0p0);
    };

    auto Apm0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_  ,ist_pp0);
    };

    auto Am00 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_  ,ist_p00);
    };

    auto A000 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_000);
    };

    auto Ap00 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_p00);
    };

    auto Amp0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_  ,ist_pp0);
    };

    auto A0p0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_0p0);
    };

    auto App0 = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_pp0);
    };

    auto Ammp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_-1,k_  ,ist_ppp);
    };

    auto A0mp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_  ,ist_0pp);
    };

    auto Apmp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_-1,k_  ,ist_ppp);
    };

    auto Am0p = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_  ,ist_p0p);
    };

    auto A00p = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_00p);
    };

    auto Ap0p = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_p0p);
    };

    auto Ampp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_-1,j_  ,k_  ,ist_ppp);
    };

    auto A0pp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_0pp);
    };

    auto Appp = [&fsten] (int i_, int j_, int k_) -> T {
        return fsten(i_  ,j_  ,k_  ,ist_ppp);
    };

    auto restrict_from_mmm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) * fsten(ii_-1,jj_-1,kk_-1,ist_inv);
        return r;
    };
    amrex::ignore_unused(restrict_from_mmm_to);

    auto restrict_from_0mm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_,jj_-2,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_,jj_-1,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_,jj_-1,kk_-2,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_,jj_-1,kk_-1,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_pmm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) * fsten(ii_+1,jj_-1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_m0m_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_-2,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_-1,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_-1,jj_,kk_-2,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_-1,jj_,kk_-1,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_00m_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_,jj_,kk_-2,ist_00p));
        T w2 = std::abs(fsten(ii_,jj_,kk_-1,ist_00p));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_p0m_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_  ,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_+1,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_+1,jj_,kk_-2,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_+1,jj_,kk_-1,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_mpm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-2,jj_+1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) * fsten(ii_-1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_0pm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_,jj_  ,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_,jj_+1,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_,jj_+1,kk_-2,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_,jj_+1,kk_-1,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_ppm_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_+1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) * fsten(ii_+1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_mm0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_-2,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_-1,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_-1,jj_-2,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_-1,jj_-1,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) * (T(1.) + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_0m0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_,jj_-2,kk_,ist_0p0));
        T w2 = std::abs(fsten(ii_,jj_-1,kk_,ist_0p0));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_pm0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_  ,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_+1,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_+1,jj_-2,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_+1,jj_-1,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) * (T(1.) + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_m00_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_-2,jj_,kk_,ist_p00));
        T w2 = std::abs(fsten(ii_-1,jj_,kk_,ist_p00));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_000_to = [] (int /*ii_*/, int /*jj_*/, int /*kk_*/) -> T {
        return T(1.);
    };

    auto restrict_from_p00_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_  ,jj_,kk_,ist_p00));
        T w2 = std::abs(fsten(ii_+1,jj_,kk_,ist_p00));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_mp0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_-2,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_-1,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_-1,jj_  ,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_-1,jj_+1,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_0p0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_,jj_  ,kk_,ist_0p0));
        T w2 = std::abs(fsten(ii_,jj_+1,kk_,ist_0p0));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_pp0_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_  ,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_+1,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_+1,jj_  ,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_+1,jj_+1,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_mmp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) * fsten(ii_-1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0mp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_,jj_-2,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_,jj_-1,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_,jj_-1,kk_  ,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_,jj_-1,kk_+1,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_pmp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-2,kk_+1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) * fsten(ii_+1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_m0p_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_-2,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_-1,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_-1,jj_,kk_  ,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_-1,jj_,kk_+1,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) * (T(1.) + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_00p_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1 = std::abs(fsten(ii_,jj_,kk_  ,ist_00p));
        T w2 = std::abs(fsten(ii_,jj_,kk_+1,ist_00p));
        if (w1 == T(0.) && w2 == T(0.)) {
            return T(0.5);
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_p0p_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_  ,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_+1,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_+1,jj_,kk_  ,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_+1,jj_,kk_+1,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_mpp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_-1,jj_+1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) * fsten(ii_-1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0pp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T w1m = std::abs(fsten(ii_,jj_  ,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) + eps_v<T>);
        T w1p = std::abs(fsten(ii_,jj_+1,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps_v<T>);
        T w2m = std::abs(fsten(ii_,jj_+1,kk_  ,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) + eps_v<T>);
        T w2p = std::abs(fsten(ii_,jj_+1,kk_+1,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps_v<T>);
        T wmm = std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp)) * (T(1.) + w1m + w2m);
        T wpm = std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) * (T(1.) + w1p + w2m);
        T wmp = std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) * (T(1.) + w1m + w2p);
        T wpp = std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) * (T(1.) + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps_v<T>);
    };

    auto restrict_from_ppp_to = [&fsten] (int ii_, int jj_, int kk_) -> T {
        T r = T(1.);
        r += std::abs(fsten(ii_  ,jj_+1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_  ,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_+1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_+1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp)) + eps_v<T>);
        r += std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp)) + eps_v<T>);
        r *= std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp)) * fsten(ii_+1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    Array3D<T,-1,1,-1,1,-1,1> p;
    Array3D<T,-1,1,-1,1,-1,1> ap;
    T cs1, cs2, cs3, cs4;

    // csten(i,j,k,ist_p00)
    int iii = ii;
    int jjj = jj;
    int kkk = kk;;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk-1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk-1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk  );
    p( 0, 0, 0) = T(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk  );
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,+1) = interp_from_pmm_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,+1) = interp_from_0mm_to(iii+2,jjj+1,kkk+1);
    ap(0,-1,-1) =
                     Ap00(iii,jjj-1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk-1) * p(-1, 0, 0);
    ap(1,-1,-1) =
                     A000(iii+1,jjj-1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk-1) * p( 0, 0, 0);
    ap(0,0,-1) =
                     Apm0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk-1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk-1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk-1) * p(-1,+1, 0);
    ap(1,0,-1) =
                     A0m0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,+1, 0);
    ap(0,1,-1) =
                     Apm0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,+1, 0);
    ap(1,1,-1) =
                     A0m0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,+1, 0);
    ap(0,-1,0) =
                     Ap0m(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj-1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,+1);
    ap(1,-1,0) =
                     A00m(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj-1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj-1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj-1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,+1);
    ap(0,0,0) =
                     Apmm(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0m(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appm(iii,jjj,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Ap00(iii,jjj,kkk) * p(-1, 0, 0)
      +              App0(iii,jjj,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj,kkk) * p(-1,-1,+1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,+1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,+1);
    ap(1,0,0) =
                     A0mm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0m(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              Apm0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj,kkk) * p(-1, 0, 0)
      +              Ap00(iii+1,jjj,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj,kkk) * p(-1,+1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,+1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,+1);
    ap(0,1,0) =
                     Apmm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,+1);
    ap(1,1,0) =
                     A0mm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              Apm0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk) * p(-1,+1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,+1);
    ap(0,-1,1) =
                     Ap0m(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,+1);
    ap(1,-1,1) =
                     A00m(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,+1);
    ap(0,0,1) =
                     Apmm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appm(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,+1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,+1);
    ap(1,0,1) =
                     A0mm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,+1);
    ap(0,1,1) =
                     Apmm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,+1);
    ap(1,1,1) =
                     A0mm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,+1);
    csten(i,j,k,ist_p00) = T(0.125) *
      ( restrict_from_0mm_to(iii,jjj,kkk) * ap( 0,-1,-1)
      + restrict_from_pmm_to(iii,jjj,kkk) * ap(+1,-1,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_0p0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk  );
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii+1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     A0p0(iii-1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk-1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk-1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(1,0,-1) =
                     Amp0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     A000(iii-1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk-1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk-1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk-1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk-1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(1,1,-1) =
                     Am00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     A0pm(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pm(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(-1,-1, 0)
      +              A0p0(iii,jjj,kkk) * p( 0,-1, 0)
      +              App0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(-1,-1,+1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,+1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(1,0,0) =
                     Ampm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     A00m(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00m(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pm(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              A000(iii,jjj+1,kkk) * p( 0,-1, 0)
      +              Ap00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              A0p0(iii,jjj+1,kkk) * p( 0, 0, 0)
      +              App0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,+1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(1,1,0) =
                     Am0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     A0pm(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pm(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,+1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(1,0,1) =
                     Ampm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     A00m(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00m(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pm(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,+1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(1,1,1) =
                     Am0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,+1);
    csten(i,j,k,ist_0p0) = T(0.125) *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_00p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj  ,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj  ,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj+1,kkk+2);
    ap(-1,-1,0) =
                     A00p(iii-1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(1,-1,0) =
                     Am0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              A00p(iii,jjj,kkk) * p( 0, 0,-1)
      +              Ap0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     A0mp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(1,1,0) =
                     Ammp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A000(iii-1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(1,-1,1) =
                     Am00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              A000(iii,jjj,kkk+1) * p( 0, 0,-1)
      +              Ap00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              A00p(iii,jjj,kkk+1) * p( 0, 0, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     A0m0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(1,1,1) =
                     Amm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,+1, 0);
    csten(i,j,k,ist_00p) = T(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_pp0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk  );
    p( 0, 0, 0) = T(1.);
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj+2,kkk+1);
    ap(0,0,-1) =
                     App0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk-1) * p(-1,-1, 0);
    ap(1,0,-1) =
                     A0p0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,1,-1) =
                     Ap00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(-1, 0, 0);
    ap(1,1,-1) =
                     A000(iii+1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,0,0) =
                     Appm(iii,jjj,kkk) * p(-1,-1,-1)
      +              App0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Appp(iii,jjj,kkk) * p(-1,-1,+1);
    ap(1,0,0) =
                     A0pm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,+1);
    ap(0,1,0) =
                     Ap0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,+1);
    ap(1,1,0) =
                     A00m(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,0,1) =
                     Appm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              App0(iii,jjj,kkk+1) * p(-1,-1,+1);
    ap(1,0,1) =
                     A0pm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,1,1) =
                     Ap0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,+1);
    ap(1,1,1) =
                     A00m(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1);
    cs1 = T(0.125) *
      ( restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mp0)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk-1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk  );
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk  );
    p( 0,-1,+1) = interp_from_0pm_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii-1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     Amp0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk-1) * p(+1,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     Am00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     Ampm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     Am0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii-1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     Ampm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     Am0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    cs2 = T(0.125) *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_pp0) = T(0.5)*(cs1 + cs2);

    // csten(i,j,k,ist_p0p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk+2);
    ap(0,-1,0) =
                     Ap0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,-1);
    ap(1,-1,0) =
                     A00p(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Apmp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,-1);
    ap(1,1,0) =
                     A0mp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Ap00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(-1, 0, 0);
    ap(1,-1,1) =
                     A000(iii+1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Apm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,+1, 0);
    ap(1,1,1) =
                     A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0);
    cs1 = T(0.125) *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_m0p)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj  ,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj-1,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj  ,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj+1,kkk+2);

    ap(-1,-1,0) =
                     Am0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(+1, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     Ammp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Am00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     Amm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    cs2 = T(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_p0p) = T(0.5)*(cs1+cs2);

    // csten(i,j,k,ist_0pp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     A0pp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(1,0,0) =
                     Ampp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     A00p(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(1,1,0) =
                     Am0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     A0p0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(1,0,1) =
                     Amp0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     A000(iii-1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(1,1,1) =
                     Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0);
    cs1 = T(0.125) *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i,j+1,k,ist_0mp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj-2,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj-2,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     A0mp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(1,-1,0) =
                     Ammp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A0m0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(1,-1,1) =
                     Amm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    cs2 = T(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    csten(i,j,k,ist_0pp) = T(0.5)*(cs1+cs2);

    // csten(i,j,k,ist_ppp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk+2);
    p( 0, 0, 0) = T(1.);
    ap(0,0,0) =
                     Appp(iii,jjj,kkk) * p(-1,-1,-1);
    ap(1,0,0) =
                     A0pp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,-1);
    ap(0,1,0) =
                     Ap0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,-1);
    ap(1,1,0) =
                     A00p(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,0,1) =
                     App0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk+1) * p(-1,-1, 0);
    ap(1,0,1) =
                     A0p0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,1,1) =
                     Ap00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(-1, 0, 0);
    ap(1,1,1) =
                     A000(iii+1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0);
    cs1 = T(0.125) *
      ( restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mpp)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     Ampp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     Am0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     Amp0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    cs2 = T(0.125) *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    // alternative: csten(i,j+1,k,ist_pmp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj-2,kkk+2);
    p( 0, 0, 0) = T(1.);
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj-1,kkk+2);
    ap(0,-1,0) =
                     Apmp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,+1,-1);
    ap(1,-1,0) =
                     A0mp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Apm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,+1, 0);
    ap(1,-1,1) =
                     A0m0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    cs3 = T(0.125) *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    // alternative: csten(i+1,j+1,k,ist_mmp)
    iii = ii+2;
    jjj = jj+2;
    kkk = kk;
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj-2,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0, 0) = T(1.);
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj-2,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     Ammp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Amm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,+1, 0);
    cs4 = T(0.125) *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1));

    csten(i,j,k,ist_ppp) = T(0.25)*(cs1+cs2+cs3+cs4);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_sten_doit (int i, int j, int k, Array4<T const> const& x,
                              Array4<std::type_identity_t<T> const> const& sten) noexcept
{
    using namespace nodelap_detail;

    return     x(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000)
        //
        +      x(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00)
        +      x(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00)
        //
        +      x(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0)
        +      x(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0)
        //
        +      x(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p)
        +      x(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p)
        //
        +      x(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0)
        +      x(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0)
        +      x(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0)
        +      x(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0)
        //
        +      x(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p)
        +      x(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p)
        +      x(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p)
        +      x(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p)
        //
        +      x(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp)
        +      x(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp)
        +      x(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp)
        +      x(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp)
        //
        +      x(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp)
        +      x(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp)
        +      x(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp)
        +      x(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp)
        +      x(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp)
        +      x(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp)
        +      x(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp)
        +      x(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_adotx_sten (int i, int j, int k, Array4<T const> const& x,
                         Array4<std::type_identity_t<T> const> const& sten, Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        return T(0.0);
    } else {
        return mlndlap_adotx_sten_doit(i,j,k,x,sten);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (int i, int j, int k, Array4<T> const& sol,
                                Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else if (sten(i,j,k,ist_000) != T(0.0)) {
        T Ax = mlndlap_adotx_sten_doit<T>(i,j,k,sol,sten);
        sol(i,j,k) += (rhs(i,j,k) - Ax) / sten(i,j,k,ist_000);
    }
}

template <typename T>
inline
void mlndlap_gauss_seidel_sten (Box const& bx, Array4<T> const& sol,
                                Array4<std::type_identity_t<T> const> const& rhs,
                                Array4<std::type_identity_t<T> const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    AMREX_LOOP_3D(bx, i, j, k,
    {
        mlndlap_gauss_seidel_sten(i,j,k,sol,rhs,sten,msk);
    });
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int i, int j, int k, Array4<T> const& fine,
                            Array4<std::type_identity_t<T> const> const& crse, Array4<std::type_identity_t<T> const> const& sten,
                            Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    if (!msk(i,j,k) && sten(i,j,k,ist_000) != T(0.0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool ieven = ic*2 == i;
        bool jeven = jc*2 == j;
        bool keven = kc*2 == k;
        T fv;
        if (ieven && jeven && keven) {
            fv = crse(ic,jc,kc);
        } else if (ieven && jeven) {
            T w1 = std::abs(sten(i,j,k-1,ist_00p));
            T w2 = std::abs(sten(i,j,k  ,ist_00p));
            if (w1 == T(0.0) && w2 == T(0.0)) {
                fv = T(0.5)*(crse(ic,jc,kc)+crse(ic,jc,kc+1));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc,kc+1)) / (w1+w2);
            }
        } else if (ieven && keven) {
            T w1 = std::abs(sten(i,j-1,k,ist_0p0));
            T w2 = std::abs(sten(i,j  ,k,ist_0p0));
            if (w1 == T(0.0) && w2 == T(0.0)) {
                fv = T(0.5)*(crse(ic,jc,kc)+crse(ic,jc+1,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc+1,kc)) / (w1+w2);
            }
        } else if (jeven && keven) {
            T w1 = std::abs(sten(i-1,j,k,ist_p00));
            T w2 = std::abs(sten(i  ,j,k,ist_p00));
            if (w1 == T(0.0) && w2 == T(0.0)) {
                fv = T(0.5)*(crse(ic,jc,kc)+crse(ic+1,jc,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic+1,jc,kc)) / (w1+w2);
            }
        } else if (ieven) {
            T w1m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j-1,k  ,ist_0pp)) + eps_v<T>);
            T w1p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i,j  ,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps_v<T>);
            T w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k-1,ist_0pp)) + eps_v<T>);
            T w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i,j-1,k  ,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps_v<T>);
            T wmm = std::abs(sten(i,j-1,k-1,ist_0pp)) * (T(1.0) + w1m + w2m);
            T wpm = std::abs(sten(i,j  ,k-1,ist_0pp)) * (T(1.0) + w1p + w2m);
            T wmp = std::abs(sten(i,j-1,k  ,ist_0pp)) * (T(1.0) + w1m + w2p);
            T wpp = std::abs(sten(i,j  ,k  ,ist_0pp)) * (T(1.0) + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic,jc+1,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic,jc+1,kc+1))
                / (wmm+wpm+wmp+wpp+eps_v<T>);
        } else if (jeven) {
            T w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i-1,j,k  ,ist_p0p)) + eps_v<T>);
            T w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps_v<T>);
            T w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k-1,ist_p0p)) + eps_v<T>);
            T w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i-1,j,k  ,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps_v<T>);
            T wmm = std::abs(sten(i-1,j,k-1,ist_p0p)) * (T(1.0) + w1m + w2m);
            T wpm = std::abs(sten(i  ,j,k-1,ist_p0p)) * (T(1.0) + w1p + w2m);
            T wmp = std::abs(sten(i-1,j,k  ,ist_p0p)) * (T(1.0) + w1m + w2p);
            T wpp = std::abs(sten(i  ,j,k  ,ist_p0p)) * (T(1.0) + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic+1,jc,kc+1))
                / (wmm+wpm+wmp+wpp+eps_v<T>);
        } else if (keven) {
            T w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i-1,j  ,k,ist_pp0)) + eps_v<T>);
            T w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps_v<T>);
            T w2m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j-1,k,ist_pp0)) + eps_v<T>);
            T w2p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i-1,j  ,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps_v<T>);
            T wmm = std::abs(sten(i-1,j-1,k,ist_pp0)) * (T(1.0) + w1m + w2m);
            T wpm = std::abs(sten(i  ,j-1,k,ist_pp0)) * (T(1.0) + w1p + w2m);
            T wmp = std::abs(sten(i-1,j  ,k,ist_pp0)) * (T(1.0) + w1m + w2p);
            T wpp = std::abs(sten(i  ,j  ,k,ist_pp0)) * (T(1.0) + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc+1,kc) + wpp*crse(ic+1,jc+1,kc))
                / (wmm+wpm+wmp+wpp+eps_v<T>);
        } else {
            T wmmm = T(1.0);
            T wpmm = T(1.0);
            T wmpm = T(1.0);
            T wppm = T(1.0);
            T wmmp = T(1.0);
            T wpmp = T(1.0);
            T wmpp = T(1.0);
            T wppp = T(1.0);

            T wtmp = std::abs(sten(i-1,j,k,ist_p00)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i-1,j  ,k  ,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wmpm += wtmp;
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p00)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp))
                + std::abs(sten(i,j  ,k  ,ist_ppp)) + eps_v<T>);
            wpmm += wtmp;
            wppm += wtmp;
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0p0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i  ,j-1,k  ,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wpmm += wtmp;
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0p0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp))
                + std::abs(sten(i  ,j,k  ,ist_ppp)) + eps_v<T>);
            wmpm += wtmp;
            wppm += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i  ,j  ,k-1,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wpmm += wtmp;
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp))
                + std::abs(sten(i  ,j  ,k,ist_ppp)) + eps_v<T>);
            wmmp += wtmp;
            wpmp += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j-1,k,ist_pp0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wmmp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_pp0)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp)) + eps_v<T>);
            wpmm += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_pp0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp)) + eps_v<T>);
            wmpm += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_pp0)) /
                ( std::abs(sten(i,j,k-1,ist_ppp))
                + std::abs(sten(i,j,k  ,ist_ppp)) + eps_v<T>);
            wppm += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j,k-1,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wmpm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_p0p)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp)) + eps_v<T>);
            wpmm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp)) + eps_v<T>);
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p0p)) /
                ( std::abs(sten(i,j-1,k,ist_ppp))
                + std::abs(sten(i,j  ,k,ist_ppp)) + eps_v<T>);
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp)) + eps_v<T>);
            wmmm += wtmp;
            wpmm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp)) + eps_v<T>);
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp)) + eps_v<T>);
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0pp)) /
                ( std::abs(sten(i-1,j,k,ist_ppp))
                + std::abs(sten(i  ,j,k,ist_ppp)) + eps_v<T>);
            wmpp += wtmp;
            wppp += wtmp;

            wmmm *= std::abs(sten(i-1,j-1,k-1,ist_ppp));
            wpmm *= std::abs(sten(i  ,j-1,k-1,ist_ppp));
            wmpm *= std::abs(sten(i-1,j  ,k-1,ist_ppp));
            wppm *= std::abs(sten(i  ,j  ,k-1,ist_ppp));
            wmmp *= std::abs(sten(i-1,j-1,k  ,ist_ppp));
            wpmp *= std::abs(sten(i  ,j-1,k  ,ist_ppp));
            wmpp *= std::abs(sten(i-1,j  ,k  ,ist_ppp));
            wppp *= std::abs(sten(i  ,j  ,k  ,ist_ppp));
            fv = (wmmm*crse(ic,jc  ,kc  ) + wpmm*crse(ic+1,jc  ,kc  )
                  + wmpm*crse(ic,jc+1,kc  ) + wppm*crse(ic+1,jc+1,kc  )
                  + wmmp*crse(ic,jc  ,kc+1) + wpmp*crse(ic+1,jc  ,kc+1)
                  + wmpp*crse(ic,jc+1,kc+1) + wppp*crse(ic+1,jc+1,kc+1))
                / (wmmm + wpmm + wmpm + wppm + wmmp + wpmp + wmpp + wppp + eps_v<T>);
        }

        fine(i,j,k) += fv;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int i, int j, int k, Array4<T> const& crse,
                              Array4<std::type_identity_t<T> const> const& fine, Array4<std::type_identity_t<T> const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = T(0.0);
    } else {

        T cv = fine(ii,jj,kk);

        // ************************************
        // Adding fine(ii-1,jj,kk)
        // ************************************

        T sten_lo = std::abs(sten(ii-2,jj,kk,ist_p00));
        T sten_hi = std::abs(sten(ii-1,jj,kk,ist_p00));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii-1,jj,kk);
        } else {
            cv += fine(ii-1,jj,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii+1,jj,kk)
        // ************************************

        sten_lo = std::abs(sten(ii  ,jj,kk,ist_p00));
        sten_hi = std::abs(sten(ii+1,jj,kk,ist_p00));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii+1,jj,kk);
        } else {
            cv += fine(ii+1,jj,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj-1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj-2,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj-1,kk,ist_0p0));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii,jj-1,kk);
        } else {
            cv += fine(ii,jj-1,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj+1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj  ,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj+1,kk,ist_0p0));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii,jj+1,kk);
        } else {
            cv += fine(ii,jj+1,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk-1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk-2,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk-1,ist_00p));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii,jj,kk-1);
        } else {
            cv += fine(ii,jj,kk-1)*sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk+1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk  ,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk+1,ist_00p));

        if (sten_lo == T(0.0) && sten_hi == T(0.0)) {
            cv += T(0.5)*fine(ii,jj,kk+1);
        } else {
            cv += fine(ii,jj,kk+1)*sten_lo  / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii-1,jj-1,kk)
        // ************************************

        // keven
        T w1m = std::abs(sten(ii-2,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-2,jj-1,kk,ist_pp0)) + eps_v<T>);
        T w1p = std::abs(sten(ii-1,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-1,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps_v<T>);
        T w2m = std::abs(sten(ii-1,jj-2,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-2,kk,ist_pp0)) + eps_v<T>);
        T w2p = std::abs(sten(ii-1,jj-1,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-1,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps_v<T>);
        T wmm = std::abs(sten(ii-2,jj-2,kk,ist_pp0)) * (T(1.0) + w1m + w2m);
        T wpm = std::abs(sten(ii-1,jj-2,kk,ist_pp0)) * (T(1.0) + w1p + w2m);
        T wmp = std::abs(sten(ii-2,jj-1,kk,ist_pp0)) * (T(1.0) + w1m + w2p);
        T wpp = std::abs(sten(ii-1,jj-1,kk,ist_pp0)) * (T(1.0) + w1p + w2p);
        cv += fine(ii-1,jj-1,kk)*wpp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii+1,jj-1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj-1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii  ,jj-1,kk,ist_pp0)) + eps_v<T>);
        w1p = std::abs(sten(ii+1,jj-1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps_v<T>);
        w2m = std::abs(sten(ii+1,jj-2,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-2,kk,ist_pp0)) + eps_v<T>);
        w2p = std::abs(sten(ii+1,jj-1,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps_v<T>);
        wmm = std::abs(sten(ii  ,jj-2,kk,ist_pp0)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj-2,kk,ist_pp0)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj-1,kk,ist_pp0)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj-1,kk,ist_pp0)) * (T(1.0) + w1p + w2p);
        cv += fine(ii+1,jj-1,kk)*wmp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii-1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii-2,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-2,jj+1,kk,ist_pp0)) + eps_v<T>);
        w1p = std::abs(sten(ii-1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps_v<T>);
        w2m = std::abs(sten(ii-1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj  ,kk,ist_pp0)) + eps_v<T>);
        w2p = std::abs(sten(ii-1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj+1,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps_v<T>);
        wmm = std::abs(sten(ii-2,jj  ,kk,ist_pp0)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj  ,kk,ist_pp0)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj+1,kk,ist_pp0)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj+1,kk,ist_pp0)) * (T(1.0) + w1p + w2p);
        cv += fine(ii-1,jj+1,kk)*wpm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii+1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj+1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps_v<T>);
        w1p = std::abs(sten(ii+1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj  ,kk,ist_pp0)) + eps_v<T>);
        w2m = std::abs(sten(ii+1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps_v<T>);
        w2p = std::abs(sten(ii+1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj+1,kk,ist_pp0)) + eps_v<T>);
        wmm = std::abs(sten(ii  ,jj  ,kk,ist_pp0)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj  ,kk,ist_pp0)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj+1,kk,ist_pp0)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj+1,kk,ist_pp0)) * (T(1.0) + w1p + w2p);
        cv += fine(ii+1,jj+1,kk)*wmm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii-1,jj,kk-1)
        // ************************************

        // jeven
        w1m = std::abs(sten(ii-2,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-2,jj,kk-1,ist_p0p)) + eps_v<T>);
        w1p = std::abs(sten(ii-1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps_v<T>);
        w2m = std::abs(sten(ii-1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-2,ist_p0p)) + eps_v<T>);
        w2p = std::abs(sten(ii-1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps_v<T>);
        wmm = std::abs(sten(ii-2,jj,kk-2,ist_p0p)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk-2,ist_p0p)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk-1,ist_p0p)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk-1,ist_p0p)) * (T(1.0) + w1p + w2p);
        cv += fine(ii-1,jj,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii+1,jj,kk-1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk-1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps_v<T>);
        w1p = std::abs(sten(ii+1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii+1,jj,kk-1,ist_p0p)) + eps_v<T>);
        w2m = std::abs(sten(ii+1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-2,ist_p0p)) + eps_v<T>);
        w2p = std::abs(sten(ii+1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps_v<T>);
        wmm = std::abs(sten(ii  ,jj,kk-2,ist_p0p)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk-2,ist_p0p)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk-1,ist_p0p)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk-1,ist_p0p)) * (T(1.0) + w1p + w2p);
        cv += fine(ii+1,jj,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii-1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii-2,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-2,jj,kk  ,ist_p0p)) + eps_v<T>);
        w1p = std::abs(sten(ii-1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps_v<T>);
        w2m = std::abs(sten(ii-1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii-2,jj,kk  ,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps_v<T>);
        w2p = std::abs(sten(ii-1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk+1,ist_p0p)) + eps_v<T>);
        wmm = std::abs(sten(ii-2,jj,kk  ,ist_p0p)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk  ,ist_p0p)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk+1,ist_p0p)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk+1,ist_p0p)) * (T(1.0) + w1p + w2p);
        cv += fine(ii-1,jj,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii+1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk+1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps_v<T>);
        w1p = std::abs(sten(ii+1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
              +std::abs(sten(ii+1,jj,kk  ,ist_p0p)) + eps_v<T>);
        w2m = std::abs(sten(ii+1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii+1,jj,kk  ,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps_v<T>);
        w2p = std::abs(sten(ii+1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk+1,ist_p0p)) + eps_v<T>);
        wmm = std::abs(sten(ii  ,jj,kk  ,ist_p0p)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk  ,ist_p0p)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk+1,ist_p0p)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk+1,ist_p0p)) * (T(1.0) + w1p + w2p);
        cv += fine(ii+1,jj,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii,jj-1,kk-1)
        // ************************************

        // ieven
        w1m = std::abs(sten(ii,jj-2,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-2,kk-1,ist_0pp)) + eps_v<T>);
        w2m = std::abs(sten(ii,jj-1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-2,ist_0pp)) + eps_v<T>);
        w1p = std::abs(sten(ii,jj-1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps_v<T>);
        w2p = std::abs(sten(ii,jj-1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps_v<T>);
        wmm = std::abs(sten(ii,jj-2,kk-2,ist_0pp)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk-2,ist_0pp)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk-1,ist_0pp)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk-1,ist_0pp)) * (T(1.0) + w1p + w2p);
        cv += fine(ii,jj-1,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii,jj+1,kk-1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps_v<T>);
        w1p = std::abs(sten(ii,jj+1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj+1,kk-1,ist_0pp)) + eps_v<T>);
        w2m = std::abs(sten(ii,jj+1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-2,ist_0pp)) + eps_v<T>);
        w2p = std::abs(sten(ii,jj+1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps_v<T>);
        wmm = std::abs(sten(ii,jj  ,kk-2,ist_0pp)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk-2,ist_0pp)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk-1,ist_0pp)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk-1,ist_0pp)) * (T(1.0) + w1p + w2p);
        cv += fine(ii,jj+1,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii,jj-1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj-2,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-2,kk  ,ist_0pp)) + eps_v<T>);
        w1p = std::abs(sten(ii,jj-1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps_v<T>);
        w2m = std::abs(sten(ii,jj-1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj-2,kk  ,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps_v<T>);
        w2p = std::abs(sten(ii,jj-1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk+1,ist_0pp)) + eps_v<T>);
        wmm = std::abs(sten(ii,jj-2,kk  ,ist_0pp)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk  ,ist_0pp)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk+1,ist_0pp)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk+1,ist_0pp)) * (T(1.0) + w1p + w2p);
        cv += fine(ii,jj-1,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine(ii,jj+1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps_v<T>);
        w1p = std::abs(sten(ii,jj+1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj+1,kk  ,ist_0pp)) + eps_v<T>);
        w2m = std::abs(sten(ii,jj+1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj+1,kk  ,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps_v<T>);
        w2p = std::abs(sten(ii,jj+1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk+1,ist_0pp)) + eps_v<T>);
        wmm = std::abs(sten(ii,jj  ,kk  ,ist_0pp)) * (T(1.0) + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk  ,ist_0pp)) * (T(1.0) + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk+1,ist_0pp)) * (T(1.0) + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk+1,ist_0pp)) * (T(1.0) + w1p + w2p);
        cv += fine(ii,jj+1,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps_v<T>);

        // ************************************
        // Adding fine at corners
        // ************************************

        T wmmm = T(1.0)
            + std::abs(sten(ii  ,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp)) + eps_v<T>);
        wmmm *= std::abs(sten(ii,jj,kk,ist_ppp));
        cv += wmmm*fine(ii+1,jj+1,kk+1)*sten(ii+1,jj+1,kk+1,ist_inv);

        T wpmm = T(1.0)
            + std::abs(sten(ii-1,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp)) + eps_v<T>);
        wpmm *= std::abs(sten(ii-1,jj,kk,ist_ppp));
        cv += wpmm*fine(ii-1,jj+1,kk+1)*sten(ii-1,jj+1,kk+1,ist_inv);

        T wmpm = T(1.0)
            + std::abs(sten(ii  ,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps_v<T>);
        wmpm *= std::abs(sten(ii  ,jj-1,kk  ,ist_ppp));
        cv += wmpm*fine(ii+1,jj-1,kk+1)*sten(ii+1,jj-1,kk+1,ist_inv);

        T wppm = T(1.0)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps_v<T>);
        wppm *= std::abs(sten(ii-1,jj-1,kk  ,ist_ppp));
        cv += wppm*fine(ii-1,jj-1,kk+1)*sten(ii-1,jj-1,kk+1,ist_inv);

        T wmmp = T(1.0)
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps_v<T>);
        wmmp *= std::abs(sten(ii  ,jj  ,kk-1,ist_ppp));
        cv += wmmp*fine(ii+1,jj+1,kk-1)*sten(ii+1,jj+1,kk-1,ist_inv);

        T wpmp = T(1.0)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps_v<T>);
        wpmp *= std::abs(sten(ii-1,jj  ,kk-1,ist_ppp));
        cv += wpmp*fine(ii-1,jj+1,kk-1)*sten(ii-1,jj+1,kk-1,ist_inv);

        T wmpp = T(1.0)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps_v<T>);
        wmpp *= std::abs(sten(ii  ,jj-1,kk-1,ist_ppp));
        cv += wmpp*fine(ii+1,jj-1,kk-1)*sten(ii+1,jj-1,kk-1,ist_inv);

        T wppp = T(1.0)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps_v<T>);
        wppp *= std::abs(sten(ii-1,jj-1,kk-1,ist_ppp));
        cv += wppp*fine(ii-1,jj-1,kk-1)*sten(ii-1,jj-1,kk-1,ist_inv);

        crse(i,j,k) = cv * T(0.125);
    }
}

#ifdef AMREX_USE_EB

/// \cond DOXYGEN_IGNORE
namespace nodelap_detail {

    constexpr int i_S_x     = 0;
    constexpr int i_S_y     = 1;
    constexpr int i_S_z     = 2;
    constexpr int i_S_x2    = 3;
    constexpr int i_S_y2    = 4;
    constexpr int i_S_z2    = 5;
    constexpr int i_S_x_y   = 6;
    constexpr int i_S_x_z   = 7;
    constexpr int i_S_y_z   = 8;
    constexpr int i_S_x2_y  = 9;
    constexpr int i_S_x2_z  = 10;
    constexpr int i_S_x_y2  = 11;
    constexpr int i_S_y2_z  = 12;
    constexpr int i_S_x_z2  = 13;
    constexpr int i_S_y_z2  = 14;
    constexpr int i_S_x2_y2 = 15;
    constexpr int i_S_x2_z2 = 16;
    constexpr int i_S_y2_z2 = 17;
    constexpr int i_S_xyz   = 18;
    constexpr int n_Sintg   = 19;

    constexpr int i_c_xmym = 0;
    constexpr int i_c_xmyb = 1;
    constexpr int i_c_xmyp = 2;
    constexpr int i_c_xbym = 3;
    constexpr int i_c_xbyb = 4;
    constexpr int i_c_xbyp = 5;
    constexpr int i_c_xpym = 6;
    constexpr int i_c_xpyb = 7;
    constexpr int i_c_xpyp = 8;
    constexpr int i_c_xmzm = 9;
    constexpr int i_c_xmzb = 10;
    constexpr int i_c_xmzp = 11;
    constexpr int i_c_xbzm = 12;
    constexpr int i_c_xbzb = 13;
    constexpr int i_c_xbzp = 14;
    constexpr int i_c_xpzm = 15;
    constexpr int i_c_xpzb = 16;
    constexpr int i_c_xpzp = 17;
    constexpr int i_c_ymzm = 18;
    constexpr int i_c_ymzb = 19;
    constexpr int i_c_ymzp = 20;
    constexpr int i_c_ybzm = 21;
    constexpr int i_c_ybzb = 22;
    constexpr int i_c_ybzp = 23;
    constexpr int i_c_ypzm = 24;
    constexpr int i_c_ypzb = 25;
    constexpr int i_c_ypzp = 26;
    constexpr int n_conn = 27;

    constexpr int i_B_x        =  0;
    constexpr int i_B_y        =  1;
    constexpr int i_B_z        =  2;
    constexpr int i_B_x_y      =  3;
    constexpr int i_B_x_z      =  4;
    constexpr int i_B_y_z      =  5;
    constexpr int i_B_xyz      =  6;
    constexpr int numSurfIntgs =  7;

}
/// \endcond

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_connection (int i, int j, int k, Array4<T> const& conn,
                             Array4<std::type_identity_t<T> const> const& intg, Array4<std::type_identity_t<T> const> const& vol,
                             Array4<EBCellFlag const> const& flag) noexcept
{
    using namespace nodelap_detail;

    if (flag(i,j,k).isCovered()) {
        for (int n = 0; n < n_conn; ++n) { conn(i,j,k,n) = T(0.); }
    } else if (flag(i,j,k).isRegular() || vol(i,j,k) >= almostone_v<T>) {
        for (int n = 0; n < n_conn; ++n) { conn(i,j,k,n) = T(1.); }
    } else {
        // Scaled by 9
        conn(i,j,k,i_c_xmym) = T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(9.)*( intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xmyb) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + T(18.)*( intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmyp) =  T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(9.)*(-intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xbym) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(18.)*(intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 36
        conn(i,j,k,i_c_xbyb) = T(2.25)*vol(i,j,k)
            + T(9.)*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + T(36.)*intg(i,j,k,i_S_x2_y2);

        // Scaled by 18
        conn(i,j,k,i_c_xbyp) =  T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(18.)*(-intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpym) = T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(9.)*(-intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xpyb) = T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + T(18.)*(-intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpyp) =  T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + T(9.)*( intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                      +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzm) = T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_x) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(9.)*(intg(i,j,k,i_S_x_z) - intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xmzb) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + T(18.)*(intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzp) = T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_x  ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(9.)*(-intg(i,j,k,i_S_x_z  ) + intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzm) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(18.)*(intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzb) = T(2.25)*vol(i,j,k)
            + T(9.)*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + T(36.)*intg(i,j,k,i_S_x2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_xbzp) = T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(18.)*(-intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzm) = T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(9.)*(-intg(i,j,k,i_S_x_z ) - intg(i,j,k,i_S_x2_z)
                     +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xpzb) = T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_x   ) + intg(i,j,k,i_S_x2   ) - intg(i,j,k,i_S_z2))
            + T(18.)*(-intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzp) = T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + T(9.)*( intg(i,j,k,i_S_x_z ) + intg(i,j,k,i_S_x2_z)
                      +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzm) = T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(9.)*(intg(i,j,k,i_S_y_z) - intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ymzb) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_y) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + T(18.)*(intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzp) = T(0.5625)*vol(i,j,k)
            + T(2.25)*(-intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(9.)*(-intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ybzm) = T(1.125)*vol(i,j,k)
            + T(4.5)*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(18.)*(intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 36
        conn(i,j,k,i_c_ybzb) = T(2.25)*vol(i,j,k)
            + T(9.)*(-intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
                     + T(36.)*intg(i,j,k,i_S_y2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_ybzp) = T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(18.)*(-intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzm) = T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_y ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(9.)*(-intg(i,j,k,i_S_y_z ) - intg(i,j,k,i_S_y2_z)
                     +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ypzb) = T(1.125)*vol(i,j,k)
            + T(4.5)*( intg(i,j,k,i_S_y   ) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + T(18.)*(-intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzp) = T(0.5625)*vol(i,j,k)
            + T(2.25)*( intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + T(9.)*( intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                      +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_eb (int i, int j, int k, Array4<T> const& sten,
                             Array4<std::type_identity_t<T> const> const& sig, Array4<std::type_identity_t<T> const> const& conn,
                             GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    using namespace nodelap_detail;

    T facx = T(1./36.)*dxinv[0]*dxinv[0];
    T facy = T(1./36.)*dxinv[1]*dxinv[1];
    T facz = T(1./36.)*dxinv[2]*dxinv[2];

    // i+1,j,k
    sten(i,j,k,ist_p00) = (
        sig(i,j  ,k  )*(T(4.)*facx*conn(i,j  ,k  ,i_c_ymzm) - T(2.)*facy*conn(i,j  ,k  ,i_c_xbzm) - T(2.)*facz*conn(i,j  ,k  ,i_c_xbym) ) +
        sig(i,j-1,k  )*(T(4.)*facx*conn(i,j-1,k  ,i_c_ypzm) - T(2.)*facy*conn(i,j-1,k  ,i_c_xbzm) - T(2.)*facz*conn(i,j-1,k  ,i_c_xbyp) ) +
        sig(i,j  ,k-1)*(T(4.)*facx*conn(i,j  ,k-1,i_c_ymzp) - T(2.)*facy*conn(i,j  ,k-1,i_c_xbzp) - T(2.)*facz*conn(i,j  ,k-1,i_c_xbym) ) +
        sig(i,j-1,k-1)*(T(4.)*facx*conn(i,j-1,k-1,i_c_ypzp) - T(2.)*facy*conn(i,j-1,k-1,i_c_xbzp) - T(2.)*facz*conn(i,j-1,k-1,i_c_xbyp) ) );

    // i,j+1,k
    sten(i,j,k,ist_0p0) = (
        sig(i  ,j,k  )*(T(-2.)*facx*conn(i  ,j,k  ,i_c_ybzm) + T(4.)*facy*conn(i  ,j,k  ,i_c_xmzm) - T(2.)*facz*conn(i  ,j,k  ,i_c_xmyb) ) +
        sig(i-1,j,k  )*(T(-2.)*facx*conn(i-1,j,k  ,i_c_ybzm) + T(4.)*facy*conn(i-1,j,k  ,i_c_xpzm) - T(2.)*facz*conn(i-1,j,k  ,i_c_xpyb) ) +
        sig(i  ,j,k-1)*(T(-2.)*facx*conn(i  ,j,k-1,i_c_ybzp) + T(4.)*facy*conn(i  ,j,k-1,i_c_xmzp) - T(2.)*facz*conn(i  ,j,k-1,i_c_xmyb) ) +
        sig(i-1,j,k-1)*(T(-2.)*facx*conn(i-1,j,k-1,i_c_ybzp) + T(4.)*facy*conn(i-1,j,k-1,i_c_xpzp) - T(2.)*facz*conn(i-1,j,k-1,i_c_xpyb) ) );

    // i,j,k+1
    sten(i,j,k,ist_00p) = (
        sig(i  ,j  ,k)*(T(-2.)*facx*conn(i  ,j  ,k,i_c_ymzb) - T(2.)*facy*conn(i  ,j  ,k,i_c_xmzb) + T(4.)*facz*conn(i  ,j  ,k,i_c_xmym) ) +
        sig(i-1,j  ,k)*(T(-2.)*facx*conn(i-1,j  ,k,i_c_ymzb) - T(2.)*facy*conn(i-1,j  ,k,i_c_xpzb) + T(4.)*facz*conn(i-1,j  ,k,i_c_xpym) ) +
        sig(i  ,j-1,k)*(T(-2.)*facx*conn(i  ,j-1,k,i_c_ypzb) - T(2.)*facy*conn(i  ,j-1,k,i_c_xmzb) + T(4.)*facz*conn(i  ,j-1,k,i_c_xmyp) ) +
        sig(i-1,j-1,k)*(T(-2.)*facx*conn(i-1,j-1,k,i_c_ypzb) - T(2.)*facy*conn(i-1,j-1,k,i_c_xpzb) + T(4.)*facz*conn(i-1,j-1,k,i_c_xpyp) ) );

    // i+1,j+1,k
    sten(i,j,k,ist_pp0) = (
        sig(i,j,k  )*(T(2.)*facx*conn(i,j,k  ,i_c_ybzm) + T(2.)*facy*conn(i,j,k  ,i_c_xbzm) - facz*conn(i,j,k  ,i_c_xbyb) ) +
        sig(i,j,k-1)*(T(2.)*facx*conn(i,j,k-1,i_c_ybzp) + T(2.)*facy*conn(i,j,k-1,i_c_xbzp) - facz*conn(i,j,k-1,i_c_xbyb) ) );

    // i+1,j,k+1
    sten(i,j,k,ist_p0p) = (
        sig(i,j,k  )*(T(2.)*facx*conn(i,j,k  ,i_c_ymzb) - facy*conn(i,j,k  ,i_c_xbzb) + T(2.)*facz*conn(i,j,k  ,i_c_xbym) ) +
        sig(i,j-1,k)*(T(2.)*facx*conn(i,j-1,k,i_c_ypzb) - facy*conn(i,j-1,k,i_c_xbzb) + T(2.)*facz*conn(i,j-1,k,i_c_xbyp) ) );

    // i,j+1,k+1
    sten(i,j,k,ist_0pp) = (
        sig(i  ,j,k)*(-facx*conn(i  ,j,k,i_c_ybzb) + T(2.)*facy*conn(i  ,j,k,i_c_xmzb) + T(2.)*facz*conn(i  ,j,k,i_c_xmyb) ) +
        sig(i-1,j,k)*(-facx*conn(i-1,j,k,i_c_ybzb) + T(2.)*facy*conn(i-1,j,k,i_c_xpzb) + T(2.)*facz*conn(i-1,j,k,i_c_xpyb) ) );

    // i+1,j+1,k+1
    sten(i,j,k,ist_ppp) = sig(i,j,k) * (facx*conn(i,j,k,i_c_ybzb) + facy*conn(i,j,k,i_c_xbzb) + facz*conn(i,j,k,i_c_xbyb) );
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_eb (int i, int j, int k, Array4<T> const& rhs, Array4<std::type_identity_t<T> const> const& vel,
                      Array4<std::type_identity_t<T> const> const& vfrac, Array4<std::type_identity_t<T> const> const& intg,
                      Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                      Box const& nodal_domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    using namespace nodelap_detail;

    T facx = T(0.25)*dxinv[0];
    T facy = T(0.25)*dxinv[1];
    T facz = T(0.25)*dxinv[2];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (!msk(i,j,k)) {

        T zero_ilo = T(1.0);
        T zero_ihi = T(1.0);
        T zero_jlo = T(1.0);
        T zero_jhi = T(1.0);
        T zero_klo = T(1.0);
        T zero_khi = T(1.0);

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
            && i == domlo.x)
        {
            zero_ilo = T(0.0);
        }
        if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
            && i == domhi.x)
        {
            zero_ihi = T(0.0);
        }
        if ((bclo[1] == LinOpBCType::Neumann || bclo[1] == LinOpBCType::inflow)
            && j == domlo.y)
        {
            zero_jlo = T(0.0);
        }
        if ((bchi[1] == LinOpBCType::Neumann || bchi[1] == LinOpBCType::inflow)
            && j == domhi.y)
        {
            zero_jhi = T(0.0);
        }
        if ((bclo[2] == LinOpBCType::Neumann || bclo[2] == LinOpBCType::inflow)
            && k == domlo.z)
        {
            zero_klo = T(0.0);
        }
        if ((bchi[2] == LinOpBCType::Neumann || bchi[2] == LinOpBCType::inflow)
            && k == domhi.z)
        {
            zero_khi = T(0.0);
        }

        rhs(i,j,k) = facx*(
            vel(i-1,j-1,k  ,0)*(    -vfrac(i-1,j-1,k  )
                               -T(2.)*intg(i-1,j-1,k  ,i_S_y)
                               +T(2.)*intg(i-1,j-1,k  ,i_S_z)
                               +T(4.)*intg(i-1,j-1,k  ,i_S_y_z))*zero_jlo*zero_khi
           +vel(i  ,j-1,k  ,0)*(     vfrac(i  ,j-1,k  )
                               +T(2.)*intg(i  ,j-1,k  ,i_S_y)
                               -T(2.)*intg(i  ,j-1,k  ,i_S_z)
                               -T(4.)*intg(i  ,j-1,k  ,i_S_y_z))*zero_jlo*zero_khi
           +vel(i-1,j  ,k  ,0)*(    -vfrac(i-1,j  ,k  )
                               +T(2.)*intg(i-1,j  ,k  ,i_S_y)
                               +T(2.)*intg(i-1,j  ,k  ,i_S_z)
                               -T(4.)*intg(i-1,j  ,k  ,i_S_y_z))*zero_jhi*zero_khi
           +vel(i  ,j  ,k  ,0)*(     vfrac(i  ,j  ,k  )
                               -T(2.)*intg(i  ,j  ,k  ,i_S_y)
                               -T(2.)*intg(i  ,j  ,k  ,i_S_z)
                               +T(4.)*intg(i  ,j  ,k  ,i_S_y_z))*zero_jhi*zero_khi
           +vel(i-1,j-1,k-1,0)*(    -vfrac(i-1,j-1,k-1)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_y)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_z)
                               -T(4.)*intg(i-1,j-1,k-1,i_S_y_z))*zero_jlo*zero_klo
           +vel(i  ,j-1,k-1,0)*(     vfrac(i  ,j-1,k-1)
                               +T(2.)*intg(i  ,j-1,k-1,i_S_y)
                               +T(2.)*intg(i  ,j-1,k-1,i_S_z)
                               +T(4.)*intg(i  ,j-1,k-1,i_S_y_z))*zero_jlo*zero_klo
           +vel(i-1,j  ,k-1,0)*(    -vfrac(i-1,j  ,k-1)
                               +T(2.)*intg(i-1,j  ,k-1,i_S_y)
                               -T(2.)*intg(i-1,j  ,k-1,i_S_z)
                               +T(4.)*intg(i-1,j  ,k-1,i_S_y_z))*zero_jhi*zero_klo
           +vel(i  ,j  ,k-1,0)*(     vfrac(i  ,j  ,k-1)
                               -T(2.)*intg(i  ,j  ,k-1,i_S_y)
                               +T(2.)*intg(i  ,j  ,k-1,i_S_z)
                               -T(4.)*intg(i  ,j  ,k-1,i_S_y_z))*zero_jhi*zero_klo )
            + facy*(
            vel(i-1,j-1,k  ,1)*(    -vfrac(i-1,j-1,k  )
                               -T(2.)*intg(i-1,j-1,k  ,i_S_x)
                               +T(2.)*intg(i-1,j-1,k  ,i_S_z)
                               +T(4.)*intg(i-1,j-1,k  ,i_S_x_z))*zero_ilo*zero_khi
           +vel(i  ,j-1,k  ,1)*(    -vfrac(i  ,j-1,k  )
                               +T(2.)*intg(i  ,j-1,k  ,i_S_x)
                               +T(2.)*intg(i  ,j-1,k  ,i_S_z)
                               -T(4.)*intg(i  ,j-1,k  ,i_S_x_z))*zero_ihi*zero_khi
           +vel(i-1,j  ,k  ,1)*(     vfrac(i-1,j  ,k  )
                               +T(2.)*intg(i-1,j  ,k  ,i_S_x)
                               -T(2.)*intg(i-1,j  ,k  ,i_S_z)
                               -T(4.)*intg(i-1,j  ,k  ,i_S_x_z))*zero_ilo*zero_khi
           +vel(i  ,j  ,k  ,1)*(     vfrac(i  ,j  ,k  )
                               -T(2.)*intg(i  ,j  ,k  ,i_S_x)
                               -T(2.)*intg(i  ,j  ,k  ,i_S_z)
                               +T(4.)*intg(i  ,j  ,k  ,i_S_x_z))*zero_ihi*zero_khi
           +vel(i-1,j-1,k-1,1)*(    -vfrac(i-1,j-1,k-1)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_x)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_z)
                               -T(4.)*intg(i-1,j-1,k-1,i_S_x_z))*zero_ilo*zero_klo
           +vel(i  ,j-1,k-1,1)*(    -vfrac(i  ,j-1,k-1)
                               +T(2.)*intg(i  ,j-1,k-1,i_S_x)
                               -T(2.)*intg(i  ,j-1,k-1,i_S_z)
                               +T(4.)*intg(i  ,j-1,k-1,i_S_x_z))*zero_ihi*zero_klo
           +vel(i-1,j  ,k-1,1)*(     vfrac(i-1,j  ,k-1)
                               +T(2.)*intg(i-1,j  ,k-1,i_S_x)
                               +T(2.)*intg(i-1,j  ,k-1,i_S_z)
                               +T(4.)*intg(i-1,j  ,k-1,i_S_x_z))*zero_ilo*zero_klo
           +vel(i  ,j  ,k-1,1)*(     vfrac(i  ,j  ,k-1)
                               -T(2.)*intg(i  ,j  ,k-1,i_S_x)
                               +T(2.)*intg(i  ,j  ,k-1,i_S_z)
                               -T(4.)*intg(i  ,j  ,k-1,i_S_x_z))*zero_ihi*zero_klo )
            + facz*(
            vel(i-1,j-1,k  ,2)*(     vfrac(i-1,j-1,k  )
                               +T(2.)*intg(i-1,j-1,k  ,i_S_x)
                               +T(2.)*intg(i-1,j-1,k  ,i_S_y)
                               +T(4.)*intg(i-1,j-1,k  ,i_S_x_y))*zero_ilo*zero_jlo
           +vel(i  ,j-1,k  ,2)*(     vfrac(i  ,j-1,k  )
                               -T(2.)*intg(i  ,j-1,k  ,i_S_x)
                               +T(2.)*intg(i  ,j-1,k  ,i_S_y)
                               -T(4.)*intg(i  ,j-1,k  ,i_S_x_y))*zero_ihi*zero_jlo
           +vel(i-1,j  ,k  ,2)*(     vfrac(i-1,j  ,k  )
                               +T(2.)*intg(i-1,j  ,k  ,i_S_x)
                               -T(2.)*intg(i-1,j  ,k  ,i_S_y)
                               -T(4.)*intg(i-1,j  ,k  ,i_S_x_y))*zero_ilo*zero_jhi
           +vel(i  ,j  ,k  ,2)*(     vfrac(i  ,j  ,k  )
                               -T(2.)*intg(i  ,j  ,k  ,i_S_x)
                               -T(2.)*intg(i  ,j  ,k  ,i_S_y)
                               +T(4.)*intg(i  ,j  ,k  ,i_S_x_y))*zero_ihi*zero_jhi
           +vel(i-1,j-1,k-1,2)*(    -vfrac(i-1,j-1,k-1)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_x)
                               -T(2.)*intg(i-1,j-1,k-1,i_S_y)
                               -T(4.)*intg(i-1,j-1,k-1,i_S_x_y))*zero_ilo*zero_jlo
           +vel(i  ,j-1,k-1,2)*(    -vfrac(i  ,j-1,k-1)
                               +T(2.)*intg(i  ,j-1,k-1,i_S_x)
                               -T(2.)*intg(i  ,j-1,k-1,i_S_y)
                               +T(4.)*intg(i  ,j-1,k-1,i_S_x_y))*zero_ihi*zero_jlo
           +vel(i-1,j  ,k-1,2)*(    -vfrac(i-1,j  ,k-1)
                               -T(2.)*intg(i-1,j  ,k-1,i_S_x)
                               +T(2.)*intg(i-1,j  ,k-1,i_S_y)
                               +T(4.)*intg(i-1,j  ,k-1,i_S_x_y))*zero_ilo*zero_jhi
           +vel(i  ,j  ,k-1,2)*(    -vfrac(i  ,j  ,k-1)
                               +T(2.)*intg(i  ,j  ,k-1,i_S_x)
                               +T(2.)*intg(i  ,j  ,k-1,i_S_y)
                               -T(4.)*intg(i  ,j  ,k-1,i_S_x_y))*zero_ihi*zero_jhi );
    } else {
        rhs(i,j,k) = T(0.);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void add_eb_flow_contribution (int i, int j, int k, Array4<T> const& rhs,
                      Array4<int const> const& msk, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv,
                      Array4<std::type_identity_t<T> const> const& bareaarr,
                      Array4<std::type_identity_t<T> const> const& sintg,
                      Array4<std::type_identity_t<T> const> const& eb_vel_dot_n) noexcept
{
    using namespace nodelap_detail;

    T fac_eb = T(0.125) * dxinv[0];

    if (!msk(i,j,k)) {
        rhs(i,j,k) += fac_eb*(
            eb_vel_dot_n(i-1,j-1,k  )*(        bareaarr(i-1,j-1,k)
                                        +T(2.)*sintg(i-1,j-1,k  ,i_B_x)
                                        +T(2.)*sintg(i-1,j-1,k  ,i_B_y)
                                        +T(4.)*sintg(i-1,j-1,k  ,i_B_x_y)
                                        -T(2.)*sintg(i-1,j-1,k  ,i_B_z)
                                        -T(4.)*sintg(i-1,j-1,k  ,i_B_y_z)
                                        -T(4.)*sintg(i-1,j-1,k  ,i_B_x_z)
                                        -T(8.)*sintg(i-1,j-1,k  ,i_B_xyz))
           +eb_vel_dot_n(i  ,j-1,k  )*(        bareaarr(i  ,j-1,k)
                                        -T(2.)*sintg(i  ,j-1,k  ,i_B_x)
                                        +T(2.)*sintg(i  ,j-1,k  ,i_B_y)
                                        -T(4.)*sintg(i  ,j-1,k  ,i_B_x_y)
                                        -T(2.)*sintg(i  ,j-1,k  ,i_B_z)
                                        -T(4.)*sintg(i  ,j-1,k  ,i_B_y_z)
                                        +T(4.)*sintg(i  ,j-1,k  ,i_B_x_z)
                                        +T(8.)*sintg(i  ,j-1,k  ,i_B_xyz))
           +eb_vel_dot_n(i-1,j  ,k  )*(        bareaarr(i-1,j  ,k)
                                        +T(2.)*sintg(i-1,j  ,k  ,i_B_x)
                                        -T(2.)*sintg(i-1,j  ,k  ,i_B_y)
                                        -T(4.)*sintg(i-1,j  ,k  ,i_B_x_y)
                                        -T(2.)*sintg(i-1,j  ,k  ,i_B_z)
                                        +T(4.)*sintg(i-1,j  ,k  ,i_B_y_z)
                                        -T(4.)*sintg(i-1,j  ,k  ,i_B_x_z)
                                        +T(8.)*sintg(i-1,j  ,k  ,i_B_xyz))
           +eb_vel_dot_n(i  ,j  ,k  )*(        bareaarr(i  ,j  ,k)
                                        -T(2.)*sintg(i  ,j  ,k  ,i_B_x)
                                        -T(2.)*sintg(i  ,j  ,k  ,i_B_y)
                                        +T(4.)*sintg(i  ,j  ,k  ,i_B_x_y)
                                        -T(2.)*sintg(i  ,j  ,k  ,i_B_z)
                                        +T(4.)*sintg(i  ,j  ,k  ,i_B_y_z)
                                        +T(4.)*sintg(i  ,j  ,k  ,i_B_x_z)
                                        -T(8.)*sintg(i  ,j  ,k  ,i_B_xyz))
           +eb_vel_dot_n(i-1,j-1,k-1)*(        bareaarr(i-1,j-1,k-1)
                                        +T(2.)*sintg(i-1,j-1,k-1,i_B_x)
                                        +T(2.)*sintg(i-1,j-1,k-1,i_B_y)
                                        +T(4.)*sintg(i-1,j-1,k-1,i_B_x_y)
                                        +T(2.)*sintg(i-1,j-1,k-1,i_B_z)
                                        +T(4.)*sintg(i-1,j-1,k-1,i_B_y_z)
                                        +T(4.)*sintg(i-1,j-1,k-1,i_B_x_z)
                                        +T(8.)*sintg(i-1,j-1,k-1,i_B_xyz))
           +eb_vel_dot_n(i  ,j-1,k-1)*(        bareaarr(i  ,j-1,k-1)
                                        -T(2.)*sintg(i  ,j-1,k-1,i_B_x)
                                        +T(2.)*sintg(i  ,j-1,k-1,i_B_y)
                                        -T(4.)*sintg(i  ,j-1,k-1,i_B_x_y)
                                        +T(2.)*sintg(i  ,j-1,k-1,i_B_z)
                                        +T(4.)*sintg(i  ,j-1,k-1,i_B_y_z)
                                        -T(4.)*sintg(i  ,j-1,k-1,i_B_x_z)
                                        -T(8.)*sintg(i  ,j-1,k-1,i_B_xyz))
           +eb_vel_dot_n(i-1,j  ,k-1)*(        bareaarr(i-1,j  ,k-1)
                                        +T(2.)*sintg(i-1,j  ,k-1,i_B_x)
                                        -T(2.)*sintg(i-1,j  ,k-1,i_B_y)
                                        -T(4.)*sintg(i-1,j  ,k-1,i_B_x_y)
                                        +T(2.)*sintg(i-1,j  ,k-1,i_B_z)
                                        -T(4.)*sintg(i-1,j  ,k-1,i_B_y_z)
                                        +T(4.)*sintg(i-1,j  ,k-1,i_B_x_z)
                                        -T(8.)*sintg(i-1,j  ,k-1,i_B_xyz))
           +eb_vel_dot_n(i  ,j  ,k-1)*(        bareaarr(i  ,j  ,k-1)
                                        -T(2.)*sintg(i  ,j  ,k-1,i_B_x)
                                        -T(2.)*sintg(i  ,j  ,k-1,i_B_y)
                                        +T(4.)*sintg(i  ,j  ,k-1,i_B_x_y)
                                        +T(2.)*sintg(i  ,j  ,k-1,i_B_z)
                                        -T(4.)*sintg(i  ,j  ,k-1,i_B_y_z)
                                        -T(4.)*sintg(i  ,j  ,k-1,i_B_x_z)
                                        +T(8.)*sintg(i  ,j  ,k-1,i_B_xyz)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                        Array4<std::type_identity_t<T> const> const& sig, Array4<std::type_identity_t<T> const> const& vfrac,
                        Array4<std::type_identity_t<T> const> const& intg, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    using namespace nodelap_detail;

    if (vfrac(i,j,k) == T(0.)) {
        u(i,j,k,0) = u(i,j,k,1) = u(i,j,k,2) = T(0.);
    } else {
        T dpdx = T(0.25)*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
        T dpdy = T(0.25)*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
        T dpdz = T(0.25)*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
                                +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));

        T dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k);

        u(i,j,k,0) -= sig(i,j,k)*dxinv[0]*(dpdx + T(0.5)*intg(i,j,k,i_S_y  )*dpp_xy +
                                                  T(0.5)*intg(i,j,k,i_S_z  )*dpp_xz +
                                                            intg(i,j,k,i_S_y_z)*dpp_xyz );
        u(i,j,k,1) -= sig(i,j,k)*dxinv[1]*(dpdy + T(0.5)*intg(i,j,k,i_S_x  )*dpp_xy +
                                                  T(0.5)*intg(i,j,k,i_S_z  )*dpp_yz +
                                                            intg(i,j,k,i_S_x_z)*dpp_xyz );
        u(i,j,k,2) -= sig(i,j,k)*dxinv[2]*(dpdz + T(0.5)*intg(i,j,k,i_S_x  )*dpp_xz +
                                                  T(0.5)*intg(i,j,k,i_S_y  )*dpp_yz +
                                                            intg(i,j,k,i_S_x_y)*dpp_xyz );
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb_c (int i, int j, int k, Array4<T> const& u, Array4<std::type_identity_t<T> const> const& p,
                        std::type_identity_t<T> sig, Array4<std::type_identity_t<T> const> const& vfrac,
                        Array4<std::type_identity_t<T> const> const& intg, GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    using namespace nodelap_detail;

    if (vfrac(i,j,k) == T(0.)) {
        u(i,j,k,0) = u(i,j,k,1) = u(i,j,k,2) = T(0.);
    } else {
        T dpdx = T(0.25)*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
        T dpdy = T(0.25)*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
                                -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
        T dpdz = T(0.25)*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
                                +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));

        T dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        T dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k);

        u(i,j,k,0) -= sig*dxinv[0]*(dpdx + T(0.5)*intg(i,j,k,i_S_y  )*dpp_xy +
                                           T(0.5)*intg(i,j,k,i_S_z  )*dpp_xz +
                                                     intg(i,j,k,i_S_y_z)*dpp_xyz );
        u(i,j,k,1) -= sig*dxinv[1]*(dpdy + T(0.5)*intg(i,j,k,i_S_x  )*dpp_xy +
                                           T(0.5)*intg(i,j,k,i_S_z  )*dpp_yz +
                                                     intg(i,j,k,i_S_x_z)*dpp_xyz );
        u(i,j,k,2) -= sig*dxinv[2]*(dpdz + T(0.5)*intg(i,j,k,i_S_x  )*dpp_xz +
                                           T(0.5)*intg(i,j,k,i_S_y  )*dpp_yz +
                                                     intg(i,j,k,i_S_x_y)*dpp_xyz );
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlndlap_rhcc_eb (int i, int j, int k, Array4<T const> const& rhcc,
                      Array4<std::type_identity_t<T> const> const& vfrac, Array4<std::type_identity_t<T> const> const& intg,
                      Array4<int const> const& msk) noexcept
{
    using namespace nodelap_detail;

    if (!msk(i,j,k)) {
        return
                          rhcc(i  ,j  ,k  ) *
         ( T(0.125) * vfrac(i  ,j  ,k  )
         + T(0.25) * (-intg(i  ,j  ,k  ,i_S_x)
                         -intg(i  ,j  ,k  ,i_S_y)
                         -intg(i  ,j  ,k  ,i_S_z))
         + T(0.5) * (  intg(i  ,j  ,k  ,i_S_x_y)
                         +intg(i  ,j  ,k  ,i_S_x_z)
                         +intg(i  ,j  ,k  ,i_S_y_z))
            +          ( -intg(i  ,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k  ) *
         ( T(0.125) * vfrac(i-1,j  ,k  )
         + T(0.25) * ( intg(i-1,j  ,k  ,i_S_x)
                         -intg(i-1,j  ,k  ,i_S_y)
                         -intg(i-1,j  ,k  ,i_S_z))
         + T(0.5) * ( -intg(i-1,j  ,k  ,i_S_x_y)
                         -intg(i-1,j  ,k  ,i_S_x_z)
                         +intg(i-1,j  ,k  ,i_S_y_z))
            +          (  intg(i-1,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k  ) *
         ( T(0.125) * vfrac(i  ,j-1,k  )
         + T(0.25) * (-intg(i  ,j-1,k  ,i_S_x)
                         +intg(i  ,j-1,k  ,i_S_y)
                         -intg(i  ,j-1,k  ,i_S_z))
         + T(0.5) * ( -intg(i  ,j-1,k  ,i_S_x_y)
                         +intg(i  ,j-1,k  ,i_S_x_z)
                         -intg(i  ,j-1,k  ,i_S_y_z))
            +          (  intg(i  ,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k  ) *
         ( T(0.125) * vfrac(i-1,j-1,k  )
         + T(0.25) * ( intg(i-1,j-1,k  ,i_S_x)
                         +intg(i-1,j-1,k  ,i_S_y)
                         -intg(i-1,j-1,k  ,i_S_z))
         + T(0.5) * (  intg(i-1,j-1,k  ,i_S_x_y)
                         -intg(i-1,j-1,k  ,i_S_x_z)
                         -intg(i-1,j-1,k  ,i_S_y_z))
            +          ( -intg(i-1,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j  ,k-1) *
         ( T(0.125) * vfrac(i  ,j  ,k-1)
         + T(0.25) * (-intg(i  ,j  ,k-1,i_S_x)
                         -intg(i  ,j  ,k-1,i_S_y)
                         +intg(i  ,j  ,k-1,i_S_z))
         + T(0.5) * (  intg(i  ,j  ,k-1,i_S_x_y)
                         -intg(i  ,j  ,k-1,i_S_x_z)
                         -intg(i  ,j  ,k-1,i_S_y_z))
            +          (  intg(i  ,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k-1) *
         ( T(0.125) * vfrac(i-1,j  ,k-1)
         + T(0.25) * ( intg(i-1,j  ,k-1,i_S_x)
                         -intg(i-1,j  ,k-1,i_S_y)
                         +intg(i-1,j  ,k-1,i_S_z))
         + T(0.5) * ( -intg(i-1,j  ,k-1,i_S_x_y)
                         +intg(i-1,j  ,k-1,i_S_x_z)
                         -intg(i-1,j  ,k-1,i_S_y_z))
            +          ( -intg(i-1,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k-1) *
         ( T(0.125) * vfrac(i  ,j-1,k-1)
         + T(0.25) * (-intg(i  ,j-1,k-1,i_S_x)
                         +intg(i  ,j-1,k-1,i_S_y)
                         +intg(i  ,j-1,k-1,i_S_z))
         + T(0.5) * ( -intg(i  ,j-1,k-1,i_S_x_y)
                         -intg(i  ,j-1,k-1,i_S_x_z)
                         +intg(i  ,j-1,k-1,i_S_y_z))
            +          ( -intg(i  ,j-1,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k-1) *
         ( T(0.125) * vfrac(i-1,j-1,k-1)
         + T(0.25) * ( intg(i-1,j-1,k-1,i_S_x)
                         +intg(i-1,j-1,k-1,i_S_y)
                         +intg(i-1,j-1,k-1,i_S_z))
         + T(0.5) * (  intg(i-1,j-1,k-1,i_S_x_y)
                         +intg(i-1,j-1,k-1,i_S_x_z)
                         +intg(i-1,j-1,k-1,i_S_y_z))
             +         (  intg(i-1,j-1,k-1,i_S_xyz)));
    } else {
        return T(0.);
    }
}

#endif


template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_sten_cpu (Box const& ndbx,
                                 Array4<AlgGid const> const& gid,
                                 Array4<int const> const& lid,
                                 AlgInt* ncols, AlgInt* cols,
                                 Real* mat, // NOLINT(readability-non-const-parameter)
                                 Array4<Real const> const& sten) noexcept
{
    constexpr int ist_000 = 0;
    constexpr int ist_p00 = 1;
    constexpr int ist_0p0 = 2;
    constexpr int ist_00p = 3;
    constexpr int ist_pp0 = 4;
    constexpr int ist_p0p = 5;
    constexpr int ist_0pp = 6;
    constexpr int ist_ppp = 7;

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            cols[nelems] = gid(i,j,k);
            mat[nelems] = sten(i,j,k,ist_000);
            AlgInt nelems_old = nelems;
            ++nelems;

            if                (gid(i-1,j-1,k-1) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k-1);
                mat[nelems] = sten(i-1,j-1,k-1,ist_ppp);
                ++nelems;
            }

            if                (gid(i,j-1,k-1) < gidmax) {
                cols[nelems] = gid(i,j-1,k-1);
                mat[nelems] = sten(i,j-1,k-1,ist_0pp);
                ++nelems;
            }

            if                (gid(i+1,j-1,k-1) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k-1);
                mat[nelems] = sten(i,j-1,k-1,ist_ppp);
                ++nelems;
            }

            if                (gid(i-1,j,k-1) < gidmax) {
                cols[nelems] = gid(i-1,j,k-1);
                mat[nelems] = sten(i-1,j,k-1,ist_p0p);
                ++nelems;
            }

            if                (gid(i,j,k-1) < gidmax) {
                cols[nelems] = gid(i,j,k-1);
                mat[nelems] = sten(i,j,k-1,ist_00p);
                ++nelems;
            }

            if                (gid(i+1,j,k-1) < gidmax) {
                cols[nelems] = gid(i+1,j,k-1);
                mat[nelems] = sten(i,j,k-1,ist_p0p);
                ++nelems;
            }

            if                (gid(i-1,j+1,k-1) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k-1);
                mat[nelems] = sten(i-1,j,k-1,ist_ppp);
                ++nelems;
            }

            if                (gid(i,j+1,k-1) < gidmax) {
                cols[nelems] = gid(i,j+1,k-1);
                mat[nelems] = sten(i,j,k-1,ist_0pp);
                ++nelems;
            }

            if                (gid(i+1,j+1,k-1) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k-1);
                mat[nelems] = sten(i,j,k-1,ist_ppp);
                ++nelems;
            }

            if                (gid(i-1,j-1,k) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k);
                mat[nelems] = sten(i-1,j-1,k,ist_pp0);
                ++nelems;
            }

            if                (gid(i,j-1,k) < gidmax) {
                cols[nelems] = gid(i,j-1,k);
                mat[nelems] = sten(i,j-1,k,ist_0p0);
                ++nelems;
            }

            if                (gid(i+1,j-1,k) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k);
                mat[nelems] = sten(i,j-1,k,ist_pp0);
                ++nelems;
            }

            if                (gid(i-1,j,k) < gidmax) {
                cols[nelems] = gid(i-1,j,k);
                mat[nelems] = sten(i-1,j,k,ist_p00);
                ++nelems;
            }

            if                (gid(i+1,j,k) < gidmax) {
                cols[nelems] = gid(i+1,j,k);
                mat[nelems] = sten(i,j,k,ist_p00);
                ++nelems;
            }

            if                (gid(i-1,j+1,k) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k);
                mat[nelems] = sten(i-1,j,k,ist_pp0);
                ++nelems;
            }

            if                (gid(i,j+1,k) < gidmax) {
                cols[nelems] = gid(i,j+1,k);
                mat[nelems] = sten(i,j,k,ist_0p0);
                ++nelems;
            }

            if                (gid(i+1,j+1,k) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k);
                mat[nelems] = sten(i,j,k,ist_pp0);
                ++nelems;
            }

            if                (gid(i-1,j-1,k+1) < gidmax) {
                cols[nelems] = gid(i-1,j-1,k+1);
                mat[nelems] = sten(i-1,j-1,k,ist_ppp);
                ++nelems;
            }

            if                (gid(i,j-1,k+1) < gidmax) {
                cols[nelems] = gid(i,j-1,k+1);
                mat[nelems] = sten(i,j-1,k,ist_0pp);
                ++nelems;
            }

            if                (gid(i+1,j-1,k+1) < gidmax) {
                cols[nelems] = gid(i+1,j-1,k+1);
                mat[nelems] = sten(i,j-1,k,ist_ppp);
                ++nelems;
            }

            if                (gid(i-1,j,k+1) < gidmax) {
                cols[nelems] = gid(i-1,j,k+1);
                mat[nelems] = sten(i-1,j,k,ist_p0p);
                ++nelems;
            }

            if                (gid(i,j,k+1) < gidmax) {
                cols[nelems] = gid(i,j,k+1);
                mat[nelems] = sten(i,j,k,ist_00p);
                ++nelems;
            }

            if                (gid(i+1,j,k+1) < gidmax) {
                cols[nelems] = gid(i+1,j,k+1);
                mat[nelems] = sten(i,j,k,ist_p0p);
                ++nelems;
            }

            if                (gid(i-1,j+1,k+1) < gidmax) {
                cols[nelems] = gid(i-1,j+1,k+1);
                mat[nelems] = sten(i-1,j,k,ist_ppp);
                ++nelems;
            }

            if                (gid(i,j+1,k+1) < gidmax) {
                cols[nelems] = gid(i,j+1,k+1);
                mat[nelems] = sten(i,j,k,ist_0pp);
                ++nelems;
            }

            if                (gid(i+1,j+1,k+1) < gidmax) {
                cols[nelems] = gid(i+1,j+1,k+1);
                mat[nelems] = sten(i,j,k,ist_ppp);
                ++nelems;
            }

            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_aa_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sig,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
    Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
    Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
    Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
    Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
    Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;

    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = sig(i-1,j-1,k-1) * fxyz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = sig(i  ,j-1,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k-1) && fm2xm2y4z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp +=     sig(i,j,k-1) * fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax) {
                    cols[nelems] = gid(i,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=     sig(i  ,j,k-1) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = sig(i-1,j  ,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=     sig(i,j  ,k-1) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = sig(i  ,j  ,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k) && fm2x4ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp +=     sig(i,j-1,k) * fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=     sig(i  ,j-1,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k) && f4xm2ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp +=     sig(i-1,j,k) * f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k) && f4xm2ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=     sig(i  ,j,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=     sig(i  ,j-1,k) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp +=     sig(i  ,j,k) * f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=     sig(i-1,j  ,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k) && fm2x4ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=     sig(i,j  ,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=     sig(i-1,j  ,k) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp +=     sig(i,j  ,k) * fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp +=     sig(i  ,j  ,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp +=     sig(i  ,j  ,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = sig(i-1,j-1,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=     sig(i,j-1,k  ) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = sig(i  ,j-1,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=     sig(i-1,j,k  ) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k+1) && fm2xm2y4z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=     sig(i,j-1,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=     sig(i-1,j,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp +=     sig(i,j,k  ) * fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax) {
                    cols[nelems] = gid(i,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp +=     sig(i  ,j-1,k  ) * f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp +=     sig(i  ,j,k  ) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = sig(i-1,j  ,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp +=     sig(i-1,j  ,k  ) * fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp +=     sig(i,j  ,k  ) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = sig(i  ,j  ,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_ha_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sx,
                               Array4<Real const> const& sy,
                               Array4<Real const> const& sz,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];

    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = sx(i-1,j-1,k-1) * facx
                    +      sy(i-1,j-1,k-1) * facy
                    +      sz(i-1,j-1,k-1) * facz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx
                        +       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        +       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx
                        +       sy(i,j-1,k-1) * facy * Real(2.0)
                        +       sz(i,j-1,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = sx(i  ,j-1,k-1) * facx
                    +      sy(i  ,j-1,k-1) * facy
                    +      sz(i  ,j-1,k-1) * facz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(2.0)
                        -       sy(i-1,j-1,k-1) * facy
                        +       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=      sx(i-1,j,k-1) * facx * Real(2.0)
                        -       sy(i-1,j,k-1) * facy
                        +       sz(i-1,j,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx * Real(2.0)
                        -       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        +       sz(i-1,j-1,k-1) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx * Real(2.0)
                        -       sy(i,j-1,k-1) * facy * Real(2.0)
                        +       sz(i,j-1,k-1) * facz * Real(4.0);

                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=    - sx(i-1,j,k-1) * facx * Real(2.0)
                        -       sy(i-1,j,k-1) * facy * Real(2.0)
                        +       sz(i-1,j,k-1) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp +=    - sx(i,j,k-1) * facx * Real(2.0)
                        -       sy(i,j,k-1) * facy * Real(2.0)
                        +       sz(i,j,k-1) * facz * Real(4.0);
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(2.0)
                        -       sy(i  ,j-1,k-1) * facy
                        +       sz(i  ,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=      sx(i  ,j,k-1) * facx * Real(2.0)
                        -       sy(i  ,j,k-1) * facy
                        +       sz(i  ,j,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = sx(i-1,j  ,k-1) * facx
                    +      sy(i-1,j  ,k-1) * facy
                    +      sz(i-1,j  ,k-1) * facz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=    - sx(i-1,j  ,k-1) * facx
                        +       sy(i-1,j  ,k-1) * facy * Real(2.0)
                        +       sz(i-1,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=    - sx(i,j  ,k-1) * facx
                        +       sy(i,j  ,k-1) * facy * Real(2.0)
                        +       sz(i,j  ,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = sx(i  ,j  ,k-1) * facx
                    +      sy(i  ,j  ,k-1) * facy
                    +      sz(i  ,j  ,k-1) * facz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(2.0)
                        +       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        -       sz(i-1,j-1,k-1) * facz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=      sx(i-1,j-1,k) * facx * Real(2.0)
                        +       sy(i-1,j-1,k) * facy * Real(2.0)
                        -       sz(i-1,j-1,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx * Real(2.0)
                        +       sy(i-1,j-1,k-1) * facy * Real(4.0)
                        -       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx * Real(2.0)
                        +       sy(i,j-1,k-1) * facy * Real(4.0)
                        -       sz(i,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=    - sx(i-1,j-1,k) * facx * Real(2.0)
                        +       sy(i-1,j-1,k) * facy * Real(4.0)
                        -       sz(i-1,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp +=    - sx(i,j-1,k) * facx * Real(2.0)
                        +       sy(i,j-1,k) * facy * Real(4.0)
                        -       sz(i,j-1,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(2.0)
                        +       sy(i  ,j-1,k-1) * facy * Real(2.0)
                        -       sz(i  ,j-1,k-1) * facz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=      sx(i  ,j-1,k) * facx * Real(2.0)
                        +       sy(i  ,j-1,k) * facy * Real(2.0)
                        -       sz(i  ,j-1,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(4.0)
                        -       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        -       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=      sx(i-1,j,k-1) * facx * Real(4.0)
                        -       sy(i-1,j,k-1) * facy * Real(2.0)
                        -       sz(i-1,j,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=      sx(i-1,j-1,k) * facx * Real(4.0)
                        -       sy(i-1,j-1,k) * facy * Real(2.0)
                        -       sz(i-1,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp +=      sx(i-1,j,k) * facx * Real(4.0)
                        -       sy(i-1,j,k) * facy * Real(2.0)
                        -       sz(i-1,j,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(4.0)
                        -       sy(i  ,j-1,k-1) * facy * Real(2.0)
                        -       sz(i  ,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=      sx(i  ,j,k-1) * facx * Real(4.0)
                        -       sy(i  ,j,k-1) * facy * Real(2.0)
                        -       sz(i  ,j,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=      sx(i  ,j-1,k) * facx * Real(4.0)
                        -       sy(i  ,j-1,k) * facy * Real(2.0)
                        -       sz(i  ,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp +=      sx(i  ,j,k) * facx * Real(4.0)
                        -       sy(i  ,j,k) * facy * Real(2.0)
                        -       sz(i  ,j,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=      sx(i-1,j  ,k-1) * facx * Real(2.0)
                        +       sy(i-1,j  ,k-1) * facy * Real(2.0)
                        -       sz(i-1,j  ,k-1) * facz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=      sx(i-1,j  ,k) * facx * Real(2.0)
                        +       sy(i-1,j  ,k) * facy * Real(2.0)
                        -       sz(i-1,j  ,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=    - sx(i-1,j  ,k-1) * facx * Real(2.0)
                        +       sy(i-1,j  ,k-1) * facy * Real(4.0)
                        -       sz(i-1,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=    - sx(i,j  ,k-1) * facx * Real(2.0)
                        +       sy(i,j  ,k-1) * facy * Real(4.0)
                        -       sz(i,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=    - sx(i-1,j  ,k) * facx * Real(2.0)
                        +       sy(i-1,j  ,k) * facy * Real(4.0)
                        -       sz(i-1,j  ,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp +=    - sx(i,j  ,k) * facx * Real(2.0)
                        +       sy(i,j  ,k) * facy * Real(4.0)
                        -       sz(i,j  ,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp +=      sx(i  ,j  ,k-1) * facx * Real(2.0)
                        +       sy(i  ,j  ,k-1) * facy * Real(2.0)
                        -       sz(i  ,j  ,k-1) * facz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp +=      sx(i  ,j  ,k) * facx * Real(2.0)
                        +       sy(i  ,j  ,k) * facy * Real(2.0)
                        -       sz(i  ,j  ,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = sx(i-1,j-1,k  ) * facx
                    +      sy(i-1,j-1,k  ) * facy
                    +      sz(i-1,j-1,k  ) * facz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=    - sx(i-1,j-1,k  ) * facx
                        +       sy(i-1,j-1,k  ) * facy * Real(2.0)
                        +       sz(i-1,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=    - sx(i,j-1,k  ) * facx
                        +       sy(i,j-1,k  ) * facy * Real(2.0)
                        +       sz(i,j-1,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = sx(i  ,j-1,k  ) * facx
                    +      sy(i  ,j-1,k  ) * facy
                    +      sz(i  ,j-1,k  ) * facz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=      sx(i-1,j-1,k  ) * facx * Real(2.0)
                        -       sy(i-1,j-1,k  ) * facy
                        +       sz(i-1,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=      sx(i-1,j,k  ) * facx * Real(2.0)
                        -       sy(i-1,j,k  ) * facy
                        +       sz(i-1,j,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=    - sx(i-1,j-1,k  ) * facx * Real(2.0)
                        -       sy(i-1,j-1,k  ) * facy * Real(2.0)
                        +       sz(i-1,j-1,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=    - sx(i,j-1,k  ) * facx * Real(2.0)
                        -       sy(i,j-1,k  ) * facy * Real(2.0)
                        +       sz(i,j-1,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=    - sx(i-1,j,k  ) * facx * Real(2.0)
                        -       sy(i-1,j,k  ) * facy * Real(2.0)
                        +       sz(i-1,j,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp +=    - sx(i,j,k  ) * facx * Real(2.0)
                        -       sy(i,j,k  ) * facy * Real(2.0)
                        +       sz(i,j,k  ) * facz * Real(4.0);
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax && tmp != Real(0.0)) {
                    cols[nelems] = gid(i,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp +=      sx(i  ,j-1,k  ) * facx * Real(2.0)
                        -       sy(i  ,j-1,k  ) * facy
                        +       sz(i  ,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp +=      sx(i  ,j,k  ) * facx * Real(2.0)
                        -       sy(i  ,j,k  ) * facy
                        +       sz(i  ,j,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = sx(i-1,j  ,k  ) * facx
                    +      sy(i-1,j  ,k  ) * facy
                    +      sz(i-1,j  ,k  ) * facz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp +=    - sx(i-1,j  ,k  ) * facx
                        +       sy(i-1,j  ,k  ) * facy * Real(2.0)
                        +       sz(i-1,j  ,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp +=    - sx(i,j  ,k  ) * facx
                        +       sy(i,j  ,k  ) * facy * Real(2.0)
                        +       sz(i,j  ,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = sx(i  ,j  ,k  ) * facx
                    +      sy(i  ,j  ,k  ) * facy
                    +      sz(i  ,j  ,k  ) * facz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

template <typename AlgInt, typename AlgGid>
void mlndlap_fillijmat_cs_cpu (Box const& ndbx,
                               Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Real sigma,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0] * sigma;
    Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1] * sigma;
    Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2] * sigma;
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
    Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
    Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
    Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
    Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
    Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;

    const Box& nddom = amrex::surroundingNodes(ccdom);

    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    AlgInt nelems = 0;
    amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
    {
        if (lid(i,j,k) >= 0)
        {
            AlgInt nelems_old = nelems;
            cols[nelems_old] = gid(i,j,k);
            Real m0 = Real(0.);
            ++nelems;

            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k-1) && fm2xm2y4z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp += fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax) {
                    cols[nelems] = gid(i,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k-1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k) && fm2x4ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp += fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    cols[nelems] = gid(i,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k) && f4xm2ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp += f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    cols[nelems] = gid(i-1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k) && f4xm2ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp += f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    cols[nelems] = gid(i+1,j,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k) && fm2x4ym2z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp += fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    cols[nelems] = gid(i,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j-1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j,k+1) && fm2xm2y4z != Real(0.0)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp += fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax) {
                    cols[nelems] = gid(i,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i-1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    cols[nelems] = gid(i+1,j+1,k+1);
                    mat[nelems] = tmp;
                    ++nelems;
                }
            }

            mat[nelems_old] = m0;
            ncols[lid(i,j,k)] = nelems - nelems_old;
        }
    });
}

#ifdef AMREX_USE_GPU

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_sten_gpu (const int ps, const int i, const int j, const int k,
                                 const int offset,
                                 Array4<AlgGid const> const& gid,
                                 Array4<int const> const& lid,
                                 AlgInt* ncols, AlgInt* cols,
                                 Real* mat, // NOLINT(readability-non-const-parameter)
                                 Array4<Real const> const& sten) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        constexpr int ist_000 = 1-1;
        constexpr int ist_p00 = 2-1;
        constexpr int ist_0p0 = 3-1;
        constexpr int ist_00p = 4-1;
        constexpr int ist_pp0 = 5-1;
        constexpr int ist_p0p = 6-1;
        constexpr int ist_0pp = 7-1;
        constexpr int ist_ppp = 8-1;

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;

        if (offset == 1 || offset == 0) {
            if                (gid(i-1,j-1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j-1,k-1);
                    mat[ps] = sten(i-1,j-1,k-1,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if                (gid(i,j-1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j-1,k-1);
                    mat[ps] = sten(i,j-1,k-1,ist_0pp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if                (gid(i+1,j-1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j-1,k-1);
                    mat[ps] = sten(i,j-1,k-1,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if                (gid(i-1,j,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j,k-1);
                    mat[ps] = sten(i-1,j,k-1,ist_p0p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if                (gid(i,j,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j,k-1);
                    mat[ps] = sten(i,j,k-1,ist_00p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if                (gid(i+1,j,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j,k-1);
                    mat[ps] = sten(i,j,k-1,ist_p0p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if                (gid(i-1,j+1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j+1,k-1);
                    mat[ps] = sten(i-1,j,k-1,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if                (gid(i,j+1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j+1,k-1);
                    mat[ps] = sten(i,j,k-1,ist_0pp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 9 || offset == 0) {
            if                (gid(i+1,j+1,k-1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j+1,k-1);
                    mat[ps] = sten(i,j,k-1,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 10 || offset == 0) {
            if                (gid(i-1,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j-1,k);
                    mat[ps] = sten(i-1,j-1,k,ist_pp0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 11 || offset == 0) {
            if                (gid(i,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j-1,k);
                    mat[ps] = sten(i,j-1,k,ist_0p0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 12 || offset == 0) {
            if                (gid(i+1,j-1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j-1,k);
                    mat[ps] = sten(i,j-1,k,ist_pp0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 13 || offset == 0) {
            if                (gid(i-1,j,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j,k);
                    mat[ps] = sten(i-1,j,k,ist_p00);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 14 || offset == 0) {
            if                (gid(i+1,j,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j,k);
                    mat[ps] = sten(i,j,k,ist_p00);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 15 || offset == 0) {
            if                (gid(i-1,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j+1,k);
                    mat[ps] = sten(i-1,j,k,ist_pp0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 16 || offset == 0) {
            if                (gid(i,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j+1,k);
                    mat[ps] = sten(i,j,k,ist_0p0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 17 || offset == 0) {
            if                (gid(i+1,j+1,k) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j+1,k);
                    mat[ps] = sten(i,j,k,ist_pp0);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 18 || offset == 0) {
            if                (gid(i-1,j-1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j-1,k+1);
                    mat[ps] = sten(i-1,j-1,k,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 19 || offset == 0) {
            if                (gid(i,j-1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j-1,k+1);
                    mat[ps] = sten(i,j-1,k,ist_0pp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 20 || offset == 0) {
            if                (gid(i+1,j-1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j-1,k+1);
                    mat[ps] = sten(i,j-1,k,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 21 || offset == 0) {
            if                (gid(i-1,j,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j,k+1);
                    mat[ps] = sten(i-1,j,k,ist_p0p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 22 || offset == 0) {
            if                (gid(i,j,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j,k+1);
                    mat[ps] = sten(i,j,k,ist_00p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 23 || offset == 0) {
            if                (gid(i+1,j,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j,k+1);
                    mat[ps] = sten(i,j,k,ist_p0p);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 24 || offset == 0) {
            if                (gid(i-1,j+1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i-1,j+1,k+1);
                    mat[ps] = sten(i-1,j,k,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 25 || offset == 0) {
            if                (gid(i,j+1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i,j+1,k+1);
                    mat[ps] = sten(i,j,k,ist_0pp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        if (offset == 26 || offset == 0) {
            if                (gid(i+1,j+1,k+1) < gidmax) {
                if (offset != 0) {
                    cols[ps] = gid(i+1,j+1,k+1);
                    mat[ps] = sten(i,j,k,ist_ppp);
                }
                ++nelems;
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = sten(i,j,k,ist_000);
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_aa_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sig,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
        Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
        Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
        Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
        Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
        Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
        Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
        Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;

        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = sig(i-1,j-1,k-1) * fxyz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = sig(i  ,j-1,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp +=     sig(i,j,k-1) * fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=     sig(i  ,j,k-1) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = sig(i-1,j  ,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=     sig(i,j  ,k-1) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 9 || offset == 0) {
            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = sig(i  ,j  ,k-1) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 10 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 11 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=     sig(i,j-1,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp +=     sig(i,j-1,k) * fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 12 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=     sig(i  ,j-1,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 13 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=     sig(i-1,j-1,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=     sig(i-1,j,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=     sig(i-1,j-1,k) * f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp +=     sig(i-1,j,k) * f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 14 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=     sig(i  ,j-1,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=     sig(i  ,j,k-1) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=     sig(i  ,j-1,k) * f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp +=     sig(i  ,j,k) * f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 15 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=     sig(i-1,j  ,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 16 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=     sig(i-1,j  ,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=     sig(i,j  ,k-1) * fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=     sig(i-1,j  ,k) * fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp +=     sig(i,j  ,k) * fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 17 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp +=     sig(i  ,j  ,k-1) * f2x2ymz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp +=     sig(i  ,j  ,k) * f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 18 || offset == 0) {
            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = sig(i-1,j-1,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 19 || offset == 0) {
            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=     sig(i,j-1,k  ) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 20 || offset == 0) {
            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = sig(i  ,j-1,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 21 || offset == 0) {
            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=     sig(i-1,j,k  ) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 22 || offset == 0) {
            if (nddom.contains(i,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=     sig(i-1,j-1,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=     sig(i,j-1,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=     sig(i-1,j,k  ) * fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp +=     sig(i,j,k  ) * fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 23 || offset == 0) {
            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp +=     sig(i  ,j-1,k  ) * f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp +=     sig(i  ,j,k  ) * f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 24 || offset == 0) {
            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = sig(i-1,j  ,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 25 || offset == 0) {
            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp +=     sig(i-1,j  ,k  ) * fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp +=     sig(i,j  ,k  ) * fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 26 || offset == 0) {
            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = sig(i  ,j  ,k  ) * fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_ha_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Array4<Real const> const& sx,
                               Array4<Real const> const& sy,
                               Array4<Real const> const& sz,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
        Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
        Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];

        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = sx(i-1,j-1,k-1) * facx
                    +      sy(i-1,j-1,k-1) * facy
                    +      sz(i-1,j-1,k-1) * facz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx
                        +       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        +       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx
                        +       sy(i,j-1,k-1) * facy * Real(2.0)
                        +       sz(i,j-1,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = sx(i  ,j-1,k-1) * facx
                    +      sy(i  ,j-1,k-1) * facy
                    +      sz(i  ,j-1,k-1) * facz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(2.0)
                        -       sy(i-1,j-1,k-1) * facy
                        +       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=      sx(i-1,j,k-1) * facx * Real(2.0)
                        -       sy(i-1,j,k-1) * facy
                        +       sz(i-1,j,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx * Real(2.0)
                        -       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        +       sz(i-1,j-1,k-1) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx * Real(2.0)
                        -       sy(i,j-1,k-1) * facy * Real(2.0)
                        +       sz(i,j-1,k-1) * facz * Real(4.0);

                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=    - sx(i-1,j,k-1) * facx * Real(2.0)
                        -       sy(i-1,j,k-1) * facy * Real(2.0)
                        +       sz(i-1,j,k-1) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp +=    - sx(i,j,k-1) * facx * Real(2.0)
                        -       sy(i,j,k-1) * facy * Real(2.0)
                        +       sz(i,j,k-1) * facz * Real(4.0);
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(2.0)
                        -       sy(i  ,j-1,k-1) * facy
                        +       sz(i  ,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=      sx(i  ,j,k-1) * facx * Real(2.0)
                        -       sy(i  ,j,k-1) * facy
                        +       sz(i  ,j,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = sx(i-1,j  ,k-1) * facx
                    +      sy(i-1,j  ,k-1) * facy
                    +      sz(i-1,j  ,k-1) * facz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=    - sx(i-1,j  ,k-1) * facx
                        +       sy(i-1,j  ,k-1) * facy * Real(2.0)
                        +       sz(i-1,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=    - sx(i,j  ,k-1) * facx
                        +       sy(i,j  ,k-1) * facy * Real(2.0)
                        +       sz(i,j  ,k-1) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 9 || offset == 0) {
            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = sx(i  ,j  ,k-1) * facx
                    +      sy(i  ,j  ,k-1) * facy
                    +      sz(i  ,j  ,k-1) * facz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 10 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(2.0)
                        +       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        -       sz(i-1,j-1,k-1) * facz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=      sx(i-1,j-1,k) * facx * Real(2.0)
                        +       sy(i-1,j-1,k) * facy * Real(2.0)
                        -       sz(i-1,j-1,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 11 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=    - sx(i-1,j-1,k-1) * facx * Real(2.0)
                        +       sy(i-1,j-1,k-1) * facy * Real(4.0)
                        -       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp +=    - sx(i,j-1,k-1) * facx * Real(2.0)
                        +       sy(i,j-1,k-1) * facy * Real(4.0)
                        -       sz(i,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=    - sx(i-1,j-1,k) * facx * Real(2.0)
                        +       sy(i-1,j-1,k) * facy * Real(4.0)
                        -       sz(i-1,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp +=    - sx(i,j-1,k) * facx * Real(2.0)
                        +       sy(i,j-1,k) * facy * Real(4.0)
                        -       sz(i,j-1,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 12 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(2.0)
                        +       sy(i  ,j-1,k-1) * facy * Real(2.0)
                        -       sz(i  ,j-1,k-1) * facz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=      sx(i  ,j-1,k) * facx * Real(2.0)
                        +       sy(i  ,j-1,k) * facy * Real(2.0)
                        -       sz(i  ,j-1,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 13 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp +=      sx(i-1,j-1,k-1) * facx * Real(4.0)
                        -       sy(i-1,j-1,k-1) * facy * Real(2.0)
                        -       sz(i-1,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp +=      sx(i-1,j,k-1) * facx * Real(4.0)
                        -       sy(i-1,j,k-1) * facy * Real(2.0)
                        -       sz(i-1,j,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp +=      sx(i-1,j-1,k) * facx * Real(4.0)
                        -       sy(i-1,j-1,k) * facy * Real(2.0)
                        -       sz(i-1,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp +=      sx(i-1,j,k) * facx * Real(4.0)
                        -       sy(i-1,j,k) * facy * Real(2.0)
                        -       sz(i-1,j,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 14 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp +=      sx(i  ,j-1,k-1) * facx * Real(4.0)
                        -       sy(i  ,j-1,k-1) * facy * Real(2.0)
                        -       sz(i  ,j-1,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp +=      sx(i  ,j,k-1) * facx * Real(4.0)
                        -       sy(i  ,j,k-1) * facy * Real(2.0)
                        -       sz(i  ,j,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp +=      sx(i  ,j-1,k) * facx * Real(4.0)
                        -       sy(i  ,j-1,k) * facy * Real(2.0)
                        -       sz(i  ,j-1,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp +=      sx(i  ,j,k) * facx * Real(4.0)
                        -       sy(i  ,j,k) * facy * Real(2.0)
                        -       sz(i  ,j,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 15 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=      sx(i-1,j  ,k-1) * facx * Real(2.0)
                        +       sy(i-1,j  ,k-1) * facy * Real(2.0)
                        -       sz(i-1,j  ,k-1) * facz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=      sx(i-1,j  ,k) * facx * Real(2.0)
                        +       sy(i-1,j  ,k) * facy * Real(2.0)
                        -       sz(i-1,j  ,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 16 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp +=    - sx(i-1,j  ,k-1) * facx * Real(2.0)
                        +       sy(i-1,j  ,k-1) * facy * Real(4.0)
                        -       sz(i-1,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp +=    - sx(i,j  ,k-1) * facx * Real(2.0)
                        +       sy(i,j  ,k-1) * facy * Real(4.0)
                        -       sz(i,j  ,k-1) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp +=    - sx(i-1,j  ,k) * facx * Real(2.0)
                        +       sy(i-1,j  ,k) * facy * Real(4.0)
                        -       sz(i-1,j  ,k) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp +=    - sx(i,j  ,k) * facx * Real(2.0)
                        +       sy(i,j  ,k) * facy * Real(4.0)
                        -       sz(i,j  ,k) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 17 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp +=      sx(i  ,j  ,k-1) * facx * Real(2.0)
                        +       sy(i  ,j  ,k-1) * facy * Real(2.0)
                        -       sz(i  ,j  ,k-1) * facz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp +=      sx(i  ,j  ,k) * facx * Real(2.0)
                        +       sy(i  ,j  ,k) * facy * Real(2.0)
                        -       sz(i  ,j  ,k) * facz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 18 || offset == 0) {
            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = sx(i-1,j-1,k  ) * facx
                    +      sy(i-1,j-1,k  ) * facy
                    +      sz(i-1,j-1,k  ) * facz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 19 || offset == 0) {
            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=    - sx(i-1,j-1,k  ) * facx
                        +       sy(i-1,j-1,k  ) * facy * Real(2.0)
                        +       sz(i-1,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=    - sx(i,j-1,k  ) * facx
                        +       sy(i,j-1,k  ) * facy * Real(2.0)
                        +       sz(i,j-1,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 20 || offset == 0) {
            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = sx(i  ,j-1,k  ) * facx
                    +      sy(i  ,j-1,k  ) * facy
                    +      sz(i  ,j-1,k  ) * facz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 21 || offset == 0) {
            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=      sx(i-1,j-1,k  ) * facx * Real(2.0)
                        -       sy(i-1,j-1,k  ) * facy
                        +       sz(i-1,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=      sx(i-1,j,k  ) * facx * Real(2.0)
                        -       sy(i-1,j,k  ) * facy
                        +       sz(i-1,j,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 22 || offset == 0) {
            if (nddom.contains(i,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp +=    - sx(i-1,j-1,k  ) * facx * Real(2.0)
                        -       sy(i-1,j-1,k  ) * facy * Real(2.0)
                        +       sz(i-1,j-1,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp +=    - sx(i,j-1,k  ) * facx * Real(2.0)
                        -       sy(i,j-1,k  ) * facy * Real(2.0)
                        +       sz(i,j-1,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp +=    - sx(i-1,j,k  ) * facx * Real(2.0)
                        -       sy(i-1,j,k  ) * facy * Real(2.0)
                        +       sz(i-1,j,k  ) * facz * Real(4.0);
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp +=    - sx(i,j,k  ) * facx * Real(2.0)
                        -       sy(i,j,k  ) * facy * Real(2.0)
                        +       sz(i,j,k  ) * facz * Real(4.0);
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 23 || offset == 0) {
            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp +=      sx(i  ,j-1,k  ) * facx * Real(2.0)
                        -       sy(i  ,j-1,k  ) * facy
                        +       sz(i  ,j-1,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp +=      sx(i  ,j,k  ) * facx * Real(2.0)
                        -       sy(i  ,j,k  ) * facy
                        +       sz(i  ,j,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 24 || offset == 0) {
            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = sx(i-1,j  ,k  ) * facx
                    +      sy(i-1,j  ,k  ) * facy
                    +      sz(i-1,j  ,k  ) * facz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 25 || offset == 0) {
            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp +=    - sx(i-1,j  ,k  ) * facx
                        +       sy(i-1,j  ,k  ) * facy * Real(2.0)
                        +       sz(i-1,j  ,k  ) * facz * Real(2.0);
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp +=    - sx(i,j  ,k  ) * facx
                        +       sy(i,j  ,k  ) * facy * Real(2.0)
                        +       sz(i,j  ,k  ) * facz * Real(2.0);
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 26 || offset == 0) {
            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = sx(i  ,j  ,k  ) * facx
                    +      sy(i  ,j  ,k  ) * facy
                    +      sz(i  ,j  ,k  ) * facz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

template <typename AlgInt, typename AlgGid>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_fillijmat_cs_gpu (const int ps, const int i, const int j, const int k,
                               const int offset,
                               Box const& /*ndbx*/, Array4<AlgGid const> const& gid,
                               Array4<int const> const& lid,
                               AlgInt* ncols, AlgInt* cols,
                               Real* mat, // NOLINT(readability-non-const-parameter)
                               Real sigma, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                               Box const& ccdom) noexcept
{
    if (lid(i,j,k) >= 0)
    {
        Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0] * sigma;
        Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1] * sigma;
        Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2] * sigma;
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
        Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
        Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
        Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
        Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
        Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;

        const Box& nddom = amrex::surroundingNodes(ccdom);

        constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
        int nelems = 0;
        Real m0 = Real(0.);

        if (offset == 1 || offset == 0) {
            if (nddom.contains(i-1,j-1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if (               gid(i-1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 2 || offset == 0) {
            if (nddom.contains(i,j-1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 3 || offset == 0) {
            if (nddom.contains(i+1,j-1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 4 || offset == 0) {
            if (nddom.contains(i-1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 5 || offset == 0) {
            if (nddom.contains(i,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k-1)) {
                    tmp += fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 6 || offset == 0) {
            if (nddom.contains(i+1,j,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 7 || offset == 0) {
            if (nddom.contains(i-1,j+1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 8 || offset == 0) {
            if (nddom.contains(i,j+1,k-1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 9 || offset == 0) {
            if (nddom.contains(i+1,j+1,k-1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k-1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k-1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 10 || offset == 0) {
            if (nddom.contains(i-1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 11 || offset == 0) {
            if (nddom.contains(i,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j-1,k)) {
                    tmp += fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 12 || offset == 0) {
            if (nddom.contains(i+1,j-1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j-1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 13 || offset == 0) {
            if (nddom.contains(i-1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j-1,k)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i-1,j,k)) {
                    tmp += f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 14 || offset == 0) {
            if (nddom.contains(i+1,j,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k-1)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j-1,k)) {
                    tmp += f4xm2ym2z;
                }
                if (ccdom.contains(i  ,j,k)) {
                    tmp += f4xm2ym2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 15 || offset == 0) {
            if (nddom.contains(i-1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i-1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 16 || offset == 0) {
            if (nddom.contains(i,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k-1)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i-1,j  ,k)) {
                    tmp += fm2x4ym2z;
                }
                if (ccdom.contains(i,j  ,k)) {
                    tmp += fm2x4ym2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 17 || offset == 0) {
            if (nddom.contains(i+1,j+1,k)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j  ,k-1)) {
                    tmp += f2x2ymz;
                }
                if (ccdom.contains(i  ,j  ,k)) {
                    tmp += f2x2ymz;
                }
                m0 -= tmp;
                if                (gid(i+1,j+1,k) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 18 || offset == 0) {
            if (nddom.contains(i-1,j-1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 19 || offset == 0) {
            if (nddom.contains(i,j-1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 20 || offset == 0) {
            if (nddom.contains(i+1,j-1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j-1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j-1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 21 || offset == 0) {
            if (nddom.contains(i-1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i-1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 22 || offset == 0) {
            if (nddom.contains(i,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j-1,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j-1,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i-1,j,k  )) {
                    tmp += fm2xm2y4z;
                }
                if (ccdom.contains(i,j,k  )) {
                    tmp += fm2xm2y4z;
                }
                m0 -= tmp;
                if                (gid(i,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 23 || offset == 0) {
            if (nddom.contains(i+1,j,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i  ,j-1,k  )) {
                    tmp += f2xmy2z;
                }
                if (ccdom.contains(i  ,j,k  )) {
                    tmp += f2xmy2z;
                }
                m0 -= tmp;
                if                (gid(i+1,j,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 24 || offset == 0) {
            if (nddom.contains(i-1,j+1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i-1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i-1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 25 || offset == 0) {
            if (nddom.contains(i,j+1,k+1)) {
                Real tmp = Real(0.);
                if (ccdom.contains(i-1,j  ,k  )) {
                    tmp += fmx2y2z;
                }
                if (ccdom.contains(i,j  ,k  )) {
                    tmp += fmx2y2z;
                }
                m0 -= tmp;
                if                (gid(i,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        if (offset == 26 || offset == 0) {
            if (nddom.contains(i+1,j+1,k+1)) {
                Real tmp = fxyz;
                m0 -= tmp;
                if                (gid(i+1,j+1,k+1) < gidmax) {
                    if (offset != 0) {
                        cols[ps] = gid(i+1,j+1,k+1);
                        mat[ps] = tmp;
                    }
                    ++nelems;
                }
            }
            if (offset != 0) { return; }
        }

        // Only offset == 0 could get this far.
        cols[ps] = gid(i,j,k);
        mat[ps] = m0;
        ncols[lid(i,j,k)] = nelems+1;
    }
}

#endif


template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_ha (int i, int j, int k, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sx,
                         Array4<std::type_identity_t<T> const> const& sy, Array4<std::type_identity_t<T> const> const& sz,
                         Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
        T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];

        T s0 = T(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                   +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                             +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                   +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                             +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                   +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        T Ax = sol(i,j,k)*s0
                 + sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                    +facy*sy(i-1,j-1,k-1)
                                    +facz*sz(i-1,j-1,k-1))
                 + sol(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                    +facy*sy(i  ,j-1,k-1)
                                    +facz*sz(i  ,j-1,k-1))
                 + sol(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                    +facy*sy(i-1,j  ,k-1)
                                    +facz*sz(i-1,j  ,k-1))
                 + sol(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                    +facy*sy(i  ,j  ,k-1)
                                    +facz*sz(i  ,j  ,k-1))
                 + sol(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                    +facy*sy(i-1,j-1,k  )
                                    +facz*sz(i-1,j-1,k  ))
                 + sol(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                    +facy*sy(i  ,j-1,k  )
                                    +facz*sz(i  ,j-1,k  ))
                 + sol(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                    +facy*sy(i-1,j  ,k  )
                                    +facz*sz(i-1,j  ,k  ))
                 + sol(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                    +facy*sy(i  ,j  ,k  )
                                    +facz*sz(i  ,j  ,k  ))
                 +sol(i  ,j-1,k-1)*(          -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                    +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                    +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                 +sol(i  ,j+1,k-1)*(          -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                    +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                    +T(2.0)*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                 +sol(i  ,j-1,k+1)*(          -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                    +T(2.0)*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                    +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                 +sol(i  ,j+1,k+1)*(          -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                    +T(2.0)*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                    +T(2.0)*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                 +sol(i-1,j  ,k-1)*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                              -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                    +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                 +sol(i+1,j  ,k-1)*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                              -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                    +T(2.0)*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                 +sol(i-1,j  ,k+1)*( T(2.0)*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                              -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                    +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                 +sol(i+1,j  ,k+1)*( T(2.0)*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                              -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                    +T(2.0)*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                 +sol(i-1,j-1,k  )*( T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                    +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                              -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                 +sol(i+1,j-1,k  )*( T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                    +T(2.0)*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                              -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                 +sol(i-1,j+1,k  )*( T(2.0)*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                    +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                              -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                 +sol(i+1,j+1,k  )*( T(2.0)*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                    +T(2.0)*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                              -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                 + T(2.0)*sol(i-1,j,k)*(T(2.0)*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                                    -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                                    -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                 + T(2.0)*sol(i+1,j,k)*(T(2.0)*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                                    -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                                    -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                 + T(2.0)*sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                    +T(2.0)*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                              -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                 + T(2.0)*sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                    +T(2.0)*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                              -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                 + T(2.0)*sol(i,j,k-1)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                              -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                    +T(2.0)*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                 + T(2.0)*sol(i,j,k+1)*(   -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                              -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                    +T(2.0)*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_aa (int i, int j, int k, Array4<T> const& sol,
                         Array4<std::type_identity_t<T> const> const& rhs, Array4<std::type_identity_t<T> const> const& sig,
                         Array4<int const> const& msk,
                         GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
        T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
        T fxyz = facx + facy + facz;
        T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
        T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
        T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
        T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
        T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
        T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

        T s0 = T(-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                  +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
        T Ax = sol(i,j,k)*s0
            + fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

        sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mlndlap_gscolor_c (int i, int j, int k, Array4<T> const& sol,
                        Array4<std::type_identity_t<T> const> const& rhs, std::type_identity_t<T> sig,
                        Array4<int const> const& msk,
                        GpuArray<std::type_identity_t<T>,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else {
        T facx = T(1.0/36.0)*dxinv[0]*dxinv[0];
        T facy = T(1.0/36.0)*dxinv[1]*dxinv[1];
        T facz = T(1.0/36.0)*dxinv[2]*dxinv[2];
        T fxyz = facx + facy + facz;
        T fmx2y2z = -facx + T(2.0)*facy + T(2.0)*facz;
        T f2xmy2z = T(2.0)*facx - facy + T(2.0)*facz;
        T f2x2ymz = T(2.0)*facx + T(2.0)*facy - facz;
        T f4xm2ym2z = T(4.0)*facx - T(2.0)*facy - T(2.0)*facz;
        T fm2x4ym2z = -T(2.0)*facx + T(4.0)*facy - T(2.0)*facz;
        T fm2xm2y4z = -T(2.0)*facx - T(2.0)*facy + T(4.0)*facz;

        T s0 = T(-4.0)*fxyz*T(8.);
        T Ax = sol(i,j,k)*s0
            + fxyz*(sol(i-1,j-1,k-1)
                  + sol(i+1,j-1,k-1)
                  + sol(i-1,j+1,k-1)
                  + sol(i+1,j+1,k-1)
                  + sol(i-1,j-1,k+1)
                  + sol(i+1,j-1,k+1)
                  + sol(i-1,j+1,k+1)
                  + sol(i+1,j+1,k+1))
            + fmx2y2z*(sol(i  ,j-1,k-1)*T(2.)
                     + sol(i  ,j+1,k-1)*T(2.)
                     + sol(i  ,j-1,k+1)*T(2.)
                     + sol(i  ,j+1,k+1)*T(2.))
            + f2xmy2z*(sol(i-1,j  ,k-1)*T(2.)
                     + sol(i+1,j  ,k-1)*T(2.)
                     + sol(i-1,j  ,k+1)*T(2.)
                     + sol(i+1,j  ,k+1)*T(2.))
            + f2x2ymz*(sol(i-1,j-1,k  )*T(2.)
                     + sol(i+1,j-1,k  )*T(2.)
                     + sol(i-1,j+1,k  )*T(2.)
                     + sol(i+1,j+1,k  )*T(2.))
            + f4xm2ym2z*(sol(i-1,j,k)*T(4.)
                       + sol(i+1,j,k)*T(4.))
            + fm2x4ym2z*(sol(i,j-1,k)*T(4.)
                       + sol(i,j+1,k)*T(4.))
            + fm2xm2y4z*(sol(i,j,k-1)*T(4.)
                       + sol(i,j,k+1)*T(4.));

        sol(i,j,k) += (rhs(i,j,k) - Ax*sig) / (s0*sig);
    }
}

}
#endif
