#ifndef AMREX_EB_STATE_REDIST_SLOPE_LIMITER_K_H_
#define AMREX_EB_STATE_REDIST_SLOPE_LIMITER_K_H_

namespace amrex {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
amrex_calc_alpha_stencil(Real q_hat, Real q_max, Real q_min, Real state) noexcept
{
#ifdef AMREX_USE_FLOAT
    constexpr Real epsilon = amrex::Real(1.e-6);
#else
    constexpr Real epsilon = 1.e-12;
#endif

    const Real sml = epsilon*amrex::max(amrex::Math::abs(q_max),amrex::Math::abs(q_min));
    Real alpha;

    if ((q_hat-state) > sml) {
        alpha = amrex::min(1.0_rt,(q_max-state)/(q_hat-state));
    } else if ((q_hat-state) < -sml) {
        alpha = amrex::min(1.0_rt,(q_min-state)/(q_hat-state));
    } else {
        alpha = 1.0_rt;
    }
    return alpha;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
amrex_apply_centroid_limiter(int i, int j, int k, int n,
                             amrex::Array4<amrex::Real const> const& state,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             AMREX_D_DECL(amrex::Real& slope_x,
                                          amrex::Real& slope_y,
                                          amrex::Real& slope_z),
                             amrex::Array4<amrex::Real const> const& ccent) noexcept
{
#ifdef AMREX_USE_FLOAT
    constexpr Real epsilon = amrex::Real(1.e-6);
#else
    constexpr Real epsilon = 1.e-12;
#endif

    AMREX_D_TERM(amrex::Real xalpha = 1.0;,
                 amrex::Real yalpha = 1.0;,
                 amrex::Real zalpha = 1.0;);

    // Load center cell flag once to avoid repeated array accesses
    auto const cell_flag = flag(i,j,k);

    // Compute the limiters needed to keep the predicted q_hat between the max and min
#if (AMREX_SPACEDIM == 2)
    int kk = 0;
#elif (AMREX_SPACEDIM == 3)
#if defined(AMREX_USE_CUDA)
#pragma unroll 1
#endif
    for(int kk(-1); kk<=1; kk++)
#endif
    {
#if defined(AMREX_USE_CUDA)
#pragma unroll 1
#endif
     for(int jj(-1); jj<=1; jj++){
#if defined(AMREX_USE_CUDA)
#pragma unroll 1
#endif
      for(int ii(-1); ii<=1; ii++){

        if (cell_flag.isConnected(ii,jj,kk) && AMREX_D_TERM(xalpha > 0.0_rt, || yalpha > 0.0_rt, || zalpha > 0.0_rt))
        {
            AMREX_D_TERM(Real delta_x = ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0) + static_cast<Real>(ii);,
                         Real delta_y = ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1) + static_cast<Real>(jj);,
                         Real delta_z = ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2) + static_cast<Real>(kk););

            // Capture which dimensions are affected before deltas go out of scope
            AMREX_D_TERM(bool const affects_x = amrex::Math::abs(delta_x) > epsilon;,
                         bool const affects_y = amrex::Math::abs(delta_y) > epsilon;,
                         bool const affects_z = amrex::Math::abs(delta_z) > epsilon;);

            Real q_hat = state(i,j,k,n) + AMREX_D_TERM(  delta_x * slope_x,
                                                       + delta_y * slope_y,
                                                       + delta_z * slope_z);

            Real q_max = amrex::max(state(i+ii,j+jj,k+kk,n),state(i,j,k,n));
            Real q_min = amrex::min(state(i+ii,j+jj,k+kk,n),state(i,j,k,n));

            if ( q_hat-q_max > amrex::Math::abs(epsilon*q_max) || q_hat-q_min < -1.0*amrex::Math::abs(epsilon*q_min) )
            {
                Real new_lim = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n));

                AMREX_D_TERM(if (affects_x) { xalpha = amrex::min(xalpha,new_lim); },
                             if (affects_y) { yalpha = amrex::min(yalpha,new_lim); },
                             if (affects_z) { zalpha = amrex::min(zalpha,new_lim); });
            }
        }
      }
     }
    }

    AMREX_D_TERM(slope_x *= xalpha;,
                 slope_y *= yalpha;,
                 slope_z *= zalpha;);
}

}

#endif
