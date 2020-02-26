/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
#define WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_

#include "Utils/WarpXConst.H"
#include <AMReX_Random.H>

#include <cmath>  // isnan() isinf()
#include <limits> // numeric_limits<float>::min()


/* \brief Update particle velocities according to
 *        F. Perez et al., Phys.Plasmas.19.083104 (2012),
 *        which is based on Nanbu's method, PhysRevE.55.4642 (1997).
 *        @param[in] LmdD is max(Debye length, minimal interparticle distance).
 *        @param[in] L is the Coulomb log. A fixed L will be used if L > 0,
 *        otherwise L will be calculated based on the algorithm.
 *        To see if there are nan or inf updated velocities,
 *        compile with USE_ASSERTION=TRUE.
*/

template <typename T_R>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumPerezElastic (
    T_R& u1x, T_R& u1y, T_R& u1z, T_R& u2x, T_R& u2y, T_R& u2z,
    T_R const n1, T_R const n2, T_R const n12,
    T_R const q1, T_R const m1, T_R const w1,
    T_R const q2, T_R const m2, T_R const w2,
    T_R const dt, T_R const L,  T_R const lmdD)
{

    // If g = u1 - u2 = 0, do not collide.
    if ( std::abs(u1x-u2x) < std::numeric_limits<T_R>::min() &&
         std::abs(u1y-u2y) < std::numeric_limits<T_R>::min() &&
         std::abs(u1z-u2z) < std::numeric_limits<T_R>::min() )
    { return; }

    T_R constexpr inv_c2 = T_R(1.0) / ( PhysConst::c * PhysConst::c );

    // Compute Lorentz factor gamma
    T_R const g1 = std::sqrt( T_R(1.0) + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_c2 );
    T_R const g2 = std::sqrt( T_R(1.0) + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_c2 );

    // Compute momenta
    T_R const p1x = u1x * m1;
    T_R const p1y = u1y * m1;
    T_R const p1z = u1z * m1;
    T_R const p2x = u2x * m2;
    T_R const p2y = u2y * m2;
    T_R const p2z = u2z * m2;

    // Compute center-of-mass (COM) velocity and gamma
    T_R const mass_g = m1 * g1 + m2 * g2;
    T_R const vcx    = (p1x+p2x) / mass_g;
    T_R const vcy    = (p1y+p2y) / mass_g;
    T_R const vcz    = (p1z+p2z) / mass_g;
    T_R const vcms   = vcx*vcx + vcy*vcy + vcz*vcz;
    T_R const gc     = T_R(1.0) / std::sqrt( T_R(1.0) - vcms*inv_c2 );

    // Compute vc dot v1 and v2
    T_R const vcDv1 = (vcx*u1x + vcy*u1y + vcz*u1z) / g1;
    T_R const vcDv2 = (vcx*u2x + vcy*u2y + vcz*u2z) / g2;

    // Compute p1 star
    T_R p1sx;
    T_R p1sy;
    T_R p1sz;
    if ( vcms > std::numeric_limits<T_R>::min() )
    {
        T_R const lorentz_tansform_factor =
            ( (gc-T_R(1.0))/vcms*vcDv1 - gc )*m1*g1;
        p1sx = p1x + vcx*lorentz_tansform_factor;
        p1sy = p1y + vcy*lorentz_tansform_factor;
        p1sz = p1z + vcz*lorentz_tansform_factor;
    }
    else // If vcms = 0, don't do Lorentz-transform.
    {
        p1sx = p1x;
        p1sy = p1y;
        p1sz = p1z;
    }
    T_R const p1sm = std::sqrt( p1sx*p1sx + p1sy*p1sy + p1sz*p1sz );

    // Compute gamma star
    T_R const g1s = ( T_R(1.0) - vcDv1*inv_c2 )*gc*g1;
    T_R const g2s = ( T_R(1.0) - vcDv2*inv_c2 )*gc*g2;

    // Compute the Coulomb log lnLmd
    T_R lnLmd;
    if ( L > T_R(0.0) ) { lnLmd = L; }
    else
    {
        // Compute b0
        T_R const b0 = std::abs(q1*q2) * inv_c2 /
               (T_R(4.0)*MathConst::pi*PhysConst::ep0) * gc/mass_g *
               ( m1*g1s*m2*g2s/(p1sm*p1sm*inv_c2) + T_R(1.0) );

        // Compute the minimal impact parameter
        T_R bmin = amrex::max(PhysConst::hbar*MathConst::pi/p1sm,b0);

        // Compute the Coulomb log lnLmd
        lnLmd = amrex::max( T_R(2.0),
                T_R(0.5)*std::log(T_R(1.0)+lmdD*lmdD/(bmin*bmin)) );
    }

    // Compute s
    T_R s = n1*n2/n12 * dt*lnLmd*q1*q1*q2*q2 /
          ( T_R(4.0) * MathConst::pi * PhysConst::ep0 * PhysConst::ep0 *
            m1*g1*m2*g2/(inv_c2*inv_c2) ) * gc*p1sm/mass_g *
            std::pow(m1*g1s*m2*g2s/(inv_c2*p1sm*p1sm) + T_R(1.0), 2.0);

    // Compute s'
    T_R const vrel = mass_g*p1sm/(m1*g1s*m2*g2s*gc);
    T_R const sp = std::pow(T_R(4.0)*MathConst::pi/T_R(3.0),T_R(1.0/3.0)) *
          n1*n2/n12 * dt * vrel * (m1+m2) /
          amrex::max( m1*std::pow(n1,T_R(2.0/3.0)),
                      m2*std::pow(n2,T_R(2.0/3.0)) );

    // Determine s
    s = amrex::min(s,sp);

    // Get random numbers
    T_R r = amrex::Random();

    // Compute scattering angle
    T_R cosXs;
    T_R sinXs;
    if ( s <= T_R(0.1) )
    {
        while ( true )
        {
            cosXs = T_R(1.0) + s * std::log(r);
            // Avoid the bug when r is too small such that cosXs < -1
            if ( cosXs >= T_R(-1.0) ) { break; }
            r = amrex::Random();
        }
    }
    else if ( s > T_R(0.1) && s <= T_R(3.0) )
    {
        T_R const Ainv = 0.0056958 + 0.9560202*s - 0.508139*s*s +
            0.47913906*s*s*s - 0.12788975*s*s*s*s + 0.02389567*s*s*s*s*s;
        cosXs = Ainv * std::log( std::exp(T_R(-1.0)/Ainv) +
                T_R(2.0) * r * std::sinh(T_R(1.0)/Ainv) );
    }
    else if ( s > T_R(3.0) && s <= T_R(6.0) )
    {
        T_R const A = T_R(3.0) * std::exp(-s);
        cosXs = T_R(1.0)/A * std::log( std::exp(-A) +
                T_R(2.0) * r * std::sinh(A) );
    }
    else
    {
        cosXs = T_R(2.0) * r - T_R(1.0);
    }
    sinXs = std::sqrt(T_R(1.0) - cosXs*cosXs);

    // Get random azimuthal angle
    T_R const phis = amrex::Random() * T_R(2.0) * MathConst::pi;
    T_R const cosphis = std::cos(phis);
    T_R const sinphis = std::sin(phis);

    // Compute post-collision momenta pfs in COM
    T_R p1fsx;
    T_R p1fsy;
    T_R p1fsz;
    // p1sp is the p1s perpendicular
    T_R p1sp = std::sqrt( p1sx*p1sx + p1sy*p1sy );
    // Make sure p1sp is not almost zero
    if ( p1sp > std::numeric_limits<T_R>::min() )
    {
        p1fsx = ( p1sx*p1sz/p1sp ) * sinXs*cosphis +
                ( p1sy*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sx           ) * cosXs;
        p1fsy = ( p1sy*p1sz/p1sp ) * sinXs*cosphis +
                (-p1sx*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sy           ) * cosXs;
        p1fsz = (-p1sp           ) * sinXs*cosphis +
                ( T_R(0.0)       ) * sinXs*sinphis +
                ( p1sz           ) * cosXs;
        // Note a negative sign is different from
        // Eq. (12) in Perez's paper,
        // but they are the same due to the random nature of phis.
    }
    else
    {
        // If the previous p1sp is almost zero
        // x->y  y->z  z->x
        // This set is equivalent to the one in Nanbu's paper
        p1sp = std::sqrt( p1sy*p1sy + p1sz*p1sz );
        p1fsy = ( p1sy*p1sx/p1sp ) * sinXs*cosphis +
                ( p1sz*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sy           ) * cosXs;
        p1fsz = ( p1sz*p1sx/p1sp ) * sinXs*cosphis +
                (-p1sy*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sz           ) * cosXs;
        p1fsx = (-p1sp           ) * sinXs*cosphis +
                ( T_R(0.0)       ) * sinXs*sinphis +
                ( p1sx           ) * cosXs;
    }

    T_R const p2fsx = -p1fsx;
    T_R const p2fsy = -p1fsy;
    T_R const p2fsz = -p1fsz;

    // Transform from COM to lab frame
    T_R p1fx;    T_R p2fx;
    T_R p1fy;    T_R p2fy;
    T_R p1fz;    T_R p2fz;
    if ( vcms > std::numeric_limits<T_R>::min() )
    {
        T_R const vcDp1fs = vcx*p1fsx + vcy*p1fsy + vcz*p1fsz;
        T_R const vcDp2fs = vcx*p2fsx + vcy*p2fsy + vcz*p2fsz;
        T_R const factor = (gc-T_R(1.0))/vcms;
        T_R const factor1 = factor*vcDp1fs + m1*g1s*gc;
        T_R const factor2 = factor*vcDp2fs + m2*g2s*gc;
        p1fx = p1fsx + vcx * factor1;
        p1fy = p1fsy + vcy * factor1;
        p1fz = p1fsz + vcz * factor1;
        p2fx = p2fsx + vcx * factor2;
        p2fy = p2fsy + vcy * factor2;
        p2fz = p2fsz + vcz * factor2;
    }
    else // If vcms = 0, don't do Lorentz-transform.
    {
        p1fx = p1fsx;
        p1fy = p1fsy;
        p1fz = p1fsz;
        p2fx = p2fsx;
        p2fy = p2fsy;
        p2fz = p2fsz;
    }

    // Rejection method
    r = amrex::Random();
    if ( w2 > r*amrex::max(w1, w2) )
    {
        u1x  = p1fx / m1;
        u1y  = p1fy / m1;
        u1z  = p1fz / m1;
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
    }
    r = amrex::Random();
    if ( w1 > r*amrex::max(w1, w2) )
    {
        u2x  = p2fx / m2;
        u2y  = p2fy / m2;
        u2z  = p2fz / m2;
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
    }

}

#endif // WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
