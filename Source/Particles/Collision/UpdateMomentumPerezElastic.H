#ifndef WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
#define WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_

#include <WarpXConst.H>
#include <AMReX_Random.H>
#include <cmath> // isnan() isinf()

/* \brief Update particle velocities according to
 *        F. Perez et al., Phys. Plasmas 19, 083104 (2012).
 *        LmdD is max(Debye length, minimal interparticle distance).
 *        L is the Coulomb log. A fixed L will be used if L > 0,
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

    // Temporary variabls
    T_R buf1; T_R buf2;

    T_R constexpr inv_c2 = 1./(PhysConst::c * PhysConst::c);

    // Compute Lorentz factor gamma
    T_R g1 = std::sqrt(1.+(u1x*u1x+u1y*u1y+u1z*u1z)*inv_c2);
    T_R g2 = std::sqrt(1.+(u2x*u2x+u2y*u2y+u2z*u2z)*inv_c2);

    // Compute momenta
    T_R p1x = u1x * m1;
    T_R p1y = u1y * m1;
    T_R p1z = u1z * m1;
    T_R p2x = u2x * m2;
    T_R p2y = u2y * m2;
    T_R p2z = u2z * m2;

    // Compute center-of-mass (COM) velocity and gamma
    buf1 = m1 * g1;
    buf2 = m2 * g2 + buf1;
    T_R const vcx = (p1x+p2x) / buf2;
    T_R const vcy = (p1y+p2y) / buf2;
    T_R const vcz = (p1z+p2z) / buf2;
    T_R const vcm = std::sqrt(vcx*vcx+vcy*vcy+vcz*vcz);
    T_R const gc  = 1./std::sqrt(1.-vcm*vcm*inv_c2);

    // Avoid 1/0 -> nan
    if ( vcm == 0. ) { return; }

    // Compute vc dot v1 and v2
    T_R const vcDv1 = (vcx*u1x + vcy*u1y + vcz*u1z)/g1;
    T_R const vcDv2 = (vcx*u2x + vcy*u2y + vcz*u2z)/g2;

    // Compute p1 star
    buf2 = ( (gc-1.)/vcm/vcm*vcDv1 - gc )*buf1;
    T_R const p1sx = p1x + buf2*vcx;
    T_R const p1sy = p1y + buf2*vcy;
    T_R const p1sz = p1z + buf2*vcz;
    T_R const p1sm = std::sqrt(p1sx*p1sx+p1sy*p1sy+p1sz*p1sz);

    // Compute gamma star
    T_R const g1s = (1.-vcDv1*inv_c2)*gc*g1;
    T_R const g2s = (1.-vcDv2*inv_c2)*gc*g2;

    // Compute the Coulomb log lnLmd
    T_R lnLmd;
    if ( L > 0. )
    { lnLmd = L; }
    else
    {
        // Compute b0 (buf1)
        buf1 = std::abs(q1*q2) * inv_c2 /
               (4.*MathConst::pi*PhysConst::ep0) * gc/(m1*g1+m2*g2) *
               ( m1*g1s*m2*g2s/p1sm/p1sm/inv_c2 + 1. );

        // Compute the minimal impact parameter (buf2)
        buf2 = amrex::max(PhysConst::hbar*MathConst::pi/p1sm,buf1);

        // Compute the Coulomb log lnLmd
        lnLmd = amrex::max( 2., 0.5*std::log(1.+lmdD*lmdD/buf2/buf2) );
    }

    // Compute s
    T_R s = n1*n2/n12 * dt*lnLmd*q1*q1*q2*q2 /
          ( 4. * MathConst::pi * PhysConst::ep0 * PhysConst::ep0 *
            m1*g1*m2*g2/inv_c2/inv_c2 ) * gc*p1sm/(m1*g1+m2*g2) *
            std::pow(m1*g1s*m2*g2s/inv_c2/p1sm/p1sm + 1., 2);

    // Compute s'
    buf1 = (m1*g1+m2*g2)*p1sm/(m1*g1s*m2*g2s*gc); // vrel
    T_R const sp = std::pow(4.*MathConst::pi/3.,1./3.) *
          n1*n2/n12 * dt * buf1 * (m1+m2) /
          amrex::max(m1*std::pow(n1,2./3.),m2*std::pow(n2,2./3.));

    // Determine s
    s = amrex::min(s,sp);

    // Get random numbers
    buf1 = amrex::Random();

    // Compute scattering angle
    T_R cosXs;
    T_R sinXs;
    if ( s <= 0.1 )
    {
        while ( true )
        {
            cosXs = 1. + s * std::log(buf1);
            // avoid the bug when buf1 is too small
            if ( cosXs >= -1.0 ) { break; }
            buf1 = amrex::Random();
        }
    }
    else if ( s > 0.1 && s <= 3. )
    {
        // buf2 is A inverse
        buf2 = 0.0056958 + 0.9560202*s - 0.508139*s*s +
               0.47913906*s*s*s - 0.12788975*s*s*s*s + 0.02389567*s*s*s*s*s;
        cosXs = buf2*std::log(std::exp(-1./buf2)+2.*buf1*std::sinh(1.0/buf2));
    }
    else if ( s > 3. && s <= 6. )
    {
        // buf2 is A
        buf2 = 3. * std::exp(-s);
        cosXs = 1./buf2*std::log(std::exp(-buf2)+2.*buf1*std::sinh(buf2));
    }
    else
    {
        cosXs = 2.*buf1 - 1.;
    }
    sinXs = std::sqrt(1. - cosXs*cosXs);

    // Get random azimuthal angle
    buf1 = amrex::Random() * 2. * MathConst::pi;

    // Compute post-collision momenta in COM
    // buf2 is the p1s perpendicular
    // p1 is now p1f star
    buf2 = std::sqrt(p1sx*p1sx + p1sy*p1sy);
    // Avoid 1/0 -> nan
    if ( buf2 == 0. ) { p1x = 0.; p1y = 0.; }
    else
    {
        p1x  = ( p1sx*p1sz/buf2) * sinXs*std::cos(buf1) +
               (-p1sy*p1sm/buf2) * sinXs*std::sin(buf1) +
               ( p1sx          ) * cosXs;
        p1y  = ( p1sy*p1sz/buf2) * sinXs*std::cos(buf1) +
               ( p1sx*p1sm/buf2) * sinXs*std::sin(buf1) +
               ( p1sy          ) * cosXs;
    }
    p1z  = (-buf2          ) * sinXs*std::cos(buf1) +
           ( 0.            ) * sinXs*std::sin(buf1) +
           ( p1sz          ) * cosXs;
    p2x  = -p1x;
    p2y  = -p1y;
    p2z  = -p1z;

    // Transform from COM to lab frame
    // buf1 is vc dot p1 (p1fs)
    // p1 is then updated to be p1f
    buf1 = vcx*p1x + vcy*p1y + vcz*p1z;
    buf2 = (gc-1.)/vcm/vcm*buf1 + m1*g1s*gc;
    p1x  = p1x + vcx * buf2;
    p1y  = p1y + vcy * buf2;
    p1z  = p1z + vcz * buf2;
    // buf1 is vc dot p2 (p2fs)
    // p2 is then updated to be p2f
    buf1 = vcx*p2x + vcy*p2y + vcz*p2z;
    buf2 = (gc-1.)/vcm/vcm*buf1 + m2*g2s*gc;
    p2x  = p2x + vcx * buf2;
    p2y  = p2y + vcy * buf2;
    p2z  = p2z + vcz * buf2;

    // Rejection method
    buf1 = amrex::Random();
    if ( w2/amrex::max(w1,w2) > buf1 )
    {
        u1x  = p1x / m1;
        u1y  = p1y / m1;
        u1z  = p1z / m1;
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
    }
    buf2 = amrex::Random();
    if ( w1/amrex::max(w1,w2) > buf2 )
    {
        u2x  = p2x / m2;
        u2y  = p2y / m2;
        u2z  = p2z / m2;
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
    }

}

#endif // WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
