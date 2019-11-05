#ifndef WARPX_PARTICLES_COLLISION_UPDATEMOMENTUM_PEREZELASTIC_H_
#define WARPX_PARTICLES_COLLISION_UPDATEMOMENTUM_PEREZELASTIC_H_

#include <AMReX_REAL.H>

/* \brief Update particle velocities according to
 *        F. Perez et al., Phys. Plasmas 19, 083104 (2012).
 *        LmdD is max(Debye length, minimal interparticle distance).
 *        L is the Coulomb log. A fixed L will be used if L > 0,
 *        otherwise L will be calculated based on the algorithm.*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumPerezElastic(
    amrex::ParticleReal& v1x, amrex::ParticleReal& v1y, amrex::ParticleReal& v1z,
    amrex::ParticleReal& v2x, amrex::ParticleReal& v2y, amrex::ParticleReal& v2z,
    const amrex::Real n1, const amrex::Real n2, const amrex::Real n12,
    const amrex::Real q1, const amrex::Real m1, const amrex::Real w1,
    const amrex::Real q2, const amrex::Real m2, const amrex::Real w2,
    const amrex::Real dt, const amrex::Real L,  const amrex::Real lmdD)
{
    // Temporary variabls
    amrex::Real buf1;
    amrex::Real buf2;

    constexpr amrex::Real inv_c2 = 1./(PhysConst::c * PhysConst::c);

    // Compute Lorentz factor gamma
    amrex::Real g1 = 1./std::sqrt(1.+(v1x*v1x+v1y*v1y+v1z*v1z)*inv_c2);
    amrex::Real g2 = 1./std::sqrt(1.+(v2x*v2x+v2y*v2y+v2z*v2z)*inv_c2);

    // Compute momenta
    buf1 = m1 * g1;
    amrex::Real p1x = v1x * buf1;
    amrex::Real p1y = v1y * buf1;
    amrex::Real p1z = v1z * buf1;
    buf2 = m2 * g2;
    amrex::Real p2x = v2x * buf2;
    amrex::Real p2y = v2y * buf2;
    amrex::Real p2z = v2z * buf2;

    // Compute center-of-mass (COM) velocity and gamma
    const amrex::Real vcx = (p1x+p2x)/(buf1+buf2);
    const amrex::Real vcy = (p1y+p2y)/(buf1+buf2);
    const amrex::Real vcz = (p1z+p2z)/(buf1+buf2);
    const amrex::Real vcm = std::sqrt(vcx*vcx+vcy*vcy+vcz*vcz);
    const amrex::Real gc  = 1./std::sqrt(1.+vcm*vcm*inv_c2);

    // Compute vc dot v1 and v2
    const amrex::Real vcDv1 = vcx*v1x + vcy*v1y + vcz*v1z;
    const amrex::Real vcDv2 = vcx*v2x + vcy*v2y + vcz*v2z;

    // Compute p1 star
    buf2 = ( (gc-1.)/vcm/vcm*vcDv1 - gc )*buf1;
    const amrex::Real p1sx = p1x + buf2*vcx;
    const amrex::Real p1sy = p1y + buf2*vcy;
    const amrex::Real p1sz = p1z + buf2*vcz;
    const amrex::Real p1sm = std::sqrt(p1sx*p1sx+p1sy*p1sy+p1sz*p1sz);

    // Compute gamma star
    const amrex::Real g1s = (1.-vcDv1*inv_c2)*gc*g1;
    const amrex::Real g2s = (1.-vcDv2*inv_c2)*gc*g2;

    // Compute the Coulomb log lnLmd
    amrex::Real lnLmd;
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
    amrex::Real s = n1*n2/n12 * dt*lnLmd*q1*q1*q2*q2 /
          ( 4. * MathConst::pi * PhysConst::ep0 * PhysConst::ep0 *
            m1*g1*m2*g2/inv_c2/inv_c2 ) * gc*p1sm/(m1*g1+m2*g2) *
            std::pow(m1*g1s*m2*g2s/inv_c2/p1sm/p1sm + 1., 2);

    // Compute s'
    buf1 = (m1*g1+m2*g2)*p1sm/(m1*g1s*m2*g2s*gc); // vrel
    const amrex::Real sp = std::pow(4.*MathConst::pi/3.,1./3.) *
          n1*n2/n12 * dt * buf1 * (m1+m2) /
          amrex::max(m1*std::pow(n1,2./3.),m2*std::pow(n2,2./3.));

    // Determine s
    s = amrex::min(s,sp);

    // Get random numbers
    buf1 = amrex::Random();

    // Compute scattering angle
    amrex::Real cosXs;
    amrex::Real sinXs;
    if ( s <= 0.01 )
    { cosXs = 1. + s * std::log(buf1); }
    else if ( s > 0.01 && s <= 3. )
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
    p1x = ( p1sx*p1sx/buf2) * sinXs*std::cos(buf1) +
          (-p1sy*p1sm/buf2) * sinXs*std::sin(buf1) +
          ( p1sx          ) * cosXs;
    p1y = ( p1sy*p1sz/buf2) * sinXs*std::cos(buf1) +
          ( p1sx*p1sm/buf2) * sinXs*std::sin(buf1) +
          ( p1sy          ) * cosXs;
    p1z = (-buf2          ) * sinXs*std::cos(buf1) +
          ( 0.            ) * sinXs*std::sin(buf1) +
          ( p1sz          ) * cosXs;
    p2x = -p1x;
    p2y = -p1y;
    p2z = -p1z;

    // Transform from COM to lab frame
    // buf1 is vc dot p1 (p1fs)
    // p1 is then updated to be p1f
    buf1 = vcx*p1x + vcy*p1y + vcz*p1z;
    buf2 = (gc-1.)/vcm/vcm*buf1 + m1*g1s*gc;
    p1x = p1x + vcx * buf2;
    p1y = p1y + vcy * buf2;
    p1z = p1z + vcz * buf2;
    // buf1 is vc dot p2 (p2fs)
    // p2 is then updated to be p2f
    buf1 = vcx*p2x + vcy*p2y + vcz*p2z;
    buf2 = (gc-1.)/vcm/vcm*buf1 + m2*g2s*gc;
    p2x = p2x + vcx * buf2;
    p2y = p2y + vcy * buf2;
    p2z = p2z + vcz * buf2;

    // Compute new gamma
    buf1 = std::sqrt(p1x*p1x+p1y*p1y+p1z*p1z);
    buf2 = std::sqrt(p2x*p2x+p2y*p2y+p2z*p2z);
    g1 = std::sqrt(1.+std::pow(buf1/m1/PhysConst::c,2));
    g2 = std::sqrt(1.+std::pow(buf2/m2/PhysConst::c,2));

    // Rejection method
    buf1 = amrex::Random();
    if ( w2/amrex::max(w1,w2) > buf1 )
    {
      buf1 = 1./m1/g1;
      v1x = p1x * buf1;
      v1y = p1y * buf1;
      v1z = p1z * buf1;
    }
    buf2 = amrex::Random();
    if ( w1/amrex::max(w1,w2) > buf2 )
    {
      buf2 = 1./m2/g2;
      v2x = p2x * buf2;
      v2y = p2y * buf2;
      v2z = p2z * buf2;
    }
}

#endif // WARPX_PARTICLES_COLLISION_UPDATEMOMENTUM_PEREZELASTIC_H_
