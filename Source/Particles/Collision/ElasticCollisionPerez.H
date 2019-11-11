#ifndef WARPX_PARTICLES_COLLISION_ELASTIC_COLLISION_PEREZ_H_
#define WARPX_PARTICLES_COLLISION_ELASTIC_COLLISION_PEREZ_H_

#include "WarpXParticleContainer.H"
#include <WarpX.H>
#include <AMReX_REAL.H>
#include <AMReX_DenseBins.H>
#include "ShuffleFisherYates.H"

/* \brief Update particle velocities according to
 *        F. Perez et al., Phys. Plasmas 19, 083104 (2012).
 *        LmdD is max(Debye length, minimal interparticle distance).
 *        L is the Coulomb log. A fixed L will be used if L > 0,
 *        otherwise L will be calculated based on the algorithm.*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumPerezElastic(
    amrex::Real u1x, amrex::Real u1y, amrex::Real u1z,
    amrex::Real u2x, amrex::Real u2y, amrex::Real u2z,
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
    amrex::Real g1 = std::sqrt(1.+(u1x*u1x+u1y*u1y+u1z*u1z)*inv_c2);
    amrex::Real g2 = std::sqrt(1.+(u2x*u2x+u2y*u2y+u2z*u2z)*inv_c2);

    // Compute momenta
    amrex::Real p1x = u1x * m1;
    amrex::Real p1y = u1y * m1;
    amrex::Real p1z = u1z * m1;
    amrex::Real p2x = u2x * m2;
    amrex::Real p2y = u2y * m2;
    amrex::Real p2z = u2z * m2;

    // Compute center-of-mass (COM) velocity and gamma
    buf1 = m1 * g1;
    buf2 = m2 * g2;
    const amrex::Real vcx = (p1x+p2x) / (buf1+buf2);
    const amrex::Real vcy = (p1y+p2y) / (buf1+buf2);
    const amrex::Real vcz = (p1z+p2z) / (buf1+buf2);
    const amrex::Real vcm = std::sqrt(vcx*vcx+vcy*vcy+vcz*vcz);
    const amrex::Real gc  = 1./std::sqrt(1.-vcm*vcm*inv_c2);

    // Compute vc dot v1 and v2
    const amrex::Real vcDv1 = (vcx*u1x + vcy*u1y + vcz*u1z)/g1;
    const amrex::Real vcDv2 = (vcx*u2x + vcy*u2y + vcz*u2z)/g2;

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
    if ( s <= 0.1 )
    { cosXs = 1. + s * std::log(buf1); }
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
    p1x  = ( p1sx*p1sx/buf2) * sinXs*std::cos(buf1) +
           (-p1sy*p1sm/buf2) * sinXs*std::sin(buf1) +
           ( p1sx          ) * cosXs;
    p1y  = ( p1sy*p1sz/buf2) * sinXs*std::cos(buf1) +
           ( p1sx*p1sm/buf2) * sinXs*std::sin(buf1) +
           ( p1sy          ) * cosXs;
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
    }
    buf2 = amrex::Random();
    if ( w1/amrex::max(w1,w2) > buf2 )
    {
      u2x  = p2x / m2;
      u2y  = p2y / m2;
      u2z  = p2z / m2;
    }
}

/* \brief Prepare information for and call
 *        UpdateMomentumPerezElastic().
 *        i1s,i2s is the start index for I1,I2 (inclusive).
 *        i1e,i2e is the start index for I1,I2 (exclusive).
 *        I1 and I2 are the index arrays.
 *        u1 and u2 are the velocity arrays (u=v*gamma),
 *        they could be either different or the same,
 *        their lengths are not needed,
 *        I1 and I2 determine all elements that will be used.
 *        w1 and w2 are arrays of weights.
 *        q1 and q2 are charges. m1 and m2 are masses.
 *        T1 and T2 are temperatures (Joule) and will be used if greater than zero,
 *        otherwise will be computed.
 *        dt is the time step length between two collision calls.
 *        L is the Coulomb log and will be used if greater than zero,
 *        otherwise will be computed.
 *        V is the volume of the corresponding cell.*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ElasticCollisionPerez(
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I1s,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I1e,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I2s,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I2e,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *I1,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *I2,
    amrex::Real* u1x, amrex::Real* u1y, amrex::Real* u1z,
    amrex::Real* u2x, amrex::Real* u2y, amrex::Real* u2z,
    const amrex::Real* w1, const amrex::Real* w2,
    const amrex::Real  q1, const amrex::Real  q2,
    const amrex::Real  m1, const amrex::Real  m2,
    const amrex::Real  T1, const amrex::Real  T2,
    const amrex::Real  dt, const amrex::Real   L, const amrex::Real V)
{

    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);
    int NI1 = I1e - I1s;
    int NI2 = I2e - I2s;

    // shuffle I1 and I2
    ShuffleFisherYates(I1, I1s, I1e);
    ShuffleFisherYates(I2, I2s, I2e);

    // get local T1t and T2t
    amrex::Real T1t; amrex::Real T2t;
    if ( T1 <= 0. )
    {
        amrex::Real vx = 0.; amrex::Real vy = 0.;
        amrex::Real vz = 0.; amrex::Real vs = 0.;
	amrex::Real gm = 0.; amrex::Real us = 0.;
        for (int i1 = I1s; i1 < I1e; ++i1)
        {
          us = ( u1x[ I1[i1] ] * u1x[ I1[i1] ] +
                 u1y[ I1[i1] ] * u1y[ I1[i1] ] +
                 u1z[ I1[i1] ] * u1z[ I1[i1] ] );
          gm = std::sqrt(1.+us*inv_c2);
          vx += u1x[ I1[i1] ] / gm;
	  vy += u1y[ I1[i1] ] / gm;
	  vz += u1z[ I1[i1] ] / gm;
          vs += us / gm / gm;
        }
        vx = vx / NI1; vy = vy / NI1;
	vz = vz / NI1; vs = vs / NI1;
        T1t = m1/3.*(vs-(vx*vx+vy*vy+vz*vz));
    }
    else { T1t = T1; }
    if ( T2 <= 0. )
    {
        amrex::Real vx = 0.; amrex::Real vy = 0.;
        amrex::Real vz = 0.; amrex::Real vs = 0.;
	amrex::Real gm = 0.; amrex::Real us = 0.;
        for (int i2 = I2s; i2 < I2e; ++i2)
        {
          us = ( u2x[ I2[i2] ] * u2x[ I2[i2] ] +
                 u2y[ I2[i2] ] * u2y[ I2[i2] ] +
                 u2z[ I2[i2] ] * u2z[ I2[i2] ] );
          gm = std::sqrt(1.+us*inv_c2);
          vx += u2x[ I2[i2] ] / gm;
	  vy += u2y[ I2[i2] ] / gm;
	  vz += u2z[ I2[i2] ] / gm;
          vs += us / gm / gm;
        }
        vx = vx / NI2; vy = vy / NI2;
	vz = vz / NI2; vs = vs / NI2;
        T2t = m2/3.*(vs-(vx*vx+vy*vy+vz*vz));
    }
    else { T2t = T2; }

    // local density
    amrex::Real n1  = 0.;
    amrex::Real n2  = 0.;
    amrex::Real n12 = 0.;
    for (int i1=I1s; i1<I1e; ++i1) { n1 += w1[ I1[i1] ]; }
    for (int i2=I2s; i2<I2e; ++i2) { n2 += w2[ I2[i2] ]; }
    n1 = n1 / V; n2 = n2 / V;
    {
      int i1 = I1s; int i2 = I2s;
      for (int k = 0; k < amrex::max(NI1,NI2); ++k)
      {
        n12 += amrex::min( w1[ I1[i1] ], w2[ I2[i2] ] );
        ++i1; if ( i1 == I1e ) { i1 = I1s; }
        ++i2; if ( i2 == I2e ) { i2 = I2s; }
      }
      n12 = n12 / V;
    }

    // compute Debye length lmdD
    amrex::Real lmdD;
    lmdD = 1./std::sqrt( n1*q1*q1/(T1t*PhysConst::ep0) +
                         n2*q2*q2/(T2t*PhysConst::ep0) );
    amrex::Real rmin =
        std::pow(4.*MathConst::pi/3.*amrex::max(n1,n2),-1./3.);
    lmdD = amrex::max(lmdD, rmin);

    // call UpdateMomentumPerezElastic()
    {
      int i1 = I1s; int i2 = I2s;
      for (int k = 0; k < amrex::max(NI1,NI2); ++k)
      {
          UpdateMomentumPerezElastic(
              u1x[ I1[i1] ], u1y[ I1[i1] ], u1z[ I1[i1] ],
              u2x[ I2[i2] ], u2y[ I2[i2] ], u2z[ I2[i2] ],
              n1, n2, n12,
              q1, m1, w1[ I1[i1] ], q2, m2, w2[ I2[i2] ],
              dt, L, lmdD);
          ++i1; if ( i1 == I1e ) { i1 = I1s; }
          ++i2; if ( i2 == I2e ) { i2 = I2s; }
      }
    }

}

#endif // WARPX_PARTICLES_COLLISION_ELASTIC_COLLISION_PEREZ_H_
