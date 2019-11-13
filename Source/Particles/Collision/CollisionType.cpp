#include "CollisionType.H"
#include "UpdateMomentumPerezElastic.H"

#include "WarpXParticleContainer.H"
#include <WarpX.H>
#include <AMReX_REAL.H>
#include <AMReX_DenseBins.H>

#include <AMReX_Particles.H>

void CollisionType::CollisionType(
    const std::vector<std::string>& species_names,
    const std::vector<std::string>  collision_name)
{

    std::vector<std::string> collision_species;

    ParmParse pp(collision_name);
    pp.getarr("species", collision_species);

    for (int i=0; i<species_names.size(); i++)
    {
        if (species_names[i] == collision_species[0])
	{ m_species1 = i; }
        if (species_names[i] == collision_species[1])
	{ m_species2 = i; }
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
void CollisionType::ElasticCollisionPerez(
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

/* \brief Shuffle array according to Fisher-Yates algorithm.
 *        Only shuffle the part between is <= i < ie, n = ie-is.*/
void CollisionType::ShuffleFisherYates(
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *array,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const is,
    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const ie)
{
    int j; int buf;
    for (int i = ie-1; i >= is+1; --i)
    {
        // get random number is <= j <= i
        while ( true ) {
            j = (int) floor(amrex::Random()*(i-is+1)) + is;
            if ( j >= is && j <= i ) { break; }
        }
        buf = array[i];
        array[i] = array[j];
        array[j] = buf;
    }
}
