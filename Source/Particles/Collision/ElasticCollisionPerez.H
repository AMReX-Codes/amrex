#ifndef WARPX_PARTICLES_COLLISION_ELASTIC_COLLISION_PEREZ_H_
#define WARPX_PARTICLES_COLLISION_ELASTIC_COLLISION_PEREZ_H_

#include "UpdateMomentumPerezElastic.H"
#include "ComputeTemperature.H"
#include <WarpXConst.H>
#include <AMReX_Random.H>

/* \brief Prepare information for and call
 *        UpdateMomentumPerezElastic().
 *        I1s,I2s is the start index for I1,I2 (inclusive).
 *        I1e,I2e is the start index for I1,I2 (exclusive).
 *        I1 and I2 are the index arrays.
 *        u1 and u2 are the velocity arrays (u=v*gamma),
 *        they could be either different or the same,
 *        their lengths are not needed,
 *        I1 and I2 determine all elements that will be used.
 *        w1 and w2 are arrays of weights.
 *        q1 and q2 are charges. m1 and m2 are masses.
 *        T1 and T2 are temperatures (Joule)
 *        and will be used if greater than zero,
 *        otherwise will be computed.
 *        dt is the time step length between two collision calls.
 *        L is the Coulomb log and will be used if greater than zero,
 *        otherwise will be computed.
 *        dV is the volume of the corresponding cell.
*/

template <typename T_index, typename T_R>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ElasticCollisionPerez (
    T_index const I1s, T_index const I1e,
    T_index const I2s, T_index const I2e,
    T_index *I1,       T_index *I2,
    T_R *u1x, T_R *u1y, T_R *u1z,
    T_R *u2x, T_R *u2y, T_R *u2z,
    T_R const *w1, T_R const *w2,
    T_R const  q1, T_R const  q2,
    T_R const  m1, T_R const  m2,
    T_R const  T1, T_R const  T2,
    T_R const  dt, T_R const   L, T_R const dV)
{

    T_R constexpr inv_c2 = 1./(PhysConst::c*PhysConst::c);
    int NI1 = I1e - I1s;
    int NI2 = I2e - I2s;

    // get local T1t and T2t
    T_R T1t; T_R T2t;
    if ( T1 <= 0. && L <= 0. )
    {
        T1t = ComputeTemperature(I1s,I1e,I1,u1x,u1y,u1z,m1);
    }
    else { T1t = T1; }
    if ( T2 <= 0. && L <= 0. )
    {
        T2t = ComputeTemperature(I2s,I2e,I2,u2x,u2y,u2z,m2);
    }
    else { T2t = T2; }

    // local density
    T_R n1  = 0.;
    T_R n2  = 0.;
    T_R n12 = 0.;
    for (int i1=I1s; i1<I1e; ++i1) { n1 += w1[ I1[i1] ]; }
    for (int i2=I2s; i2<I2e; ++i2) { n2 += w2[ I2[i2] ]; }
    n1 = n1 / dV; n2 = n2 / dV;
    {
      int i1 = I1s; int i2 = I2s;
      for (int k = 0; k < amrex::max(NI1,NI2); ++k)
      {
        n12 += amrex::min( w1[ I1[i1] ], w2[ I2[i2] ] );
        ++i1; if ( i1 == I1e ) { i1 = I1s; }
        ++i2; if ( i2 == I2e ) { i2 = I2s; }
      }
      n12 = n12 / dV;
    }

    // compute Debye length lmdD
    T_R lmdD;
    lmdD = T_R(1.)/std::sqrt( n1*q1*q1/(T1t*PhysConst::ep0) +
                         n2*q2*q2/(T2t*PhysConst::ep0) );
    T_R rmin =
        std::pow(T_R(4.)*MathConst::pi/T_R(3.)*amrex::max(n1,n2),T_R(-1./3.));
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
