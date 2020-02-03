/* Copyright 2019-2020 Glenn Richardson, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WarpX_QED_K_h
#define WarpX_QED_K_h

#include <AMReX_FArrayBox.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

/**
 * calc_M calculates the Magnetization field of the vacuum at a specific point and returns it as a three component vector
 * \param[in] arr This is teh empty array that will be filled with the components of the M-field
 * \param[in] ex The x-component of the E-field at the point at whicht the M-field is to be calculated
 * \param[in] ey The y-component of the E-field at the point at whicht the M-field is to be calculated
 * \param[in] ez The z-component of the E-field at the point at whicht the M-field is to be calculated
 * \param[in] bx The x-component of the B-field at the point at whicht the M-field is to be calculated
 * \param[in] by The y-component of the B-field at the point at whicht the M-field is to be calculated
 * \param[in] bz The z-component of the B-field at the point at whicht the M-field is to be calculated
 * \param[in] xi The quantum parameter being used for the simulation
 * \param[in] c2 The speed of light squared
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void calc_M(Real arr [], Real ex, Real ey, Real ez, Real bx, Real by, Real bz, Real xi, Real c2)
{
    const Real ee = ex*ex+ey*ey+ez*ez;
    const Real bb = bx*bx+by*by+bz*bz;
    const Real eb = ex*bx+ey*by+ez*bz;
    arr[0] = -2*xi*c2*( 2*bx*(ee-c2*bb) - 7*ex*eb );
    arr[1] = -2*xi*c2*( 2*by*(ee-c2*bb) - 7*ey*eb );
    arr[2] = -2*xi*c2*( 2*bz*(ee-c2*bb) - 7*ez*eb );
};


/**
 * warpx_hybrid_QED_push uses an FDTD scheme to calculate QED corrections to
 * Maxwell's equations and preforms a half timestep correction to the E-fields
 *
 * \param[in] Ex This function modifies the Ex field at the end
 * \param[in] Ey This function modifies the Ey field at the end
 * \param[in] Ez This function modifies the Ez field at the end
 * \param[in] Bx The QED corrections are non-linear functions of B. However,
 *            they do not modify B itself
 * \param[in] By The QED corrections are non-linear functions of B. However,
 *            they do not modify B itself
 * \param[in] Bz The QED corrections are non-linear functions of B. However,
 *            they do not modify B itself
 * \param[in] tempEx Since the corrections to E at a given node are non-linear functions
 *            of the values of E on the surronding nodes, temp arrays are used so that
 *            modifications to one node do not influence the corrections to the surronding nodes
 * \param[in] tempEy Since the corrections to E at a given node are non-linear functions
 *            of the values of E on the surronding nodes, temp arrays are used so that modifications to
 *            one node do not influence the corrections to the surronding nodes
 * \param[in] tempEz Since the corrections to E at a given node are non-linear functions
 *            of the values of E on the surronding nodes, temp arrays are used so that modifications to
 *            one node do not influence the corrections to the surronding nodes
 * \param[in] dx The x spatial step, used for calculating curls
 * \param[in] dy The y spatial step, used for calculating curls
 * \param[in] dz The z spatial step, used for calulating curls
 * \param[in] dt The temporal step, used for the half push/correction to the E-fields at the end of the function
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_hybrid_QED_push (int j, int k, int l, Array4<Real> const& Ex, Array4<Real>
                            const& Ey, Array4<Real> const& Ez, Array4<Real> const& Bx,
                            Array4<Real> const& By, Array4<Real const> const& Bz,
                            Array4<Real> const& tmpEx, Array4<Real> const& tmpEy,
                            Array4<Real> const& tmpEz, Real dx, Real dy, Real dz, Real dt)
{


// Defining constants to be used in the calculations

constexpr amrex::Real c2 = PhysConst::c * PhysConst::c;
const amrex::Real xi = WarpX::quantum_xi;
const amrex::Real dxi = 1./dx;
const amrex::Real dzi = 1./dz;

#if (AMREX_SPACEDIM == 3)
const amrex::Real dyi = 1./dy;

    // Picking out points for stencil to be used in curl function of M

    amrex::Real Mpx [3] = {0.,0.,0.};
    amrex::Real Mnx [3] = {0.,0.,0.};
    amrex::Real Mpy [3] = {0.,0.,0.};
    amrex::Real Mny [3] = {0.,0.,0.};
    amrex::Real Mpz [3] = {0.,0.,0.};
    amrex::Real Mnz [3] = {0.,0.,0.};

    // Calcualting the M-field at the chosen stencil points

    calc_M(Mpx, tmpEx(j+1,k,l), tmpEy(j+1,k,l), tmpEz(j+1,k,l),
                Bx(j+1,k,l), By(j+1,k,l), Bz(j+1,k,l), xi, c2);
    calc_M(Mnx, tmpEx(j-1,k,l), tmpEy(j-1,k,l), tmpEz(j-1,k,l),
                Bx(j-1,k,l), By(j-1,k,l), Bz(j-1,k,l), xi, c2);
    calc_M(Mpy, tmpEx(j,k+1,l), tmpEy(j,k+1,l), tmpEz(j,k+1,l),
                Bx(j,k+1,l), By(j,k+1,l), Bz(j,k+1,l), xi, c2);
    calc_M(Mny, tmpEx(j,k-1,l), tmpEy(j,k-1,l), tmpEz(j,k-1,l),
                Bx(j,k-1,l), By(j,k-1,l), Bz(j,k-1,l), xi, c2);
    calc_M(Mpz, tmpEx(j,k,l+1), tmpEy(j,k,l+1), tmpEz(j,k,l+1),
                Bx(j,k,l+1), By(j,k,l+1), Bz(j,k,l+1), xi, c2);
    calc_M(Mnz, tmpEx(j,k,l-1), tmpEy(j,k,l-1), tmpEz(j,k,l-1),
                Bx(j,k,l-1), By(j,k,l-1), Bz(j,k,l-1), xi, c2);

    // Calculating necessary curls

    const amrex::Real VxM[3] = {
        0.5*( (Mpy[2]-Mny[2])*dyi - (Mpz[1]-Mnz[1])*dzi ),
        0.5*( (Mpz[0]-Mnz[0])*dzi - (Mpx[2]-Mnx[2])*dxi ),
        0.5*( (Mpx[1]-Mnx[1])*dxi - (Mpy[0]-Mny[0])*dyi ),
    };

    const amrex::Real VxE[3] = {
        0.5*( (tmpEz(j,k+1,l)-tmpEz(j,k-1,l) )*dyi - (tmpEy(j,k,l+1)-tmpEy(j,k,l-1) )*dzi ),
        0.5*( (tmpEx(j,k,l+1)-tmpEx(j,k,l-1) )*dzi - (tmpEz(j+1,k,l)-tmpEz(j-1,k,l) )*dxi ),
        0.5*( (tmpEy(j+1,k,l)-tmpEy(j-1,k,l) )*dxi - (tmpEx(j,k+1,l)-tmpEx(j,k-1,l) )*dyi ),
    };

    const amrex::Real VxB[3] = {
        0.5*( (Bz(j,k+1,l)-Bz(j,k-1,l) )*dyi - (By(j,k,l+1)-By(j,k,l-1) )*dzi ),
        0.5*( (Bx(j,k,l+1)-Bx(j,k,l-1) )*dzi - (Bz(j+1,k,l)-Bz(j-1,k,l) )*dxi ),
        0.5*( (By(j+1,k,l)-By(j-1,k,l) )*dxi - (Bx(j,k+1,l)-Bx(j,k-1,l) )*dyi ),
    };

    // Defining comapct values for QED corrections

    const amrex::Real ex = tmpEx(j,k,l);
    const amrex::Real ey = tmpEy(j,k,l);
    const amrex::Real ez = tmpEz(j,k,l);
    const amrex::Real bx = Bx(j,k,l);
    const amrex::Real by = By(j,k,l);
    const amrex::Real bz = Bz(j,k,l);
    const amrex::Real ee = ex*ex + ey*ey + ez*ez;
    const amrex::Real bb = bx*bx + by*by + bz*bz;
    const amrex::Real eb = ex*bx + ey*by + ez*bz;
    const amrex::Real EVxE = ex*VxE[0] + ey*VxE[1] + ez*VxE[2];
    const amrex::Real BVxE = bx*VxE[0] + by*VxE[1] + bz*VxE[2];
    const amrex::Real EVxB = ex*VxB[0] + ey*VxB[1] + ez*VxB[2];
    const amrex::Real BVxB = bx*VxB[0] + by*VxB[1] + bz*VxB[2];

    const amrex::Real beta = 4*xi*( ee - c2*bb ) + PhysConst::ep0;

    const amrex::Real Alpha[3] = {
        2*xi*c2*( -7*bx*EVxE - 7*VxE[0]*eb + 4*ex*BVxE ) + VxM[0],
        2*xi*c2*( -7*by*EVxE - 7*VxE[1]*eb + 4*ey*BVxE ) + VxM[1],
        2*xi*c2*( -7*bz*EVxE - 7*VxE[2]*eb + 4*ez*BVxE ) + VxM[2]
    };

    const amrex::Real Omega[3] = {
        Alpha[0] + 2*xi*c2*( 4*ex*EVxB + 2*VxB[0]*( ee - c2*bb ) + 7*c2*bx*BVxB ),
        Alpha[1] + 2*xi*c2*( 4*ey*EVxB + 2*VxB[1]*( ee - c2*bb ) + 7*c2*by*BVxB ),
        Alpha[2] + 2*xi*c2*( 4*ez*EVxB + 2*VxB[2]*( ee - c2*bb ) + 7*c2*bz*BVxB )
    };

    // Calcualting matrix values for the QED correction algorithm

    const amrex::Real a00 = beta + xi*( 8*ex*ex + 14*c2*bx*bx );

    const amrex::Real a11 = beta + xi*( 8*ey*ey + 14*c2*by*by );

    const amrex::Real a22 = beta + xi*( 8*ez*ez + 14*c2*bz*bz );

    const amrex::Real a01 = xi*( 2*ex*ey + 14*c2*bx*by );

    const amrex::Real a02 = xi*( 2*ex*ez + 14*c2*bx*bz );

    const amrex::Real a12 = xi*( 2*ez*ey + 14*c2*bz*by );

    const amrex::Real detA = a00*( a11*a22 - a12*a12 ) - a01*( a01*a22 - a02*a12 )+a02*( a01*a12 - a02*a11 );

    // Calcualting the rows of the inverse matrix using the general 3x3 inverse form

    const amrex::Real invAx[3] = {a22*a11 - a12*a12, a12*a02 - a22*a01, a12*a01 - a11*a02};

    const amrex::Real invAy[3] = {a02*a12 - a22*a01, a00*a22 - a02*a02, a01*a02 - a12*a00};

    const amrex::Real invAz[3] = {a12*a01 - a02*a11, a02*a01 - a12*a00, a11*a00 - a01*a01};

    // Calcualting the final QED corrections by mutliplying the Omega vector with the inverse matrix

    const amrex::Real dEx = (-1/detA)*(invAx[0]*Omega[0] +
                                       invAx[1]*Omega[1] +
                                       invAx[2]*Omega[2]);

    const amrex::Real dEy = (-1/detA)*(invAy[0]*Omega[0] +
                                       invAy[1]*Omega[1] +
                                       invAy[2]*Omega[2]);

    const amrex::Real dEz = (-1/detA)*(invAz[0]*Omega[0] +
                                       invAz[1]*Omega[1] +
                                       invAz[2]*Omega[2]);

   // Adding the QED corrections to the origional fields

    Ex(j,k,l) = Ex(j,k,l) + 0.5*dt*dEx;

    Ey(j,k,l) = Ey(j,k,l) + 0.5*dt*dEy;

    Ez(j,k,l) = Ez(j,k,l) + 0.5*dt*dEz;


// 2D case - follows naturally from 3D case
#else

    // Picking out points for stencil to be used in curl function of M

    amrex::Real Mpx [3] = {0.,0.,0.};
    amrex::Real Mnx [3] = {0.,0.,0.};
    amrex::Real Mpz [3] = {0.,0.,0.};
    amrex::Real Mnz [3] = {0.,0.,0.};

    // Calcualting the M-field at the chosen stencil points

    calc_M(Mpx, tmpEx(j+1,k,0), tmpEy(j+1,k,0), tmpEz(j+1,k,0),
                Bx(j+1,k,0), By(j+1,k,0), Bz(j+1,k,0), xi, c2);
    calc_M(Mnx, tmpEx(j-1,k,0), tmpEy(j-1,k,0), tmpEz(j-1,k,0),
                Bx(j-1,k,0), By(j-1,k,0), Bz(j-1,k,0), xi, c2);
    calc_M(Mpz, tmpEx(j,k+1,0), tmpEy(j,k+1,0), tmpEz(j,k+1,0),
                Bx(j,k+1,0), By(j,k+1,0), Bz(j,k+1,0), xi, c2);
    calc_M(Mnz, tmpEx(j,k-1,0), tmpEy(j,k-1,0), tmpEz(j,k-1,0),
                Bx(j,k-1,0), By(j,k-1,0), Bz(j,k-1,0), xi, c2);

    // Calculating necessary curls

    const amrex::Real VxM[3] = {
        -0.5*(Mpz[1]-Mnz[1])*dzi,
        0.5*( (Mpz[0]-Mnz[0])*dzi - (Mpx[2]-Mnx[2])*dxi ),
        0.5*(Mpx[1]-Mnx[1])*dxi,
    };

    const amrex::Real VxE[3] = {
        -0.5*(tmpEy(j,k+1,0)-tmpEy(j,k-1,0) )*dzi,
        0.5*( (tmpEx(j,k+1,0)-tmpEx(j,k-1,0) )*dzi - (tmpEz(j+1,k,0)-tmpEz(j-1,k,0) )*dxi ),
        0.5*(tmpEy(j+1,k,0)-tmpEy(j-1,k,0) )*dxi,
    };

    const amrex::Real VxB[3] = {
        -0.5*(By(j,k+1,0)-By(j,k-1,0) )*dzi,
        0.5*( (Bx(j,k+1,0)-Bx(j,k-1,0) )*dzi - (Bz(j+1,k,0)-Bz(j-1,k,0) )*dxi ),
        0.5*(By(j+1,k,0)-By(j-1,k,0) )*dxi,
    };

    // Defining comapct values for QED corrections

    const amrex::Real ex = tmpEx(j,k,0);
    const amrex::Real ey = tmpEy(j,k,0);
    const amrex::Real ez = tmpEz(j,k,0);
    const amrex::Real bx = Bx(j,k,0);
    const amrex::Real by = By(j,k,0);
    const amrex::Real bz = Bz(j,k,0);
    const amrex::Real ee = ex*ex + ey*ey + ez*ez;
    const amrex::Real bb = bx*bx + by*by + bz*bz;
    const amrex::Real eb = ex*bx + ey*by + ez*bz;
    const amrex::Real EVxE = ex*VxE[0] + ey*VxE[1] + ez*VxE[2];
    const amrex::Real BVxE = bx*VxE[0] + by*VxE[1] + bz*VxE[2];
    const amrex::Real EVxB = ex*VxB[0] + ey*VxB[1] + ez*VxB[2];
    const amrex::Real BVxB = bx*VxB[0] + by*VxB[1] + bz*VxB[2];

    const amrex::Real beta = 4*xi*( ee - c2*bb ) + PhysConst::ep0;

    const amrex::Real Alpha[3] = {
        2*xi*c2*( -7*bx*EVxE - 7*VxE[0]*eb + 4*ex*BVxE ) + VxM[0],
        2*xi*c2*( -7*by*EVxE - 7*VxE[1]*eb + 4*ey*BVxE ) + VxM[1],
        2*xi*c2*( -7*bz*EVxE - 7*VxE[2]*eb + 4*ez*BVxE ) + VxM[2]
    };

    const amrex::Real Omega[3] = {
        Alpha[0] + 2*xi*c2*( 4*ex*EVxB + 2*VxB[0]*( ee - c2*bb ) + 7*c2*bx*BVxB ),
        Alpha[1] + 2*xi*c2*( 4*ey*EVxB + 2*VxB[1]*( ee - c2*bb ) + 7*c2*by*BVxB ),
        Alpha[2] + 2*xi*c2*( 4*ez*EVxB + 2*VxB[2]*( ee - c2*bb ) + 7*c2*bz*BVxB )
    };

    // Calcualting matrix values for the QED correction algorithm

    const amrex::Real a00 = beta + xi*( 8*ex*ex + 14*c2*bx*bx );

    const amrex::Real a11 = beta + xi*( 8*ey*ey + 14*c2*by*by );

    const amrex::Real a22 = beta + xi*( 8*ez*ez + 14*c2*bz*bz );

    const amrex::Real a01 = xi*( 2*ex*ey + 14*c2*bx*by );

    const amrex::Real a02 = xi*( 2*ex*ez + 14*c2*bx*bz );

    const amrex::Real a12 = xi*( 2*ez*ey + 14*c2*bz*by );

    const amrex::Real detA = a00*( a11*a22 - a12*a12 ) - a01*( a01*a22 - a02*a12 ) + a02*( a01*a12 - a02*a11 );

    // Calcualting inverse matrix values using the general 3x3 form

    const amrex::Real invAx[3] = {a22*a11 - a12*a12, a12*a02 - a22*a01, a12*a01 - a11*a02};

    const amrex::Real invAy[3] = {a02*a12 - a22*a01, a00*a22 - a02*a02, a01*a02 - a12*a00};

    const amrex::Real invAz[3] = {a12*a01 - a02*a11, a02*a01 - a12*a00, a11*a00 - a01*a01};

    // Calcualting the final QED corrections by mutliplying the Omega vector with the inverse matrix

    const amrex::Real dEx = (-1/detA)*(invAx[0]*Omega[0] +
                                       invAx[1]*Omega[1] +
                                       invAx[2]*Omega[2]);

    const amrex::Real dEy = (-1/detA)*(invAy[0]*Omega[0] +
                                       invAy[1]*Omega[1] +
                                       invAy[2]*Omega[2]);

    const amrex::Real dEz = (-1/detA)*(invAz[0]*Omega[0] +
                                       invAz[1]*Omega[1] +
                                       invAz[2]*Omega[2]);

   // Adding the QED corrections to the origional fields

    Ex(j,k,0) = Ex(j,k,0) + 0.5*dt*dEx;

    Ey(j,k,0) = Ey(j,k,0) + 0.5*dt*dEy;

    Ez(j,k,0) = Ez(j,k,0) + 0.5*dt*dEz;

#endif

}

#endif
