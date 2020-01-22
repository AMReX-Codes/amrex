#ifndef WarpX_QED_K_h
#define WarpX_QED_K_h

#include <AMReX_FArrayBox.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

/**
 *warpx_hybrid_QED_push uses an FDTD scheme to calculate QED corrections to Maxwell's equations and
 *preforms a half timestep correction to the E-fields
 *\param[in] Ex This function modifies the Ex field at the end
 *\param[in] Ey This function modifies the Ey field at the end
 *\param[in] Ez This function modifies the Ez field at the end
 *\param[in] Bx The QED corrections are non-linear functions of B. However,
 *they do not modify B itself
 *\param[in] By The QED corrections are non-linear functions of B. However,
 *they do not modify B itself
 *\param[in] Bz The QED corrections are non-linear functions of B. However,
 * they do not modify B itself
 *\param[in] tempEx Since the corrections to E at a given node are non-linear functions
 *of the values of E on the surronding nodes, temp arrays are used so that modifications to
 *one node do not influence the corrections to the surrodning nodes
 *\param[in] tempEy Since the corrections to E at a given node are non-linear functions
 *of the values of E on the surronding nodes, temp arrays are used so that modifications to
 *one node do not influence the corrections to the surrodning nodes
 *\param[in] tempEz Since the corrections to E at a given node are non-linear functions
 * of the values of E on the surronding nodes, temp arrays are used so that modifications to
 *one node do not influence the corrections to the surrodning nodes
 *\param[in] dx The x spatial step, used for calculating curls
 *\param[in] dy The y spatial step, used for calculating curls
 *\param[in] dz The z spatial step, used for calulating curls
 *\param[in] dt The temporal step, used for the half push/correction to the E-fields at the end of the function
 */

//AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_hybrid_QED_push (int j, int k, int l, Array4<Real> const& Ex, Array4<Real>
    const& Ey, Array4<Real> const& Ez, Array4<Real> const& Bx, Array4<Real> const& By,
    Array4<Real const> const& Bz, Array4<Real> const& tmpEx, Array4<Real>
    const& tmpEy, Array4<Real> const& tmpEz, Real dx, Real dy, Real dz, Real dt)
{

const amrex::Real c2 = PhysConst::c * PhysConst::c;
const amrex::Real xi = PhysConst::xi;

#if (AMREX_SPACEDIM == 3)

    // Picking out points for stencil to be used in curl function of M
    const amrex::Real Mpx[3] = {
        -2*xi*c2*
        (
            2*Bx(j+1,k,l)*
            (
                pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2)
                - c2*
                (
                    pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)
                )
            )
            - 7*tmpEx(j+1,k,l)*
            (
               tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l)
            )
        ),

        -2*xi*c2*
        (
            2*By(j+1,k,l)*
            (
                pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2)
                - c2*
                (
                    pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)
                )
            )
            - 7*tmpEy(j+1,k,l)*
            (
                tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j+1,k,l)*
            (
                pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2)
                - c2*
                (
                    pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)
                )
            )
            - 7*tmpEz(j+1,k,l)*
            (
                tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l)
            )
        ),
    };

    const amrex::Real Mnx[3] = {
        -2*xi*c2*
        (
            2*Bx(j-1,k,l)*
            (
                pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2)
                - c2*
                (
                    pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)
                )
            )
            - 7*tmpEx(j-1,k,l)*
            (
                tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l)
            )
        ),

        -2*xi*c2*
        (
            2*By(j-1,k,l)*
            (
                pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2)
                - c2*
                (
                    pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)
                )
            )
            - 7*tmpEy(j-1,k,l)*
            (
                tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j-1,k,l)*
            (
                pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2)
                - c2*
                (
                    pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)
                )
            )
            - 7*tmpEz(j-1,k,l)*
            (
                tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l)
            )
        ),
    };

    const amrex::Real Mpy[3] = {
        -2*xi*c2*
        (
            2*Bx(j,k+1,l)*
            (
                pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2)
                - c2*
                (
                    pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)
                )
            )
            - 7*tmpEx(j,k+1,l)*
            (
                tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k+1,l)*
            (
                pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2)
                - c2*
                (
                    pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)
                )
            )
            - 7*tmpEy(j,k+1,l)*
            (
                tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k+1,l)*
            (
                pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2)
                - c2*
                (
                    pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2))
            )
            - 7*tmpEz(j,k+1,l)*
            (
                 tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l)
            )
        ),
    };

    const amrex::Real Mny[3] = {
        -2*xi*c2*
        (
            2*Bx(j,k-1,l)*
            (
                pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2)
                - c2*
                (
                    pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)
                )
            )
            - 7*tmpEx(j,k-1,l)*
            (
                tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k-1,l)*
            (
                pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2)
                - c2*
                (
                    pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)
                )
            )
            - 7*tmpEy(j,k-1,l)*
            (
                tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k-1,l)*
            (
                pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2)
                - c2*
                (
                    pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)
                )
            )
            - 7*tmpEz(j,k-1,l)*
            (
                tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l)
            )
        ),
    };

    const amrex::Real Mpz[3] = {

        -2*xi*c2*
        (
            2*Bx(j,k,l+1)*
            (
                pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2)
                - c2*
                (
                    pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)
                )
            )
            - 7*tmpEx(j,k,l+1)*
            (
                tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k,l+1)*
            (
                pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2)
                - c2*
                (
                    pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)
                )
            )
            - 7*tmpEy(j,k,l+1)*
            (
                tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k,l+1)*
            (
                pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2)
                - c2*
                (
                    pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)
                )
            )
            - 7*tmpEz(j,k,l+1)*
            (
                tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1)
            )
        ),
    };

    const amrex::Real Mnz[3] = {

        -2*xi*c2*
        (
            2*Bx(j,k,l-1)*(
                pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2)
                - c2*
                (
                    pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)
                )
            )
            - 7*tmpEx(j,k,l-1)*
            (
                tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k,l-1)*
            (
                pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2)
                - c2*
                (
                    pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)
                )
            )
            - 7*tmpEy(j,k,l-1)*
            (
                tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k,l-1)*
            (
                pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2)
                - c2*
                (
                    pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)
                )
            )
            - 7*tmpEz(j,k,l-1)*
            (
                tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1)
            )
        ),
    };

    // Calculating necessary curls
    const amrex::Real VxM[3] = {
        0.5*( (Mpy[2]-Mny[2])/dy - (Mpz[1]-Mnz[1])/dz ),
        0.5*( (Mpz[0]-Mnz[0])/dz - (Mpx[2]-Mnx[2])/dx ),
        0.5*( (Mpx[1]-Mnx[1])/dx - (Mpy[0]-Mny[0])/dy ),
    };

    const amrex::Real VxE[3] = {
        0.5*( (tmpEz(j,k+1,l)-tmpEz(j,k-1,l) )/dy - (tmpEy(j,k,l+1)-tmpEy(j,k,l-1) )/dz ),
        0.5*( (tmpEx(j,k,l+1)-tmpEx(j,k,l-1) )/dz - (tmpEz(j+1,k,l)-tmpEz(j-1,k,l) )/dx ),
        0.5*( (tmpEy(j+1,k,l)-tmpEy(j-1,k,l) )/dx - (tmpEx(j,k+1,l)-tmpEx(j,k-1,l) )/dy ),
    };

    const amrex::Real VxB[3] = {
        0.5*( (Bz(j,k+1,l)-Bz(j,k-1,l) )/dy - (By(j,k,l+1)-By(j,k,l-1) )/dz ),
        0.5*( (Bx(j,k,l+1)-Bx(j,k,l-1) )/dz - (Bz(j+1,k,l)-Bz(j-1,k,l) )/dx ),
        0.5*( (By(j+1,k,l)-By(j-1,k,l) )/dx - (Bx(j,k+1,l)-Bx(j,k-1,l) )/dy ),
    };

    // Defining comapct values for QED corrections
    const amrex::Real beta =
        4*(xi)*
        (
            pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
            - c2*
            (
                pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)
            )
        )
        +PhysConst::ep0;

    const amrex::Real Alpha[3] = {

        2*xi*c2*
        (
            (-7)*Bx(j,k,l)*
            (
                tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2]
            )
            - 7*VxE[0]*
            (
                tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l)
            )
            + 4*tmpEx(j,k,l)*
            (
                Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2]
            )
        )
        + VxM[0],

        2*xi*c2*
        (
            (-7)*By(j,k,l)*
            (
                tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2]
            )
            - 7*VxE[1]*
            (
                tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l)
            )
            + 4*tmpEy(j,k,l)*
            (
                Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2]
            )
        )
        + VxM[1],

        2*xi*c2*
        (
            (-7)*Bz(j,k,l)*
            (
                tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2]
            )
            - 7*VxE[2]*
            (
                tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l)
            )
            + 4*tmpEz(j,k,l)*
            (
                Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2]
            )
        )
        + VxM[2],
    };

    const amrex::Real Omega[3] = {

        Alpha[0] + 2*xi*c2*
            (
                4*tmpEx(j,k,l)*
                (
                    tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1]+tmpEz(j,k,l)*VxB[2]
                )
                + 2*VxB[0]*
                (
                    pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
                    - c2*
                    (
                        pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)
                    )
                )
                + 7*c2*Bx(j,k,l)*
                (
                    VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l)
                )
            ),

        Alpha[1] + 2*xi*c2*
            (
                4*tmpEy(j,k,l)*
                (
                    tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1]+tmpEz(j,k,l)*VxB[2]
                )
                + 2*VxB[1]*
                (
                    pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
                    - c2*
                    (
                        pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)
                    )
                )
                + 7*c2*By(j,k,l)*
                (
                    VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l)
                )
            ),

        Alpha[2] + 2*xi*c2*
            (
                4*tmpEz(j,k,l)*
                (
                    tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1]+tmpEz(j,k,l)*VxB[2]
                )
                + 2*VxB[2]*
                (
                    pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
                    - c2*
                    (
                        pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)
                    )
                )
                + 7*c2*Bz(j,k,l)*
                (
                    VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l)
                )
            ),
    };

    // Defining matrix for the algorithm
    const amrex::Real a00 = beta + (xi) *
        ( 8*pow(tmpEx(j,k,l),2) + 14*c2*pow(Bx(j,k,l),2) );

    const amrex::Real a11 = beta +(xi) *
        ( 8*pow(tmpEy(j,k,l),2) + 14*c2*pow(By(j,k,l),2) );

    const amrex::Real a22 = beta + (xi) *
        ( 8*pow(tmpEz(j,k,l),2) + 14*c2*pow(Bz(j,k,l),2) );

    const amrex::Real a01 = (xi) *
         ( 2*tmpEx(j,k,l)*tmpEy(j,k,l) + 14*c2*Bx(j,k,l)*By(j,k,l) );

    const amrex::Real a02 = (xi) *
        ( 2*tmpEx(j,k,l)*tmpEz(j,k,l) + 14*c2*Bx(j,k,l)*Bz(j,k,l) );

    const amrex::Real a12 = (xi/PhysConst::ep0) *
        ( 2*tmpEz(j,k,l)*tmpEy(j,k,l) + 14*c2*Bz(j,k,l)*By(j,k,l) );

    const amrex::Real detA = a00*(a11*a22-pow(a12,2))-a01*(a01*a22-a02*a12)+a02*(a01*a12-a02*a11);

    // Inverting the matrix
    const amrex::Real invAx[3] = {a22*a11-pow(a12,2), a12*a02-a22*a01, a12*a01-a11*a02};
    const amrex::Real invAy[3] = {a02*a12-a22*a01, a00*a22-pow(a02,2), a01*a02-a12*a00};
    const amrex::Real invAz[3] = {a12*a01-a02*a11, a02*a01-a12*a00, a11*a00-pow(a01,2)};

    // Performing push
    Ex(j,k,l) = Ex(j,k,l) - 0.5*dt*(invAx[0]*Omega[0]+invAx[1]*Omega[1]+invAx[2]*Omega[2])/detA;

    Ey(j,k,l) = Ey(j,k,l) - 0.5*dt*(invAy[0]*Omega[0]+invAy[1]*Omega[1]+invAy[2]*Omega[2])/detA;

    Ez(j,k,l) = Ez(j,k,l) - 0.5*dt*(invAz[0]*Omega[0]+invAz[1]*Omega[1]+invAz[2]*Omega[2])/detA;


// 2D case - follows natrually from 3D case
#else
    const amrex::Real Mpx[3] = {

        -2*xi*c2*
        (
            2*Bx(j+1,k,0)*
            (
                pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2)
                - c2*
                (
                    pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)
                )
            )
            - 7*tmpEx(j+1,k,0)*
            (
                tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0)
            )
        ),

        -2*xi*c2*
        (
            2*By(j+1,k,0)*
            (
                pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2)
                - c2*
                (
                    pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)
                )
            )
            - 7*tmpEy(j+1,k,0)*
            (
                tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j+1,k,0)*
            (
                pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2)
                -c2*
                (
                    pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)
                )
            )
            - 7*tmpEz(j+1,k,0)*
            (
                tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0)
            )
        ),
    };

    const amrex::Real Mnx[3] = {

        -2*xi*c2*
        (
            2*Bx(j-1,k,0)*
            (
                pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
                - c2*
                (
                    pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)
                )
            )
            - 7*tmpEx(j-1,k,0)*
            (
                tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0)
            )
        ),

        -2*xi*c2*
        (
            2*By(j-1,k,0)*
            (
                pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
                - c2*
                (
                    pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)
                )
            )
            -7*tmpEy(j-1,k,0)*
            (
                tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j-1,k,0)*
            (
                pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
                - c2*
                (
                    pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)
                )
            )
            - 7*tmpEz(j-1,k,0)*
            (
                tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0)
            )
        ),
    };

     const amrex::Real Mpz[3] = {

        -2*xi*c2*
        (
            2*Bx(j,k+1,0)*
            (
                pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
                - c2*
                (
                    pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)
                )
            )
            - 7*tmpEx(j,k+1,0)*
            (
                tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k+1,0)*
            (
                pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
                - c2*
                (
                    pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)
                )
            )
            - 7*tmpEy(j,k+1,0)*
            (
                tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k+1,0)*
            (
                pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
                - c2*
                (
                    pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)
                )
            )
            - 7*tmpEz(j,k+1,0)*
            (
                tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0)
            )
        ),
     };

    const amrex::Real Mnz[3] = {

        -2*xi*c2*
        (
            2*Bx(j,k-1,0)*
            (
                pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
                - c2*
                (
                    pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)
                )
            )
            - 7*tmpEx(j,k-1,0)*
            (
                tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0)
            )
        ),

        -2*xi*c2*
        (
            2*By(j,k-1,0)*
            (
                pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
                - c2*
                (
                    pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)
                )
            )
            - 7*tmpEy(j,k-1,0)*
            (
                tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0)
            )
        ),

        -2*xi*c2*
        (
            2*Bz(j,k-1,0)*
            (
                pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
                - c2*
                (
                    pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)
                )
            )
            - 7*tmpEz(j,k-1,0)*
            (
                tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0)
            )
        ),
    };

    const amrex::Real VxM[3] = {
        -0.5*(Mpz[1]-Mnz[1])/dz,
        0.5*( (Mpz[0]-Mnz[0])/dz - (Mpx[2]-Mnx[2])/dx ),
        0.5*(Mpx[1]-Mnx[1])/dx,
    };

    const amrex::Real VxE[3] = {
        -0.5*(tmpEy(j,k+1,0)-tmpEy(j,k-1,0) )/dz,
        0.5*( (tmpEx(j,k+1,0)-tmpEx(j,k-1,0) )/dz - (tmpEz(j+1,k,0)-tmpEz(j-1,k,0) )/dx ),
        0.5*(tmpEy(j+1,k,0)-tmpEy(j-1,k,0) )/dx,
    };

    const amrex::Real VxB[3] = {
        -0.5*(By(j,k+1,0)-By(j,k-1,0) )/dz,
        0.5*( (Bx(j,k+1,0)-Bx(j,k-1,0) )/dz - (Bz(j+1,k,0)-Bz(j-1,k,0) )/dx ),
        0.5*(By(j+1,k,0)-By(j-1,k,0) )/dx,
    };

    const amrex::Real beta =
        4*(xi)*
        (
            pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
            - c2*
            (
                pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)
            )
        )
        + PhysConst::ep0;

    const amrex::Real Alpha[3] = {

        2*xi*c2*
            (
                (-7)*Bx(j,k,0)*
                (
                    tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2]
                )
                - 7*VxE[0]*
                (
                    tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0)
                )
                + 4*tmpEx(j,k,0)*
                (
                    Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2]
                )
            )
        + VxM[0],

        2*xi*c2*
            (
                (-7)*By(j,k,0)*
                (
                    tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2]
                )
                - 7*VxE[1]*
                (
                    tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0)
                )
                + 4*tmpEy(j,k,0)*
                (
                    Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2]
                )
            )
            + VxM[1],

        2*xi*c2*
            (
                (-7)*Bz(j,k,0)*
                (
                    tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2]
                )
                - 7*VxE[2]*
                (
                    tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0)
                )
                + 4*tmpEz(j,k,0)*
                (
                    Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2]
                )
            )
            + VxM[2],
    };

    const amrex::Real Omega[3] = {

        Alpha[0] + 2*xi*c2*
            (
                4*tmpEx(j,k,0)*
                (
                    tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2]
                )
                +2*VxB[0]*
                (
                    pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
                    - c2*
                    (
                        pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)
                    )
                )
                + 7*c2*Bx(j,k,0)*
                (
                    VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0)
                )
            ),

        Alpha[1] + 2*xi*c2*
            (
                4*tmpEy(j,k,0)*
                (
                    tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2]
                )
                +2*VxB[1]*
                (
                    pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
                    - c2*
                    (
                        pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)
                    )
                )
                + 7*c2*By(j,k,0)*
                (
                    VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0)
                )
            ),

        Alpha[2] + 2*xi*c2*
            (
                4*tmpEz(j,k,0)*
                (
                    tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2]
                )
                +2*VxB[2]*
                (
                    pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
                    - c2*
                    (
                        pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)
                    )
                )
                + 7*c2*Bz(j,k,0)*
                (
                    VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0)
                )
            ),
    };

    const amrex::Real a00 = beta + (xi) * ( 8*pow(tmpEx(j,k,0),2) + 14*c2*pow(Bx(j,k,0),2) );

    const amrex::Real a11 = beta + (xi) * ( 8*pow(tmpEy(j,k,0),2) + 14*c2*pow(By(j,k,0),2) );

    const amrex::Real a22 = beta + (xi) * ( 8*pow(tmpEz(j,k,0),2) + 14*c2*pow(Bz(j,k,0),2) );

    const amrex::Real a01 = (xi) * ( 2*tmpEx(j,k,0)*tmpEy(j,k,0) + 14*c2*Bx(j,k,0)*By(j,k,0) );

    const amrex::Real a02 = (xi) * ( 2*tmpEx(j,k,0)*tmpEz(j,k,0) + 14*c2*Bx(j,k,0)*Bz(j,k,0) );

    const amrex::Real a12 = (xi) * ( 2*tmpEz(j,k,0)*tmpEy(j,k,0) + 14*c2*Bz(j,k,0)*By(j,k,0) );

    const amrex::Real detA = a00*(a11*a22-pow(a12,2))-a01*(a01*a22-a02*a12)+a02*(a01*a12-a02*a11);

    const amrex::Real invAx[3] = {a22*a11-pow(a12,2), a12*a02-a22*a01, a12*a01-a11*a02};

    const amrex::Real invAy[3] = {a02*a12-a22*a01, a00*a22-pow(a02,2), a01*a02-a12*a00};

    const amrex::Real invAz[3] = {a12*a01-a02*a11, a02*a01-a12*a00, a11*a00-pow(a01,2)};

    const amrex::Real dEx = (-1/detA)*(invAx[0]*Omega[0] +
        invAx[1]*Omega[1] + invAx[2]*Omega[2]);

    const amrex::Real dEy = (-1/detA)*(invAy[0]*Omega[0] +
        invAy[1]*Omega[1] + invAy[2]*Omega[2]);

    const amrex::Real dEz = (-1/detA)*(invAz[0]*Omega[0] +
        invAz[1]*Omega[1] + invAz[2]*Omega[2]);

    Ex(j,k,0) = Ex(j,k,0) + 0.5*dt*dEx;

    Ey(j,k,0) = Ey(j,k,0) + 0.5*dt*dEy;

    Ez(j,k,0) = Ez(j,k,0) + 0.5*dt*dEz;

#endif

}

#endif
