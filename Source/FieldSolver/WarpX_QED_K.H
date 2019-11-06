#ifndef WarpX_QED_K_h
#define WarpX_QED_K_h

#include <AMREX_FArrayBox.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

// AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_hybrid_QED_push (int j, int k, int l, Array4<Real> const& Ex, Array4<Real>
    const& Ey, Array4<Real> const& Ez, Array4<Real> const& Bx, Array4<Real> const& By,
    Array4<Real const> const& Bz, Array4<Real> const& tmpEx, Array4<Real>
    const& tmpEy, Array4<Real> const& tmpEz, Real dtsdx, Real dtsdy, Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)

    // Picking out points for stencil to be used in curl function
    const amrex::Real Mpx[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2) - pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*Bx(j+1,k,l) -7*(tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l))*tmpEx(j+1,k,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2) -pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*By(j+1,k,l) -7*(tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l))*tmpEy(j+1,k,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,l),2)+pow(tmpEy(j+1,k,l),2)+pow(tmpEz(j+1,k,l),2) -pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*Bz(j+1,k,l) -7*(tmpEx(j+1,k,l)*Bx(j+1,k,l)+tmpEy(j+1,k,l)*By(j+1,k,l)+tmpEz(j+1,k,l)*Bz(j+1,k,l))*tmpEz(j+1,k,l))};

    const amrex::Real Mnx[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2) - pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*Bx(j-1,k,l) -7*(tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l))*tmpEx(j-1,k,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2) - pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*By(j-1,k,l) -7*(tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l))*tmpEy(j-1,k,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,l),2)+pow(tmpEy(j-1,k,l),2)+pow(tmpEz(j-1,k,l),2) - pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*Bz(j-1,k,l) -7*(tmpEx(j-1,k,l)*Bx(j-1,k,l)+tmpEy(j-1,k,l)*By(j-1,k,l)+tmpEz(j-1,k,l)*Bz(j-1,k,l))*tmpEz(j-1,k,l))};

    const amrex::Real Mpy[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2) - pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*Bx(j,k+1,l)
        -7*(tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l))*tmpEx(j,k+1,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*By(j,k+1,l)
        -7*(tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l))*tmpEy(j,k+1,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,l),2)+pow(tmpEy(j,k+1,l),2)+pow(tmpEz(j,k+1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*Bz(j,k+1,l)
        -7*(tmpEx(j,k+1,l)*Bx(j,k+1,l)+tmpEy(j,k+1,l)*By(j,k+1,l)+tmpEz(j,k+1,l)*Bz(j,k+1,l))*tmpEz(j,k+1,l))};

    const amrex::Real Mny[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*Bx(j,k-1,l)
        -7*(tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l))*tmpEx(j,k-1,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*By(j,k-1,l)
        -7*(tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l))*tmpEy(j,k-1,l)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,l),2)+pow(tmpEy(j,k-1,l),2)+pow(tmpEz(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*Bz(j,k-1,l)
        -7*(tmpEx(j,k-1,l)*Bx(j,k-1,l)+tmpEy(j,k-1,l)*By(j,k-1,l)+tmpEz(j,k-1,l)*Bz(j,k-1,l))*tmpEz(j,k-1,l))};

    const amrex::Real Mpz[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*Bx(j,k,l+1)
        -7*(tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1))*tmpEx(j,k,l+1)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*By(j,k,l+1)
        -7*(tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1))*tmpEy(j,k,l+1)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l+1),2)+pow(tmpEy(j,k,l+1),2)+pow(tmpEz(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*Bz(j,k,l+1)
        -7*(tmpEx(j,k,l+1)*Bx(j,k,l+1)+tmpEy(j,k,l+1)*By(j,k,l+1)+tmpEz(j,k,l+1)*Bz(j,k,l+1))*tmpEz(j,k,l+1))};

    const amrex::Real Mnz[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*Bx(j,k,l-1)
        -7*(tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1))*tmpEx(j,k,l-1)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*By(j,k,l-1)
        -7*(tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1))*tmpEy(j,k,l-1)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k,l-1),2)+pow(tmpEy(j,k,l-1),2)+pow(tmpEz(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*Bz(j,k,l-1)
        -7*(tmpEx(j,k,l-1)*Bx(j,k,l-1)+tmpEy(j,k,l-1)*By(j,k,l-1)+tmpEz(j,k,l-1)*Bz(j,k,l-1))*tmpEz(j,k,l-1))};

    // Calculating necessary curls
    const amrex::Real VxM[3] = {
        0.5*(dtsdy*(Mpz[2]-Mnz[2]) - dtsdz*( Mpy[1]-Mny[1] ) ),
        0.5*( dtsdz*(Mpx[0]-Mnx[0]) - dtsdx*( Mpz[2]-Mnz[2] ) ),
        0.5*( dtsdx*(Mpy[1]-Mny[1]) - dtsdy*( Mpx[0]-Mnx[0] ) )};

    const amrex::Real VxE[3] = {
        0.5*(dtsdy*(tmpEz(j,k+1,l)-tmpEz(j,k-1,l)) - dtsdz*(tmpEy(j,k,l+1)-tmpEy(j,k,l-1))),
        0.5*(dtsdz*(tmpEx(j,k,l+1)-tmpEx(j,k,l-1)) - dtsdx*(tmpEz(j+1,k,l)-tmpEz(j-1,k,l))),
        0.5*(dtsdx*(tmpEy(j+1,k,l)-tmpEy(j-1,k,l)) - dtsdy*(tmpEx(j,k+1,l)-tmpEx(j,k-1,l)))};

    const amrex::Real VxB[3] = {
        0.5*(dtsdy*(Bz(j,k+1,l)-Bz(j,k-1,l)) - dtsdz*(By(j,k,l+1)-By(j,k,l-1))),
        0.5*(dtsdz*(Bx(j,k,l+1)-Bx(j,k,l-1)) - dtsdx*(Bz(j+1,k,l)-Bz(j-1,k,l))),
        0.5*(dtsdx*(By(j+1,k,l)-By(j-1,k,l)) - dtsdy*(Bx(j,k+1,l)-Bx(j,k-1,l)))};

    // Defining comapct values for QED corrections
    const amrex::Real beta = 4*(PhysConst::xi/PhysConst::ep0)*(pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))+PhysConst::ep0;

    const amrex::Real Alpha[3] = {

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2])*Bx(j,k,l)
            -7*(tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l))*VxE[0]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*tmpEx(j,k,l)) + VxM[0],

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2])*By(j,k,l)
            -7*(tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l))*VxE[1]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*tmpEy(j,k,l)) + VxM[1],

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,l)*VxE[0]+tmpEy(j,k,l)*VxE[1]+tmpEz(j,k,l)*VxE[2])*Bz(j,k,l)
            -7*(tmpEx(j,k,l)*Bx(j,k,l)+tmpEy(j,k,l)*By(j,k,l)+tmpEz(j,k,l)*Bz(j,k,l))*VxE[2]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*tmpEz(j,k,l)) + VxM[2]};

    const amrex::Real Omega[3] = {

        Alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1]+tmpEz(j,k,l)*VxB[2])*tmpEx(j,k,l)
            +2*(pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[0]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*Bx(j,k,l)),

        Alpha[1] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1]+tmpEz(j,k,l)*VxB[2])*tmpEy(j,k,l)
            +2*(pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[1]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*By(j,k,l)),

        Alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(tmpEx(j,k,l)*VxB[0]+tmpEy(j,k,l)*VxB[1] +
            tmpEz(j,k,l)*VxB[2])*tmpEz(j,k,l)+2*(pow(tmpEx(j,k,l),2)+pow(tmpEy(j,k,l),2)+pow(tmpEz(j,k,l),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[2]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*Bz(j,k,l))};

    // Defining matrix for the algorithm
    const amrex::Real a00 = beta + (PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(tmpEx(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(Bx(j,k,l),2) );

    const amrex::Real a11 = beta +(PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(tmpEy(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(By(j,k,l),2) );

    const amrex::Real a22 = beta + (PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(tmpEz(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(Bz(j,k,l),2) );

    const amrex::Real a01 = (PhysConst::xi/PhysConst::ep0) *
         ( 8*tmpEx(j,k,l)*tmpEy(j,k,l) + 14*pow(PhysConst::c,2)*Bx(j,k,l)*By(j,k,l) );

    const amrex::Real a02 = (PhysConst::xi/PhysConst::ep0) *
        ( 8*tmpEx(j,k,l)*tmpEz(j,k,l) + 14*pow(PhysConst::c,2)*Bx(j,k,l)*Bz(j,k,l) );

    const amrex::Real a12 = (PhysConst::xi/PhysConst::ep0) *
        ( 8*tmpEz(j,k,l)*tmpEy(j,k,l) + 14*pow(PhysConst::c,2)*Bz(j,k,l)*By(j,k,l) );

    const amrex::Real detA = a00*(a11*a22-pow(a12,2))-a01*(a01*a22-a02*a12)+a02*(a01*a12-a02*a11);

    // Inverting the matrix
    const amrex::Real invAx[3] = {a22*a11-pow(a12,2), a12*a02-a22*a01, a12*a01-a11*a02};
    const amrex::Real invAy[3] = {a02*a12-a22*a01, a00*a22-pow(a02,2), a01*a02-a12*a00};
    const amrex::Real invAz[3] = {a12*a01-a02*a11, a02*a01-a12*a00, a11*a00-pow(a01,2)};
    
    // Performing push
    Ex(j,k,l) = Ex(j,k,l) - (invAx[0]*Omega[0]+
        invAx[1]*Omega[1]+invAx[2]*Omega[2])/detA;
    
    Ey(j,k,l) = Ey(j,k,l) - (invAy[0]*Omega[0]+
        invAy[1]*Omega[1]+invAy[2]*Omega[2])/detA;
        
    Ez(j,k,l) = Ez(j,k,l) - (invAz[0]*Omega[0]+
        invAz[1]*Omega[1]+invAz[2]*Omega[2])/detA;



#else
    const amrex::Real Mpx[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2) -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*Bx(j+1,k,0) -7*(tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0))*tmpEx(j+1,k,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2) -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*By(j+1,k,0) -7*(tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0))*tmpEy(j+1,k,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j+1,k,0),2)+pow(tmpEy(j+1,k,0),2)+pow(tmpEz(j+1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*Bz(j+1,k,0)
        -7*(tmpEx(j+1,k,0)*Bx(j+1,k,0)+tmpEy(j+1,k,0)*By(j+1,k,0)+tmpEz(j+1,k,0)*Bz(j+1,k,0))*tmpEz(j+1,k,0))};

    const amrex::Real Mnx[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*Bx(j-1,k,0)
        -7*(tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0))*tmpEx(j-1,k,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*By(j-1,k,0)
        -7*(tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0))*tmpEy(j-1,k,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j-1,k,0),2)+pow(tmpEy(j-1,k,0),2)+pow(tmpEz(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*Bz(j-1,k,0)
        -7*(tmpEx(j-1,k,0)*Bx(j-1,k,0)+tmpEy(j-1,k,0)*By(j-1,k,0)+tmpEz(j-1,k,0)*Bz(j-1,k,0))*tmpEz(j-1,k,0))};

     const amrex::Real Mpz[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*Bx(j,k+1,0)
        -7*(tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0))*tmpEx(j,k+1,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*By(j,k+1,0)
        -7*(tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0))*tmpEy(j,k+1,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k+1,0),2)+pow(tmpEy(j,k+1,0),2)+pow(tmpEz(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*Bz(j,k+1,0)
        -7*(tmpEx(j,k+1,0)*Bx(j,k+1,0)+tmpEy(j,k+1,0)*By(j,k+1,0)+tmpEz(j,k+1,0)*Bz(j,k+1,0))*tmpEz(j,k+1,0))};

    const amrex::Real Mnz[3] = {

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*Bx(j,k-1,0)
        -7*(tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0))*tmpEx(j,k-1,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*By(j,k-1,0)
        -7*(tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0))*tmpEy(j,k-1,0)),

        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(tmpEx(j,k-1,0),2)+pow(tmpEy(j,k-1,0),2)+pow(tmpEz(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*Bz(j,k-1,0)
        -7*(tmpEx(j,k-1,0)*Bx(j,k-1,0)+tmpEy(j,k-1,0)*By(j,k-1,0)+tmpEz(j,k-1,0)*Bz(j,k-1,0))*tmpEz(j,k-1,0))};

    const amrex::Real VxM[3] = {
        -0.5*dtsdz*( Mpz[1]-Mnz[1] ),
        0.5*( dtsdz*(Mpz[0]-Mnz[0]) - dtsdx*(Mpx[2]-Mnx[2]) ),
        0.5*dtsdx*(Mpx[1]-Mnx[1]) };

    const amrex::Real VxE[3] = {
        -0.5*dtsdz*(tmpEy(j,k+1,0)-tmpEy(j,k-1,0)),
        0.5*(dtsdz*(tmpEx(j,k+1,0)-tmpEx(j,k-1,0)) - dtsdx*(tmpEz(j+1,k,0)-tmpEz(j-1,k,0))),
        0.5*dtsdx*(tmpEy(j+1,k,0)-tmpEy(j-1,k,0))};

    const amrex::Real VxB[3] = {
        -0.5*dtsdz*(By(j,k+1,0)-By(j,k-1,0)),
        0.5*(dtsdz*(Bx(j,k+1,0)-Bx(j,k-1,0)) - dtsdx*(Bz(j+1,k,0)-Bz(j-1,k,0))),
        0.5*dtsdx*(By(j+1,k,0)-By(j-1,k,0))};

    const amrex::Real beta =
        4*(PhysConst::xi/PhysConst::ep0)*(pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)-
        pow(PhysConst::c,2)*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))+PhysConst::ep0;

    const amrex::Real Alpha[3] = {

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2])*Bx(j,k,0)
            - 7*(tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0))*VxE[0]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*tmpEx(j,k,0)) + VxM[0],

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2])*By(j,k,0)
            - 7*(tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0))*VxE[1]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*tmpEy(j,k,0)) + VxM[1],

        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(tmpEx(j,k,0)*VxE[0]+tmpEy(j,k,0)*VxE[1]+tmpEz(j,k,0)*VxE[2])*Bz(j,k,0)
            - 7*(tmpEx(j,k,0)*Bx(j,k,0)+tmpEy(j,k,0)*By(j,k,0)+tmpEz(j,k,0)*Bz(j,k,0))*VxE[2]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*tmpEz(j,k,0)) + VxM[2]};

    const amrex::Real Omega[3] = {

        Alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2])*tmpEx(j,k,0)
            +2*(pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[0]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*Bx(j,k,0)),

        Alpha[1] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2])*tmpEy(j,k,0)
            +2*(pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[1]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*By(j,k,0)),

        Alpha[2] + 2*PhysConst::xi*pow(PhysConst::c,2)*
            (4*(tmpEx(j,k,0)*VxB[0]+tmpEy(j,k,0)*VxB[1]+tmpEz(j,k,0)*VxB[2])*tmpEz(j,k,0)
            +2*(pow(tmpEx(j,k,0),2)+pow(tmpEy(j,k,0),2)+pow(tmpEz(j,k,0),2)
            -pow(PhysConst::c,2)*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[2]
            +7*pow(PhysConst::c,2)*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*Bz(j,k,0))};

    const amrex::Real a00 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(tmpEx(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(Bx(j,k,0),2) );
    
    const amrex::Real a11 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(tmpEy(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(By(j,k,0),2) );
    
    const amrex::Real a22 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(tmpEz(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(Bz(j,k,0),2) );
    
    const amrex::Real a01 = (PhysConst::xi/PhysConst::ep0) * ( 8*tmpEx(j,k,0)*tmpEy(j,k,0) + 14*pow(PhysConst::c,2)*Bx(j,k,0)*By(j,k,0) );
    
    const amrex::Real a02 = (PhysConst::xi/PhysConst::ep0) * ( 8*tmpEx(j,k,0)*tmpEz(j,k,0) + 14*pow(PhysConst::c,2)*Bx(j,k,0)*Bz(j,k,0) );
    
    const amrex::Real a12 = (PhysConst::xi/PhysConst::ep0) * ( 8*tmpEz(j,k,0)*tmpEy(j,k,0) + 14*pow(PhysConst::c,2)*Bz(j,k,0)*By(j,k,0) );

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

    Ex(j,k,0) = Ex(j,k,0) + 0.5*dEx;
    
    Ey(j,k,0) = Ey(j,k,0) + 0.5*dEy;
        
    Ez(j,k,0) = Ez(j,k,0) + 0.5*dEz;

#endif

}

#endif
