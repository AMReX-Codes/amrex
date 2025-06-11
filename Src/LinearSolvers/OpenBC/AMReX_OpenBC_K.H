#ifndef AMREX_OPENBC_K_H_
#define AMREX_OPENBC_K_H_

#include <AMReX_OpenBC.H>
#include <AMReX_LOUtil_K.H>

namespace amrex::openbc {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void scale_moments (openbc::Moments::array_type& mom)
{ // p!*q! in the order of 0!*0!, 1!*0!, ..., 7!*0!, 0!*1!, 1!*1!, 2!*1!, ..., 6!*1!, 0!*2!, ..., 0!*7!.
    mom[ 2] *= Real(0.5);
    mom[ 3] *= Real(1./6.);
    mom[ 4] *= Real(1./24.);
    mom[ 5] *= Real(1./120.);
    mom[ 6] *= Real(1./720.);
    mom[ 7] *= Real(1./5040.);
    mom[10] *= Real(0.5);
    mom[11] *= Real(1./6.);
    mom[12] *= Real(1./24.);
    mom[13] *= Real(1./120.);
    mom[14] *= Real(1./720.);
    mom[15] *= Real(0.5);
    mom[16] *= Real(0.5);
    mom[17] *= Real(0.25);
    mom[18] *= Real(1./12.);
    mom[19] *= Real(1./48.);
    mom[20] *= Real(1./240.);
    mom[21] *= Real(1./6.);
    mom[22] *= Real(1./6.);
    mom[23] *= Real(1./12.);
    mom[24] *= Real(1./36.);
    mom[25] *= Real(1./144.);
    mom[26] *= Real(1./24.);
    mom[27] *= Real(1./24.);
    mom[28] *= Real(1./48.);
    mom[29] *= Real(1./144.);
    mom[30] *= Real(1./120.);
    mom[31] *= Real(1./120.);
    mom[32] *= Real(1./240.);
    mom[33] *= Real(1./720.);
    mom[34] *= Real(1./720.);
    mom[35] *= Real(1./5040.);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real block_potential (openbc::Moments const& mom, Real xb, Real yb, Real zb)
{
    constexpr Real oneover4pi = Real(1.)/Real(4.*3.1415926535897932);

    xb -= mom.x;
    yb -= mom.y;
    zb -= mom.z;
    Real ri = Real(1.)/std::sqrt(xb*xb+yb*yb+zb*zb);
    Real ri2 = ri*ri;
    Real ri3 = ri2*ri;
    Real ri4 = ri3*ri;
    Real xr, yr;
    if (mom.face.coordDir() == 0) {
        xr = yb*ri;
        yr = zb*ri;
    } else if (mom.face.coordDir() == 1) {
        xr = xb*ri;
        yr = zb*ri;
    } else {
        xr = xb*ri;
        yr = yb*ri;
    }
    Real xr2 = xr *xr;
    Real xr4 = xr2*xr2;
    Real xr6 = xr4*xr2;
    Real yr2 = yr *yr;
    Real yr4 = yr2*yr2;
    Real yr6 = yr4*yr2;
    Real phi = ri * mom.mom[0]
        + ri2*(xr*mom.mom[1] + yr*mom.mom[8])
        + ri3*((Real(3.) * xr2 - Real(1.)) * mom.mom[2] +
               (Real(3.) * xr * yr       ) * mom.mom[9] +
               (Real(3.) * yr2 - Real(1.)) * mom.mom[15])
        + ri4 * (xr * (Real(15.) * xr2 - Real(9.)) * mom.mom[3] +
                 yr * (Real(15.) * xr2 - Real(3.)) * mom.mom[10] +
                 xr * (Real(15.) * yr2 - Real(3.)) * mom.mom[16] +
                 yr * (Real(15.) * yr2 - Real(9.)) * mom.mom[21])
        + ri4*ri * ((Real(105.) * xr4 - Real(90.) * xr2 + Real(9.)) * mom.mom[4] +
                    (xr * yr * (Real(105.) * xr2 - Real(45.))) * mom.mom[11] +
                    (Real(105.) * xr2 * yr2 - Real(15.) * xr2 - Real(15.) * yr2 + Real(3.)) * mom.mom[17] +
                    (xr * yr * (Real(105.) * yr2 - Real(45.))) * mom.mom[22] +
                    (Real(105.) * yr4 - Real(90.) * yr2 + Real(9.)) * mom.mom[26])
        + ri4*ri2 * (xr * (Real(945.)*xr4 - Real(1050.)*xr2 + Real(225.)) * mom.mom[5] +
                     yr * (Real(945.)*xr4 - Real(630.)*xr2 + Real(45.)) * mom.mom[12] +
                     xr * (Real(945.)*xr2*yr2 - Real(105.)*xr2 - Real(315.)*yr2 + Real(45.)) * mom.mom[18] +
                     yr * (Real(945.)*xr2*yr2 - Real(315.)*xr2 - Real(105.)*yr2 + Real(45.)) * mom.mom[23] +
                     xr * (Real(945.)*yr4 - Real(630.)*yr2 + Real(45.)) * mom.mom[27] +
                     yr * (Real(945.)*yr4 - Real(1050.)*yr2 + Real(225.)) * mom.mom[30])
        + ri4*ri3 * (Real(45.) * (Real(231.)*xr6 - Real(315.)*xr4 + Real(105.)*xr2 - Real(5.)) * mom.mom[6] +
                     Real(315.)*xr*yr * (Real(33.)*xr4 - Real(30.)*xr2 + Real(5.)) * mom.mom[13] +
                     Real(45.) * (Real(231.)*xr4*yr2 - Real(21.)*xr4 - Real(126.)*xr2*yr2 + Real(14.)*xr2 + Real(7.)*yr2 - Real(1.)) * mom.mom[19] +
                     Real(945.)*xr*yr * (Real(11.)*xr2*yr2 - Real(3.)*xr2 - Real(3.)*yr2 + Real(1.)) * mom.mom[24] +
                     Real(45.) * (Real(231.)*xr2*yr4 - Real(126.)*xr2*yr2 + Real(7.)*xr2 - Real(21.)*yr4 + Real(14.)*yr2 - Real(1.)) * mom.mom[28] +
                     Real(315.)*xr*yr * (Real(33.)*yr4 - Real(30.)*yr2 + Real(5.)) * mom.mom[31] +
                     Real(45.) * (Real(231.)*yr6 - Real(315.)*yr4 + Real(105.)*yr2 - Real(5.)) * mom.mom[33])
        + ri4*ri4*(Real(315.)*xr*(Real(429.)*xr6 - Real(693.)*xr4 + Real(315.)*xr2 - Real(35.)) * mom.mom[7] +
                   Real(315.)*yr*(Real(429.)*xr6 - Real(495.)*xr4 + Real(135.)*xr2 - Real(5.)) * mom.mom[14] +
                   Real(315.)*xr*(Real(429.)*xr4*yr2 - Real(33.)*xr4 - Real(330.)*xr2*yr2 + Real(30.)*xr2 + Real(45.)*yr2 - Real(5.)) * mom.mom[20] +
                   Real(945.)*yr*(Real(143.)*xr4*yr2 - Real(33.)*xr4 - Real(66.)*xr2*yr2 + Real(18.)*xr2 + Real(3.)*yr2 - Real(1.)) * mom.mom[25] +
                   Real(945.)*xr*(Real(143.)*xr2*yr4 - Real(66.)*xr2*yr2 + Real(3.)*xr2 - Real(33.)*yr4 + Real(18.)*yr2 - Real(1.)) * mom.mom[29] +
                   Real(315.)*yr*(Real(429.)*xr2*yr4 - Real(330.)*xr2*yr2 + Real(45.)*xr2 - Real(33.)*yr4 + Real(30.)*yr2 - Real(5.)) * mom.mom[32] +
                   Real(315.)*xr*(Real(429.)*yr6 - Real(495.)*yr4 + Real(135.)*yr2 - Real(5.)) * mom.mom[34] +
                   Real(315.)*yr*(Real(429.)*yr6 - Real(693.)*yr4 + Real(315.)*yr2 - Real(35.)) * mom.mom[35]);
    return phi*(-oneover4pi);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void interp_coef (int i, int ii, Real* AMREX_RESTRICT c, int crse_ratio)
{
    static_assert(openbc::P == 3, "openbc::P is assumed to be 3 here");
    Real xint = (static_cast<Real>(ii-i*crse_ratio) + Real(0.5))/static_cast<Real>(crse_ratio);
    constexpr Real x[] = {-3._rt, -2._rt, -1._rt, 0._rt, 1._rt, 2._rt, 3._rt, 4._rt};
    poly_interp_coeff<8>(xint, x, c);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real interpccx (int ii, int j, int k, Array4<Real const> const& phi, int crse_ratio)
{
    int i = amrex::coarsen(ii,crse_ratio);
    Real c[8];
    interp_coef(i,ii,c,crse_ratio);

    Real p = Real(0.);
    for (int n = 0; n < 8; ++n) {
        p += c[n] * phi(i-3+n,j,k);
    }
    return p;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real interpccy (int i, int jj, int k, Array4<Real const> const& phi, int crse_ratio)
{
    int j = amrex::coarsen(jj,crse_ratio);
    Real c[8];
    interp_coef(j,jj,c,crse_ratio);

    Real p = Real(0.);
    for (int n = 0; n < 8; ++n) {
        p += c[n] * phi(i,j-3+n,k);
    }
    return p;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real interpccz (int i, int j, int kk, Array4<Real const> const& phi, int crse_ratio)
{
    int k = amrex::coarsen(kk,crse_ratio);
    Real c[8];
    interp_coef(k,kk,c,crse_ratio);

    Real p = Real(0.);
    for (int n = 0; n < 8; ++n) {
        p += c[n] * phi(i,j,k-3+n);
    }
    return p;
}

}

#endif
