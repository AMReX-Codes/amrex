#ifndef MY_TEST_K_H_
#define MY_TEST_K_H_

#include <AMReX_FArrayBox.H>

inline void init (amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real R2,
                  amrex::Real& u, amrex::Real& v, amrex::Real& w,
                  amrex::Real& urhs, amrex::Real& vrhs, amrex::Real& wrhs,
                  amrex::Real& eta)
{
    constexpr amrex::Real pi = 3.1415926535897932;
    amrex::Real x2 = x*x;
    amrex::Real y2 = y*y;
    amrex::Real r2 = x2+y2;
    amrex::Real r2R2 = r2/R2;
    amrex::Real R4 = R2*R2;
    amrex::Real sinx = std::sin(x);
    amrex::Real cosx = std::cos(x);
    amrex::Real siny = std::sin(y);
    amrex::Real cosy = std::cos(y);
    amrex::Real sinz = std::sin(z);
    amrex::Real cosz = std::cos(z);

    amrex::Real f1 = std::sin(pi*r2R2);
    amrex::Real f2 = std::sin(2.*pi*r2R2);
    amrex::Real f3 = std::sin(3.*pi*r2R2);

    u = sinx * (1. + siny) * cosz * f1;
    v = (1. + sinx) * siny * sinz * f2;
    w = cosx * cosy * sinz * f3;
    eta = 2. + sinx * siny * sinz;

    amrex::Real df1dx = std::cos(pi*r2R2) * (2.*pi)/R2 * x;
    amrex::Real df1dy = std::cos(pi*r2R2) * (2.*pi)/R2 * y;
    amrex::Real df2dx = std::cos(2.*pi*r2R2) * (4.*pi)/R2 * x;
    amrex::Real df2dy = std::cos(2.*pi*r2R2) * (4.*pi)/R2 * y;
    amrex::Real df3dx = std::cos(3.*pi*r2R2) * (6.*pi)/R2 * x;
    amrex::Real df3dy = std::cos(3.*pi*r2R2) * (6.*pi)/R2 * y;

    amrex::Real df1dx2 = -std::sin(pi*r2R2) * (4.*pi*pi)/R4 * x2
        + std::cos(pi*r2R2) * (2.*pi)/R2;
    amrex::Real df1dy2 = -std::sin(pi*r2R2) * (4.*pi*pi)/R4 * y2
        + std::cos(pi*r2R2) * (2.*pi)/R2;
    amrex::Real df1dxdy = -std::sin(pi*r2R2) * (4.*pi*pi)/R4 * x * y;

    amrex::Real df2dx2 = -std::sin(2.*pi*r2R2) * (16.*pi*pi)/R4 * x2
        + std::cos(2.*pi*r2R2) * (4.*pi)/R2;
    amrex::Real df2dy2 = -std::sin(2.*pi*r2R2) * (16.*pi*pi)/R4 * y2
        + std::cos(2.*pi*r2R2) * (4.*pi)/R2;
    amrex::Real df2dxdy = -std::sin(2.*pi*r2R2) * (16.*pi*pi)/R4 * x * y;

    amrex::Real df3dx2 = -std::sin(3.*pi*r2R2) * (36.*pi*pi)/R4 * x2
        + std::cos(3.*pi*r2R2) * (6.*pi)/R2;
    amrex::Real df3dxdy = -std::sin(3.*pi*r2R2) * (36.*pi*pi)/R4 * x * y;
    amrex::Real df3dy2 = -std::sin(3.*pi*r2R2) * (36.*pi*pi)/R4 * y2
        + std::cos(3.*pi*r2R2) * (6.*pi)/R2;

    amrex::Real detadx = cosx * siny * sinz;
    amrex::Real detady = sinx * cosy * sinz;
    amrex::Real detadz = sinx * siny * cosz;

    amrex::Real dudx = cosx * (1.+siny) * cosz * f1
        + sinx * (1.+siny) * cosz * df1dx;
    amrex::Real dudy = sinx * cosy * cosz * f1
        + sinx * (1.+siny) * cosz * df1dy;
    amrex::Real dudz = -sinx * (1.+siny) * sinz * f1;

    amrex::Real dvdx = cosx * siny * sinz * f2
        + (1.+sinx) * siny * sinz * df2dx;
    amrex::Real dvdy = (1.+sinx) * cosy * sinz * f2
        + (1.+sinx) * siny * sinz * df2dy;
    amrex::Real dvdz = (1.+sinx) * siny * cosz * f2;

    amrex::Real dwdx = -sinx * cosy * sinz * f3
        + cosx * cosy * sinz * df3dx;
    amrex::Real dwdy = -cosx * siny * sinz * f3
        + cosx * cosy * sinz * df3dy;
    amrex::Real dwdz = cosx * cosy * cosz * f3;

    amrex::Real dudx2 = -sinx * (1.+siny) * cosz * f1
        + cosx * (1.+siny) * cosz * df1dx
        + cosx * (1.+siny) * cosz * df1dx
        + sinx * (1.+siny) * cosz * df1dx2;
    amrex::Real dudy2 = -sinx * siny * cosz * f1
        + sinx * cosy * cosz * df1dy
        + sinx * cosy * cosz * df1dy
        + sinx * (1.+siny) * cosz * df1dy2;
    amrex::Real dudz2 = -sinx * (1.+siny) * cosz * f1;
    amrex::Real dudxdy = cosx * cosy * cosz * f1
        + cosx * (1.+siny) * cosz * df1dy
        + sinx * cosy * cosz * df1dx
        + sinx * (1.+siny) * cosz * df1dxdy;
    amrex::Real dudxdz = -cosx * (1.+siny) * sinz * f1
        - sinx * (1.+siny) * sinz * df1dx;

    amrex::Real dvdx2 = -sinx * siny * sinz * f2
        + cosx * siny * sinz * df2dx
        + cosx * siny * sinz * df2dx
        + (1.+sinx) * siny * sinz * df2dx2;
    amrex::Real dvdy2 = -(1.+sinx) * siny * sinz * f2
        + (1.+sinx) * cosy * sinz * df2dy
        + (1.+sinx) * cosy * sinz * df2dy
        + (1.+sinx) * siny * sinz * df2dy2;
    amrex::Real dvdz2 = -(1.+sinx) * siny * sinz * f2;
    amrex::Real dvdxdy = cosx * cosy * sinz * f2
        + cosx * siny * sinz * df2dy
        + (1.+sinx) * cosy * sinz * df2dx
        + (1.+sinx) * siny * sinz * df2dxdy;
    amrex::Real dvdydz = (1.+sinx) * cosy * cosz * f2
        + (1.+sinx) * siny * cosz * df2dy;

    amrex::Real dwdx2 = -cosx * cosy * sinz * f3
        - sinx * cosy * sinz * df3dx
        - sinx * cosy * sinz * df3dx
        + cosx * cosy * sinz * df3dx2;
    amrex::Real dwdy2 = -cosx * cosy * sinz * f3
        - cosx * siny * sinz * df3dy
        - cosx * siny * sinz * df3dy
        + cosx * cosy * sinz * df3dy2;
    amrex::Real dwdz2 = -cosx * cosy * sinz * f3;
    amrex::Real dwdxdz = -sinx * cosy * cosz * f3
        + cosx * cosy * cosz * df3dx;
    amrex::Real dwdydz = -cosx * siny * cosz * f3
        + cosx * cosy * cosz * df3dy;

    amrex::Real dtauxxdx = detadx*((4./3.)*dudx - (2./3.)*dvdy - (2./3.)*dwdz)
        + eta*((4./3.)*dudx2 - (2./3.)*dvdxdy - (2./3.)*dwdxdz);
    amrex::Real dtauyydy = detady*((4./3.)*dvdy - (2./3.)*dudx - (2./3.)*dwdz)
        + eta*((4./3.)*dvdy2 - (2./3.)*dudxdy - (2./3.)*dwdydz);
    amrex::Real dtauzzdz = detadz*((4./3.)*dwdz - (2./3.)*dudx - (2./3.)*dvdy)
        + eta*((4./3.)*dwdz2 - (2./3.)*dudxdz - (2./3)*dvdydz);
    amrex::Real dtauxydx = detadx*(dvdx + dudy) + eta*(dvdx2 + dudxdy);
    amrex::Real dtauxydy = detady*(dvdx + dudy) + eta*(dvdxdy + dudy2);
    amrex::Real dtauxzdx = detadx*(dwdx + dudz) + eta*(dwdx2 + dudxdz);
    amrex::Real dtauxzdz = detadz*(dwdx + dudz) + eta*(dwdxdz + dudz2);
    amrex::Real dtauyzdy = detady*(dwdy + dvdz) + eta*(dwdy2 + dvdydz);
    amrex::Real dtauyzdz = detadz*(dwdy + dvdz) + eta*(dwdydz + dvdz2);

    urhs = -(dtauxxdx + dtauxydy + dtauxzdz);
    vrhs = -(dtauxydx + dtauyydy + dtauyzdz);
    wrhs = -(dtauxzdx + dtauyzdy + dtauzzdz);
}

#endif
