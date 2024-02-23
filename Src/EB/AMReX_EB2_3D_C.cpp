#include <AMReX_EB2_C.H>

namespace amrex::EB2 {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_covered (const int i, const int j, const int k,
                  Array4<EBCellFlag> const& cell,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm) noexcept
{
    vfrac(i,j,k) = 0.0_rt;
    vcent(i,j,k,0) = 0.0_rt;
    vcent(i,j,k,1) = 0.0_rt;
    vcent(i,j,k,2) = 0.0_rt;
    bcent(i,j,k,0) = -1.0_rt;
    bcent(i,j,k,1) = -1.0_rt;
    bcent(i,j,k,2) = -1.0_rt;
    bnorm(i,j,k,0) = 0.0_rt;
    bnorm(i,j,k,1) = 0.0_rt;
    bnorm(i,j,k,2) = 0.0_rt;
    barea(i,j,k) = 0.0_rt;
    cell(i,j,k).setCovered();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_eb_data (const int i, const int j, const int k,
                  Array4<EBCellFlag> const& cell, Array4<Real> const& apx,
                  Array4<Real> const& apy, Array4<Real> const& apz,
                  Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                  Array4<Real const> const& fcz, Array4<Real const> const& m2x,
                  Array4<Real const> const& m2y, Array4<Real const> const& m2z,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Real small_volfrac,
                  bool& is_small_cell, bool& is_multicut) noexcept
{
    Real axm = apx(i,j,k);
    Real axp = apx(i+1,j,k);
    Real aym = apy(i,j,k);
    Real ayp = apy(i,j+1,k);
    Real azm = apz(i,j,k);
    Real azp = apz(i,j,k+1);

    // Check for small cell first
    if (((axm == 0.0_rt && axp == 0.0_rt) &&
         (aym == 0.0_rt && ayp == 0.0_rt) &&
         (azm == 0.0_rt || azp == 0.0_rt)) ||
        ((axm == 0.0_rt && axp == 0.0_rt) &&
         (aym == 0.0_rt || ayp == 0.0_rt) &&
         (azm == 0.0_rt && azp == 0.0_rt)) ||
        ((axm == 0.0_rt || axp == 0.0_rt) &&
         (aym == 0.0_rt && ayp == 0.0_rt) &&
         (azm == 0.0_rt && azp == 0.0_rt))) {
        set_covered(i, j, k, cell, vfrac, vcent, barea, bcent, bnorm);
        is_small_cell = true;
        return;
    }

    // Check for multiple cuts
    // We know there are no multiple cuts on faces by now.
    // We need to check the case that there are two cuts
    // at the opposite corners.
    bool multi_cuts = (axm >= 0.5_rt && axm < 1.0_rt &&
                       axp >= 0.5_rt && axp < 1.0_rt &&
                       aym >= 0.5_rt && aym < 1.0_rt &&
                       ayp >= 0.5_rt && ayp < 1.0_rt &&
                       azm >= 0.5_rt && azm < 1.0_rt &&
                       azp >= 0.5_rt && azp < 1.0_rt);

    if (multi_cuts) {
        set_covered(i, j, k, cell, vfrac, vcent, barea, bcent, bnorm);
        is_multicut = true;
        return;
    }

    Real dapx = axm - axp;
    Real dapy = aym - ayp;
    Real dapz = azm - azp;
    Real apnorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
    if (apnorm == 0.0_rt) {
        bool maybe_multi_cuts = (axm == 0.0_rt && axp == 0.0_rt) ||
                                (aym == 0.0_rt && ayp == 0.0_rt) ||
                                (azm == 0.0_rt && azp == 0.0_rt);
        if (maybe_multi_cuts) {
            set_covered(i, j, k, cell, vfrac, vcent, barea, bcent, bnorm);
            is_multicut = true;
            return;
        } else {
            amrex::Abort("amrex::EB2:build_cells: apnorm==0");
        }
    }
    Real apnorminv = 1.0_rt/apnorm;
    Real nx = dapx * apnorminv;
    Real ny = dapy * apnorminv;
    Real nz = dapz * apnorminv;
    bnorm(i,j,k,0) = nx;
    bnorm(i,j,k,1) = ny;
    bnorm(i,j,k,2) = nz;
    barea(i,j,k) = nx*dapx + ny*dapy + nz*dapz;

    Real aax = 0.5_rt*(axm+axp);
    Real aay = 0.5_rt*(aym+ayp);
    Real aaz = 0.5_rt*(azm+azp);
    Real B0 = aax + aay + aaz;
    Real Bx = -nx*aax + ny*(aym*fcy(i,j,k,0)-ayp*fcy(i,j+1,k,0))
                      + nz*(azm*fcz(i,j,k,0)-azp*fcz(i,j,k+1,0));
    Real By = -ny*aay + nx*(axm*fcx(i,j,k,0)-axp*fcx(i+1,j,k,0))
                      + nz*(azm*fcz(i,j,k,1)-azp*fcz(i,j,k+1,1));
    Real Bz = -nz*aaz + nx*(axm*fcx(i,j,k,1)-axp*fcx(i+1,j,k,1))
                      + ny*(aym*fcy(i,j,k,1)-ayp*fcy(i,j+1,k,1));

    vfrac(i,j,k) = 0.5_rt*(B0 + nx*Bx + ny*By + nz*Bz);

    // remove small cell
    if (vfrac(i,j,k) < small_volfrac) {
        set_covered(i, j, k, cell, vfrac, vcent, barea, bcent, bnorm);
        is_small_cell = true;
        return;
    }

    Real bainv = 1.0_rt/barea(i,j,k);
    bcent(i,j,k,0) = bainv * (Bx + nx*vfrac(i,j,k));
    bcent(i,j,k,1) = bainv * (By + ny*vfrac(i,j,k));
    bcent(i,j,k,2) = bainv * (Bz + nz*vfrac(i,j,k));

    Real b1 = 0.5_rt*(axp-axm) + 0.5_rt*(ayp*fcy(i,j+1,k,0) + aym*fcy(i,j,k,0)) + 0.5_rt*(azp*fcz(i,j,k+1,0) + azm*fcz(i,j,k,0));
    Real b2 = 0.5_rt*(axp*fcx(i+1,j,k,0) + axm*fcx(i,j,k,0)) + 0.5_rt*(ayp-aym) + 0.5_rt*(azp*fcz(i,j,k+1,1) + azm*fcz(i,j,k,1));
    Real b3 = 0.5_rt*(axp*fcx(i+1,j,k,1) + axm*fcx(i,j,k,1)) + 0.5_rt*(ayp*fcy(i,j+1,k,1) + aym*fcy(i,j,k,1)) + 0.5_rt*(azp-azm);
    Real b4 = -nx*0.25_rt*(axp-axm) - ny*(m2y(i,j+1,k,0) - m2y(i,j,k,0)) - nz*(m2z(i,j,k+1,0) - m2z(i,j,k,0));
    Real b5 = -nx*(m2x(i+1,j,k,0) - m2x(i,j,k,0)) - ny*0.25_rt*(ayp-aym) - nz*(m2z(i,j,k+1,1) - m2z(i,j,k,1));
    Real b6 = -nx*(m2x(i+1,j,k,1) - m2x(i,j,k,1)) - ny*(m2y(i,j+1,k,1) - m2y(i,j,k,1)) - nz*0.25_rt*(azp-azm);
    Real b7 = -nx*0.5_rt*(axp*fcx(i+1,j,k,0) + axm*fcx(i,j,k,0)) - ny*0.5_rt*(ayp*fcy(i,j+1,k,0) + aym*fcy(i,j,k,0)) - nz*(m2z(i,j,k+1,2) - m2z(i,j,k,2));
    Real b8 = -nx*0.5_rt*(axp*fcx(i+1,j,k,1) + axm*fcx(i,j,k,1)) - ny*(m2y(i,j+1,k,2) - m2y(i,j,k,2)) - nz*0.5_rt*(azp*fcz(i,j,k+1,0) + azm*fcz(i,j,k,0));
    Real b9 = -nx*(m2x(i+1,j,k,2) - m2x(i,j,k,2)) - ny*0.5_rt*(ayp*fcy(i,j+1,k,1) + aym*fcy(i,j,k,1)) - nz*0.5_rt*(azp*fcz(i,j,k+1,1) + azm*fcz(i,j,k,1));

    Real ny2 = ny*ny;
    Real ny3 = ny2*ny;
    Real ny4 = ny3*ny;
    Real nz2 = nz*nz;
    Real nz3 = nz2*nz;
    Real nz4 = nz3*nz;
    Real nz5 = nz4*nz;

    Real Sx = (5._rt*(b1*(5._rt - 3._rt*ny2) + 2._rt*b4*nx*(5._rt - 3._rt*ny2) +
                   ny*(nx*(b2 + 2._rt*b5*ny) + b7*(6._rt - 4._rt*ny2))) +
               (2._rt*b8*(15._rt - 11._rt*ny2 + ny4) +
                nx*(b3*(5._rt - 2._rt*ny2) - 2._rt*b9*ny*(-5._rt + ny2)))*nz +
               (-22._rt*b7*ny - 2._rt*nx*(15._rt*b4 - 5._rt*b6 + b2*ny) +
                ny2*((16._rt*b4 - 4._rt*(b5 + b6))*nx + 10._rt*b7*ny) +
                b1*(-15._rt + 8._rt*ny2))*nz2 +
               2._rt*(-(b9*nx*ny) + 5._rt*b8*(-2._rt + ny2))*nz3 +
               2._rt*b7*ny*nz4);

    Real Sy = (5._rt*(2._rt*b7*nx*(1._rt + 2._rt*ny2) + b2*(2._rt + 3._rt*ny2) +
                   ny*(b1*nx - 2._rt*b4*(-1._rt + ny2) + b5*(4._rt + 6._rt*ny2))) +
               (2._rt*b9*(5._rt + 9._rt*ny2 + ny4) +
                ny*(2._rt*b8*nx*(4._rt + ny2) + b3*(3._rt + 2._rt*ny2)))*nz +
               (2._rt*b7*nx*(4._rt - 5._rt*ny2) - 8._rt*b2*(-1._rt + ny2) +
                2._rt*ny*(-7._rt*b4 + 8._rt*b5 + 3._rt*b6 - b1*nx +
                       2._rt*(b4 - 4._rt*b5 + b6)*ny2))*nz2 +
               2._rt*(b3*ny + b9*(4._rt - 3._rt*ny2))*nz3 +
               (-8._rt*(b2 + b7*nx) + 4._rt*(b4 - 4._rt*b5 + b6)*ny)*nz4 -
               8._rt*b9*nz5);

    Real Sz = (-2._rt*(b3 + b8*nx + b9*ny)*(-5._rt - 4._rt*ny2 + 4._rt*ny4) +
               (5._rt*(2._rt*b4 + 4._rt*b6 + b1*nx) + (3._rt*b2 + 8._rt*b7*nx)*ny -
                2._rt*(7._rt*b4 - 3._rt*b5 - 8._rt*b6 + b1*nx)*ny2 +
                2._rt*b2*ny3 + 4._rt*(b4 + b5 - 4._rt*b6)*ny4)*nz +
               (b3*(15._rt - 8._rt*ny2) - 6._rt*b9*ny*(-3._rt + ny2) -
                10._rt*b8*nx*(-2._rt + ny2))*nz2 +
               2._rt*(-5._rt*b4 + 15._rt*b6 + (b2 + b7*nx)*ny +
                   2._rt*(b4 + b5 - 4._rt*b6)*ny2)*nz3 + 2._rt*b9*ny*nz4);

    Real den = 1._rt / (10._rt*(5._rt + 4._rt*nz2 - 4._rt*nz4 + 2._rt*ny4*(-2._rt + nz2) +
                                  2._rt*ny2*(2._rt - 3._rt*nz2 + nz4)) * (vfrac(i,j,k)+1.e-30_rt) );

    vcent(i,j,k,0) = Sx * den;
    vcent(i,j,k,1) = Sy * den;
    vcent(i,j,k,2) = Sz * den;


}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cut_face_2d (Real& areafrac, Real& centx, Real& centy,
                  Real& Sx2, Real& Sy2, Real& Sxy,
                  Real axm, Real axp, Real aym, Real ayp,
                  Real bcx, Real bcy) noexcept
{
#ifdef AMREX_USE_FLOAT
    constexpr Real small = 1.e-5_rt;
    constexpr Real tiny  = 1.e-6_rt;
#else
    constexpr Real small = 1.e-14;
    constexpr Real tiny  = 1.e-15;
#endif
    Real apnorm = std::hypot(axm-axp,aym-ayp);
    Real nx = (axm-axp) * (1.0_rt/apnorm); // pointing to the wall
    Real ny = (aym-ayp) * (1.0_rt/apnorm);

    Real nxabs = std::abs(nx);
    Real nyabs = std::abs(ny);

    if (nxabs < tiny || nyabs > 1.0_rt-tiny) {
        areafrac = 0.5_rt*(axm+axp);
        if (areafrac > 1.0_rt-small) {
            areafrac = 1.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = Sy2 = 1.0_rt/12._rt;
            Sxy = 0.0_rt;
        } else if (areafrac < small) {
            areafrac = 0.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = 0.0_rt;
            Sy2 = 0.0_rt;
            Sxy = 0.0_rt;
        } else {
            centx = 0.0_rt;
            centy = (0.125_rt*(ayp-aym) + ny*0.5_rt*bcy*bcy)/areafrac;
            Sx2 = (1.0_rt/24._rt)*(axm+axp);
            Sy2 = (1.0_rt/24._rt)*(ayp+aym) + ny*(1.0_rt/3._rt)*(bcy*bcy*bcy);
            Sxy = 0.0_rt;
        }
    } else if (nyabs < tiny || nxabs > 1.0_rt-tiny) {
        areafrac = 0.5_rt*(aym+ayp);
        if (areafrac > 1.0_rt-small) {
            areafrac = 1.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = Sy2 = 1.0_rt/12._rt;
            Sxy = 0.0_rt;
        } else if (areafrac < small) {
            areafrac = 0.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = 0.0_rt;
            Sy2 = 0.0_rt;
            Sxy = 0.0_rt;
        } else {
            centx = (0.125_rt*(axp-axm) + nx*0.5_rt*bcx*bcx)/areafrac;
            centy = 0.0_rt;
            Sx2 = (1.0_rt/24._rt)*(axp+axm) + nx*(1.0_rt/3._rt)*(bcx*bcx*bcx);
            Sy2 = (1.0_rt/24._rt)*(ayp+aym);
            Sxy = 0.0_rt;
        }
    } else {
        Real signx = (nx > 0.0_rt) ? 1.0_rt : -1.0_rt;
        Real x_ym = (-0.5_rt + aym)*signx;
        Real x_yp = (-0.5_rt + ayp)*signx;
        Real aa = nxabs/ny;
        Real dx = x_ym - x_yp;
        Real dx2 = dx * (x_ym + x_yp);
        Real dx3 = dx * (x_ym*x_ym + x_ym*x_yp + x_yp*x_yp);
        Real dx4 = dx * (x_ym + x_yp) * (x_ym*x_ym + x_yp*x_yp);
        Real af1 = 0.5_rt*(axm+axp) + aa*0.5_rt*dx2;
        centx = 0.125_rt*(axp-axm) + aa*(1.0_rt/6._rt)*dx3;
        Sx2 = (1.0_rt/24._rt)*(axm+axp) + aa*(1.0_rt/12._rt)*dx4;

        Real signy = (ny > 0.0_rt) ? 1.0_rt : -1.0_rt;
        Real y_xm = (-0.5_rt + axm)*signy;
        Real y_xp = (-0.5_rt + axp)*signy;
        aa = nyabs/nx;
        Real dy = y_xm - y_xp;
        Real dy2 = dy * (y_xm + y_xp);
        Real dy3 = dy * (y_xm*y_xm + y_xm*y_xp + y_xp*y_xp);
        Real dy4 = dy * (y_xm + y_xp) * (y_xm*y_xm + y_xp*y_xp);
        Real af2 = 0.5_rt*(aym+ayp) + aa*0.5_rt*dy2;
        centy = (1.0_rt/8._rt)*(ayp-aym) + aa*(1.0_rt/6._rt)*dy3;
        Sy2 = (1.0_rt/24._rt)*(aym+ayp) + aa*(1.0_rt/12._rt)*dy4;

        Real S_b = (nxabs < nyabs)
            ? (Sx2 - (1.0_rt/24._rt) - signx*(1.0_rt/6._rt)*(x_ym*x_ym*x_ym+x_yp*x_yp*x_yp)) / ny
            : (Sy2 - (1.0_rt/24._rt) - signy*(1.0_rt/6._rt)*(y_xm*y_xm*y_xm+y_xp*y_xp*y_xp)) / nx;
        Sxy = (nxabs < nyabs)
            ? -signy*(1.0_rt/16._rt)*dy2 + 0.5_rt*nx*S_b
            : -signx*(1.0_rt/16._rt)*dx2 + 0.5_rt*ny*S_b;

        areafrac = 0.5_rt*(af1+af2);
        if (areafrac > 1.0_rt-small) {
            areafrac = 1.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = Sy2 = 1.0_rt/12._rt;
            Sxy = 0.0_rt;
        } else if (areafrac < small) {
            areafrac = 0.0_rt;
            centx = 0.0_rt;
            centy = 0.0_rt;
            Sx2 = 0.0_rt;
            Sy2 = 0.0_rt;
            Sxy = 0.0_rt;
        } else {
            centx *= 1.0_rt/areafrac;
            centy *= 1.0_rt/areafrac;
            centx = amrex::min(amrex::max(centx,Real(-0.5)),Real(0.5));
            centy = amrex::min(amrex::max(centy,Real(-0.5)),Real(0.5));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_eb_cell (int i, int j, int k,
                  Array4<EBCellFlag> const& cell, Array4<Real> const& apx,
                  Array4<Real> const& apy, Array4<Real> const& apz,
                  Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                  Array4<Real const> const& fcz, Array4<Real const> const& m2x,
                  Array4<Real const> const& m2y, Array4<Real const> const& m2z,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Real small_volfrac,
                  bool& is_small_cell, bool& is_multicut) noexcept
{
    if (cell(i,j,k).isRegular()) {
        vfrac(i,j,k) = 1.0_rt;
        vcent(i,j,k,0) = 0.0_rt;
        vcent(i,j,k,1) = 0.0_rt;
        vcent(i,j,k,2) = 0.0_rt;
        bcent(i,j,k,0) = -1.0_rt;
        bcent(i,j,k,1) = -1.0_rt;
        bcent(i,j,k,2) = -1.0_rt;
        bnorm(i,j,k,0) = 0.0_rt;
        bnorm(i,j,k,1) = 0.0_rt;
        bnorm(i,j,k,2) = 0.0_rt;
        barea(i,j,k) = 0.0_rt;
    } else if (cell(i,j,k).isCovered()) {
        vfrac(i,j,k) = 0.0_rt;
        vcent(i,j,k,0) = 0.0_rt;
        vcent(i,j,k,1) = 0.0_rt;
        vcent(i,j,k,2) = 0.0_rt;
        bcent(i,j,k,0) = -1.0_rt;
        bcent(i,j,k,1) = -1.0_rt;
        bcent(i,j,k,2) = -1.0_rt;
        bnorm(i,j,k,0) = 0.0_rt;
        bnorm(i,j,k,1) = 0.0_rt;
        bnorm(i,j,k,2) = 0.0_rt;
        barea(i,j,k) = 0.0_rt;
    } else {
        set_eb_data(i, j , k, cell, apx, apy, apz, fcx, fcy, fcz, m2x, m2y, m2z,
                    vfrac, vcent, barea, bcent, bnorm, small_volfrac,
                    is_small_cell, is_multicut);
    }
}

}

int build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                 Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                 Array4<Type_t> const& fz, Array4<Type_t const> const& ex,
                 Array4<Type_t const> const& ey, Array4<Type_t const> const& ez,
                 Array4<Real> const& levset, Array4<Real const> const& interx,
                 Array4<Real const> const& intery, Array4<Real const> const& interz,
                 Array4<Real> const& apx, Array4<Real> const& apy,
                 Array4<Real> const& apz, Array4<Real> const& fcx,
                 Array4<Real> const& fcy, Array4<Real> const& fcz,
                 Array4<Real> const& m2x, Array4<Real> const& m2y,
                 Array4<Real> const& m2z,
                 GpuArray<Real,AMREX_SPACEDIM> const& dx,
                 GpuArray<Real,AMREX_SPACEDIM> const& problo,
                 bool cover_multiple_cuts) noexcept
{
    Gpu::Buffer<int> nmulticuts = {0};
    int* hp = nmulticuts.hostData();
    int* dp = nmulticuts.data();

#ifdef AMREX_USE_FLOAT
    constexpr Real small = 1.e-5_rt;
#else
    constexpr Real small = 1.e-14;
#endif
    const Real dxinv = 1.0_rt/dx[0];
    const Real dyinv = 1.0_rt/dx[1];
    const Real dzinv = 1.0_rt/dx[2];

    const Box& xbx = amrex::grow(amrex::surroundingNodes(bx,0),1);
    AMREX_HOST_DEVICE_FOR_3D ( xbx, i, j, k,
    {
        if (fx(i,j,k) == Type::regular) {
            apx(i,j,k) = 1.0_rt;
            fcx(i,j,k,0) = 0.0_rt;
            fcx(i,j,k,1) = 0.0_rt;
            m2x(i,j,k,0) = (1.0_rt/12._rt);
            m2x(i,j,k,1) = (1.0_rt/12._rt);
            m2x(i,j,k,2) = 0.0_rt;
        } else if (fx(i,j,k) == Type::covered) {
            apx(i,j,k) = 0.0_rt;
            fcx(i,j,k,0) = 0.0_rt;
            fcx(i,j,k,1) = 0.0_rt;
            m2x(i,j,k,0) = 0.0_rt;
            m2x(i,j,k,1) = 0.0_rt;
            m2x(i,j,k,2) = 0.0_rt;
        } else {
            int ncuts = 0;
            Real bcy = 0.0_rt;
            Real bcz = 0.0_rt;

            Real lym;
            if (ey(i,j,k) == Type::regular) {
                lym = 1.0_rt;
            } else if (ey(i,j,k) == Type::covered) {
                lym = 0.0_rt;
            } else  {
                ++ncuts;
                Real cut = (intery(i,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcy  += cut;
                lym = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
                lym = amrex::min(amrex::max(Real(0.0),lym),Real(1.0));
            }

            Real lyp;
            if (ey(i,j,k+1) == Type::regular) {
                lyp = 1.0_rt;
            } else if (ey(i,j,k+1) == Type::covered) {
                lyp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (intery(i,j,k+1)-(problo[1]+j*dx[1]))*dyinv;
                bcy += cut;
                bcz += 1.0_rt;
                lyp = (levset(i,j,k+1) < 0.0_rt) ? cut : 1.0_rt-cut;
                lyp = amrex::min(amrex::max(Real(0.0),lyp),Real(1.0));
            }

            Real lzm;
            if (ez(i,j,k) == Type::regular) {
                lzm = 1.0_rt;
            } else if (ez(i,j,k) == Type::covered) {
                lzm = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interz(i,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcz += cut;
                lzm = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
                lzm = amrex::min(amrex::max(Real(0.0),lzm),Real(1.0));
            }

            Real lzp;
            if (ez(i,j+1,k) == Type::regular) {
                lzp = 1.0_rt;
            } else if (ez(i,j+1,k) == Type::covered) {
                lzp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interz(i,j+1,k)-(problo[2]+k*dx[2]))*dzinv;
                bcy += 1.0_rt;
                bcz += cut;
                lzp = (levset(i,j+1,k) < 0.0_rt) ? cut : 1.0_rt-cut;
            }

            if (ncuts > 2) {
                Gpu::Atomic::Add(dp,1);
            }

            if ((ncuts > 2) || (lym <= small && lyp <= small && lzm <= small && lzp <= small)) {
                apx(i,j,k) = 0.0_rt;
                fcx(i,j,k,0) = 0.0_rt;
                fcx(i,j,k,1) = 0.0_rt;
                m2x(i,j,k,0) = 0.0_rt;
                m2x(i,j,k,1) = 0.0_rt;
                m2x(i,j,k,2) = 0.0_rt;
            } else if (lym == lyp && lzm == lzp) {
                apx(i,j,k) = 1.0_rt;
                fcx(i,j,k,0) = 0.0_rt;
                fcx(i,j,k,1) = 0.0_rt;
                m2x(i,j,k,0) = (1.0_rt/12._rt);
                m2x(i,j,k,1) = (1.0_rt/12._rt);
                m2x(i,j,k,2) = 0.0_rt;
            } else {
                bcy = 0.5_rt*bcy - 0.5_rt;
                bcz = 0.5_rt*bcz - 0.5_rt;
                cut_face_2d(apx(i,j,k),fcx(i,j,k,0),fcx(i,j,k,1), // NOLINT(readability-suspicious-call-argument)
                            m2x(i,j,k,0),m2x(i,j,k,1),m2x(i,j,k,2),
                            lzm,lzp,lym,lyp,bcy,bcz);
            }

            if (apx(i,j,k) == 0.0_rt) {
                fx(i,j,k) = Type::covered;
            } else if (apx(i,j,k) == 1.0_rt) {
                fx(i,j,k) = Type::regular;
            }
        }
    });

    const Box& ybx = amrex::grow(amrex::surroundingNodes(bx,1),1);
    AMREX_HOST_DEVICE_FOR_3D ( ybx, i, j, k,
    {
        if (fy(i,j,k) == Type::regular) {
            apy(i,j,k) = 1.0_rt;
            fcy(i,j,k,0) = 0.0_rt;
            fcy(i,j,k,1) = 0.0_rt;
            m2y(i,j,k,0) = 1.0_rt/12._rt;
            m2y(i,j,k,1) = 1.0_rt/12._rt;
            m2y(i,j,k,2) = 0.0_rt;
        } else if (fy(i,j,k) == Type::covered) {
            apy(i,j,k) = 0.0_rt;
            fcy(i,j,k,0) = 0.0_rt;
            fcy(i,j,k,1) = 0.0_rt;
            m2y(i,j,k,0) = 0.0_rt;
            m2y(i,j,k,1) = 0.0_rt;
            m2y(i,j,k,2) = 0.0_rt;
        } else {
            int ncuts = 0;
            Real bcx = 0.0_rt;
            Real bcz = 0.0_rt;

            Real lxm;
            if (ex(i,j,k) == Type::regular) {
                lxm = 1.0_rt;
            } else if (ex(i,j,k) == Type::covered) {
                lxm = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                lxm = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
                lxm = amrex::min(amrex::max(Real(0.0),lxm),Real(1.0));
            }

            Real lxp;
            if (ex(i,j,k+1) == Type::regular) {
                lxp = 1.0_rt;
            } else if (ex(i,j,k+1) == Type::covered) {
                lxp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k+1)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                bcz += 1.0_rt;
                lxp = (levset(i,j,k+1) < 0.0_rt) ? cut : 1.0_rt-cut;
                lxp = amrex::min(amrex::max(Real(0.0),lxp),Real(1.0));
            }

            Real lzm;
            if (ez(i,j,k) == Type::regular) {
                lzm = 1.0_rt;
            } else if (ez(i,j,k) == Type::covered) {
                lzm = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interz(i,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcz += cut;
                lzm = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
            }

            Real lzp;
            if (ez(i+1,j,k) == Type::regular) {
                lzp = 1.0_rt;
            } else if (ez(i+1,j,k) == Type::covered) {
                lzp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interz(i+1,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcx += 1.0_rt;
                bcz += cut;
                lzp = (levset(i+1,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
            }

            if (ncuts > 2) {
                Gpu::Atomic::Add(dp,1);
            }

            if ((ncuts > 2) || (lxm <= small && lxp <= small && lzm <= small && lzp <= small)) {
                apy(i,j,k) = 0.0_rt;
                fcy(i,j,k,0) = 0.0_rt;
                fcy(i,j,k,1) = 0.0_rt;
                m2y(i,j,k,0) = 0.0_rt;
                m2y(i,j,k,1) = 0.0_rt;
                m2y(i,j,k,2) = 0.0_rt;
            } else if (lxm == lxp && lzm == lzp) {
                apy(i,j,k) = 1.0_rt;
                fcy(i,j,k,0) = 0.0_rt;
                fcy(i,j,k,1) = 0.0_rt;
                m2y(i,j,k,0) = 1.0_rt/12._rt;
                m2y(i,j,k,1) = 1.0_rt/12._rt;
                m2y(i,j,k,2) = 0.0_rt;
            } else {
                bcx = 0.5_rt*bcx - 0.5_rt;
                bcz = 0.5_rt*bcz - 0.5_rt;
                cut_face_2d(apy(i,j,k),fcy(i,j,k,0),fcy(i,j,k,1), // NOLINT(readability-suspicious-call-argument)
                            m2y(i,j,k,0),m2y(i,j,k,1),m2y(i,j,k,2),
                            lzm,lzp,lxm,lxp,bcx,bcz);
            }

            if (apy(i,j,k) == 0.0_rt) {
                fy(i,j,k) = Type::covered;
            } else if (apy(i,j,k) == 1.0_rt) {
                fy(i,j,k) = Type::regular;
            }
        }
    });

    const Box& zbx = amrex::grow(amrex::surroundingNodes(bx,2),1);
    AMREX_HOST_DEVICE_FOR_3D ( zbx, i, j, k,
    {
        if (fz(i,j,k) == Type::regular) {
            apz(i,j,k) = 1.0_rt;
            fcz(i,j,k,0) = 0.0_rt;
            fcz(i,j,k,1) = 0.0_rt;
            m2z(i,j,k,0) = 1.0_rt/12._rt;
            m2z(i,j,k,1) = 1.0_rt/12._rt;
            m2z(i,j,k,2) = 0.0_rt;
        } else if (fz(i,j,k) == Type::covered) {
            apz(i,j,k) = 0.0_rt;
            fcz(i,j,k,0) = 0.0_rt;
            fcz(i,j,k,1) = 0.0_rt;
            m2z(i,j,k,0) = 0.0_rt;
            m2z(i,j,k,1) = 0.0_rt;
            m2z(i,j,k,2) = 0.0_rt;
        } else {
            int ncuts = 0;
            Real bcx = 0.0_rt;
            Real bcy = 0.0_rt;

            Real lxm;
            if (ex(i,j,k) == Type::regular) {
                lxm = 1.0_rt;
            } else if (ex(i,j,k) == Type::covered) {
                lxm = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                lxm = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
                lxm = amrex::min(amrex::max(Real(0.0),lxm),Real(1.0));
            }

            Real lxp;
            if (ex(i,j+1,k) == Type::regular) {
                lxp = 1.0_rt;
            } else if (ex(i,j+1,k) == Type::covered) {
                lxp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (interx(i,j+1,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                bcy += 1.0_rt;
                lxp = (levset(i,j+1,k) < 0.0_rt) ? cut : 1.0_rt-cut;
                lxp = amrex::min(amrex::max(Real(0.0),lxp),Real(1.0));
            }

            Real lym;
            if (ey(i,j,k) == Type::regular) {
                lym = 1.0_rt;
            } else if (ey(i,j,k) == Type::covered) {
                lym = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (intery(i,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcy += cut;
                lym = (levset(i,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
            }

            Real lyp;
            if (ey(i+1,j,k) == Type::regular) {
                lyp = 1.0_rt;
            } else if (ey(i+1,j,k) == Type::covered) {
                lyp = 0.0_rt;
            } else {
                ++ncuts;
                Real cut = (intery(i+1,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcx += 1.0_rt;
                bcy += cut;
                lyp = (levset(i+1,j,k) < 0.0_rt) ? cut : 1.0_rt-cut;
            }

            if (ncuts > 2) {
                Gpu::Atomic::Add(dp,1);
            }

            if ((ncuts > 2) || (lxm <= small && lxp <= small && lym <= small && lyp <= small)) {
                apz(i,j,k) = 0.0_rt;
                fcz(i,j,k,0) = 0.0_rt;
                fcz(i,j,k,1) = 0.0_rt;
                m2z(i,j,k,0) = 0.0_rt;
                m2z(i,j,k,1) = 0.0_rt;
                m2z(i,j,k,2) = 0.0_rt;
            } else if (lxm == lxp && lym == lyp) {
                apz(i,j,k) = 1.0_rt;
                fcz(i,j,k,0) = 0.0_rt;
                fcz(i,j,k,1) = 0.0_rt;
                m2z(i,j,k,0) = 1.0_rt/12._rt;
                m2z(i,j,k,1) = 1.0_rt/12._rt;
                m2z(i,j,k,2) = 0.0_rt;
            } else {
                bcx = 0.5_rt*bcx - 0.5_rt;
                bcy = 0.5_rt*bcy - 0.5_rt;
                cut_face_2d(apz(i,j,k),fcz(i,j,k,0),fcz(i,j,k,1), // NOLINT(readability-suspicious-call-argument)
                            m2z(i,j,k,0),m2z(i,j,k,1),m2z(i,j,k,2),
                            lym,lyp,lxm,lxp,bcx,bcy);
            }

            if (apz(i,j,k) == 0.0_rt) {
                fz(i,j,k) = Type::covered;
            } else if (apz(i,j,k) == 1.0_rt) {
                fz(i,j,k) = Type::regular;
            }
        }
    });

    const Box& bxg1 = amrex::grow(bx,1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,k).isSingleValued()) {
            if (fx(i,j,k) == Type::covered && fx(i+1,j,k) == Type::covered &&
                fy(i,j,k) == Type::covered && fy(i,j+1,k) == Type::covered &&
                fz(i,j,k) == Type::covered && fz(i,j,k+1) == Type::covered)
            {
                cell(i,j,k).setCovered();
            }
            else if (fx(i,j,k) == Type::regular && fx(i+1,j,k) == Type::regular &&
                     fy(i,j,k) == Type::regular && fy(i,j+1,k) == Type::regular &&
                     fz(i,j,k) == Type::regular && fz(i,j,k+1) == Type::regular)
            {
                cell(i,j,k).setRegular();
            }
        }
    });

    nmulticuts.copyToHost();

    if (*hp > 0) {
        if (cover_multiple_cuts) {
            Box const& nbxg1 = amrex::surroundingNodes(bxg1);
            AMREX_HOST_DEVICE_FOR_3D(nbxg1, i, j, k,
            {
                if (levset(i,j,k) < Real(0.0)) {
                    bool zero_levset =
                        (xbx.contains(i  ,j-1,k-1)
                         &&        fx(i  ,j-1,k-1) == Type::covered) ||
                        (xbx.contains(i  ,j  ,k-1)
                         &&        fx(i  ,j  ,k-1) == Type::covered) ||
                        (xbx.contains(i  ,j-1,k  )
                         &&        fx(i  ,j-1,k  ) == Type::covered) ||
                        (xbx.contains(i  ,j  ,k  )
                         &&        fx(i  ,j  ,k  ) == Type::covered) ||
                        (ybx.contains(i-1,j  ,k-1)
                         &&        fy(i-1,j  ,k-1) == Type::covered) ||
                        (ybx.contains(i  ,j  ,k-1)
                         &&        fy(i  ,j  ,k-1) == Type::covered) ||
                        (ybx.contains(i-1,j  ,k  )
                         &&        fy(i-1,j  ,k  ) == Type::covered) ||
                        (ybx.contains(i  ,j  ,k  )
                         &&        fy(i  ,j  ,k  ) == Type::covered) ||
                        (zbx.contains(i-1,j-1,k  )
                         &&        fz(i-1,j-1,k  ) == Type::covered) ||
                        (zbx.contains(i  ,j-1,k  )
                         &&        fz(i  ,j-1,k  ) == Type::covered) ||
                        (zbx.contains(i-1,j  ,k  )
                         &&        fz(i-1,j  ,k  ) == Type::covered) ||
                        (zbx.contains(i  ,j  ,k  )
                         &&        fz(i  ,j  ,k  ) == Type::covered);
                    if (zero_levset) {
                        levset(i,j,k) = Real(0.0);
                    }
                }
            });
        } else {
            amrex::Abort("amrex::EB2::build_faces: more than 2 cuts not supported");
        }
    }

    return *hp;
}

void build_cells (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Type_t> const& fz, Array4<Real> const& apx,
                  Array4<Real> const& apy, Array4<Real> const& apz,
                  Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                  Array4<Real const> const& fcz, Array4<Real const> const& m2x,
                  Array4<Real const> const& m2y, Array4<Real const> const& m2z,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Array4<EBCellFlag> const& ctmp,
                  Array4<Real> const& levset, Real small_volfrac, Geometry const& geom,
                  bool extend_domain_face, bool cover_multiple_cuts,
                  int& nsmallcells, int& nmulticuts) noexcept
{
    Gpu::Buffer<int> n_smallcell_multicuts = {0,0};
    int* hp = n_smallcell_multicuts.hostData();
    int* dp = n_smallcell_multicuts.data();

    const Box& bxg1 = amrex::grow(bx,1);
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bxg1, i, j, k,
    {
        bool is_small_cell = false;
        bool is_multicut = false;
        set_eb_cell(i, j, k, cell, apx, apy, apz, fcx, fcy, fcz, m2x, m2y, m2z,
                    vfrac, vcent, barea, bcent, bnorm, small_volfrac,
                    is_small_cell, is_multicut);
        if (is_small_cell) {
            Gpu::Atomic::Add(dp, 1);
        }
        if (is_multicut) {
            Gpu::Atomic::Add(dp+1, 1);
        }
    });

    // set cells in the extended region to covered if the
    // corresponding cell on the domain face is covered
    if(extend_domain_face) {

       Box gdomain = geom.Domain();
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           if (geom.isPeriodic(idim)) {
               gdomain.setSmall(idim, std::min(gdomain.smallEnd(idim), bxg1.smallEnd(idim)));
               gdomain.setBig(idim, std::max(gdomain.bigEnd(idim), bxg1.bigEnd(idim)));
           }
       }

       if (! gdomain.contains(bxg1)) {
          AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
          {
              const auto & dlo = gdomain.loVect();
              const auto & dhi = gdomain.hiVect();

              // find the cell(ii,jj,kk) on the corr. domain face
              // this would have already been set to correct value
              bool in_extended_domain = false;
              int ii = i;
              int jj = j;
              int kk = k;
              if(i < dlo[0]) {
                  in_extended_domain = true;
                  ii = dlo[0];
              }
              else if(i > dhi[0]) {
                  in_extended_domain = true;
                  ii = dhi[0];
              }

              if(j < dlo[1]) {
                  in_extended_domain = true;
                  jj = dlo[1];
              }
              else if(j > dhi[1]) {
                  in_extended_domain = true;
                  jj = dhi[1];
              }

              if(k < dlo[2]) {
                  in_extended_domain = true;
                  kk = dlo[2];
              }
              else if(k > dhi[2]) {
                  in_extended_domain = true;
                  kk = dhi[2];
              }

              // set cell in extendable region to covered if necessary
              if( in_extended_domain && (! cell(i,j,k).isCovered())
                  && cell(ii,jj,kk).isCovered() )
              {
                  Gpu::Atomic::Add(dp, 1);
                  set_covered(i,j,k,cell,vfrac,vcent,barea,bcent,bnorm);
              }
          });
       }
    }

    n_smallcell_multicuts.copyToHost();
    nsmallcells += hp[0];
    nmulticuts  += hp[1];

    Box const& nbxg1 = amrex::surroundingNodes(bxg1);
    Box const& bxg1x = amrex::surroundingNodes(bxg1,0);
    Box const& bxg1y = amrex::surroundingNodes(bxg1,1);
    Box const& bxg1z = amrex::surroundingNodes(bxg1,2);
    AMREX_HOST_DEVICE_FOR_3D(nbxg1, i, j, k,
    {
        if (levset(i,j,k) < Real(0.0)) {
            bool zero_levset =
                (bxg1.contains(i-1,j-1,k-1)
                 &&       cell(i-1,j-1,k-1).isCovered()) ||
                (bxg1.contains(i  ,j-1,k-1)
                 &&       cell(i  ,j-1,k-1).isCovered()) ||
                (bxg1.contains(i-1,j  ,k-1)
                 &&       cell(i-1,j  ,k-1).isCovered()) ||
                (bxg1.contains(i  ,j  ,k-1)
                 &&       cell(i  ,j  ,k-1).isCovered()) ||
                (bxg1.contains(i-1,j-1,k  )
                 &&       cell(i-1,j-1,k  ).isCovered()) ||
                (bxg1.contains(i  ,j-1,k  )
                 &&       cell(i  ,j-1,k  ).isCovered()) ||
                (bxg1.contains(i-1,j  ,k  )
                 &&       cell(i-1,j  ,k  ).isCovered()) ||
                (bxg1.contains(i  ,j  ,k  )
                 &&       cell(i  ,j  ,k  ).isCovered()) ||
                (bxg1x.contains(i  ,j-1,k-1)
                 &&          fx(i  ,j-1,k-1) == Type::covered) ||
                (bxg1x.contains(i  ,j  ,k-1)
                 &&          fx(i  ,j  ,k-1) == Type::covered) ||
                (bxg1x.contains(i  ,j-1,k  )
                 &&          fx(i  ,j-1,k  ) == Type::covered) ||
                (bxg1x.contains(i  ,j  ,k  )
                 &&          fx(i  ,j  ,k  ) == Type::covered) ||
                (bxg1y.contains(i-1,j  ,k-1)
                 &&          fy(i-1,j  ,k-1) == Type::covered) ||
                (bxg1y.contains(i  ,j  ,k-1)
                 &&          fy(i  ,j  ,k-1) == Type::covered) ||
                (bxg1y.contains(i-1,j  ,k  )
                 &&          fy(i-1,j  ,k  ) == Type::covered) ||
                (bxg1y.contains(i  ,j  ,k  )
                 &&          fy(i  ,j  ,k  ) == Type::covered) ||
                (bxg1z.contains(i-1,j-1,k  )
                 &&          fz(i-1,j-1,k  ) == Type::covered) ||
                (bxg1z.contains(i  ,j-1,k  )
                 &&          fz(i  ,j-1,k  ) == Type::covered) ||
                (bxg1z.contains(i-1,j  ,k  )
                 &&          fz(i-1,j  ,k  ) == Type::covered) ||
                (bxg1z.contains(i  ,j  ,k  )
                 &&          fz(i  ,j  ,k  ) == Type::covered);
            if (zero_levset) {
                levset(i,j,k) = Real(0.0);
            }
        }
    });

    if (nsmallcells > 0 || nmulticuts > 0) {
        if (!cover_multiple_cuts && nmulticuts > 0) {
            amrex::Abort("amrex::EB2::build_cells: multi-cuts not supported");
        }
        return;
    } else {
        set_connection_flags(bx, bxg1, cell, ctmp, fx, fy, fz);
    }
}

void set_connection_flags (Box const& bx,
                           Box const& bxg1, Array4<EBCellFlag> const& cell,
                           Array4<EBCellFlag> const& ctmp, Array4<Type_t> const& fx,
                           Array4<Type_t> const& fy, Array4<Type_t> const& fz) noexcept
{
    // Build neighbors.  By default all 26 neighbors are already set.
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,k).isCovered()) {
            cell(i,j,k).setDisconnected();
        } else {
            auto flg = cell(i,j,k);

            if (fx(i,j,k) == Type::covered) {
                flg.setDisconnected(-1,0,0);
            }
            if (fx(i+1,j,k) == Type::covered) {
                flg.setDisconnected(1,0,0);
            }
            if (fy(i,j,k) == Type::covered) {
                flg.setDisconnected(0,-1,0);
            }
            if (fy(i,j+1,k) == Type::covered) {
                flg.setDisconnected(0,1,0);
            }
            if (fz(i,j,k) == Type::covered) {
                flg.setDisconnected(0,0,-1);
            }
            if (fz(i,j,k+1) == Type::covered) {
                flg.setDisconnected(0,0,1);
            }

            // x-y
            if ((fx(i,j,k) == Type::covered || fy(i-1,j,k) == Type::covered) &&
                (fx(i,j-1,k) == Type::covered || fy(i,j,k) == Type::covered))
            {
                flg.setDisconnected(-1,-1,0);
            }

            if ((fx(i+1,j,k) == Type::covered || fy(i+1,j,k) == Type::covered) &&
                (fx(i+1,j-1,k) == Type::covered || fy(i,j,k) == Type::covered))
            {
                flg.setDisconnected(1,-1,0);
            }

            if ((fx(i,j,k) == Type::covered || fy(i-1,j+1,k) == Type::covered) &&
                (fx(i,j+1,k) == Type::covered || fy(i,j+1,k) == Type::covered))
            {
                flg.setDisconnected(-1,1,0);
            }

            if ((fx(i+1,j,k) == Type::covered || fy(i+1,j+1,k) == Type::covered) &&
                (fx(i+1,j+1,k) == Type::covered || fy(i,j+1,k) == Type::covered))
            {
                flg.setDisconnected(1,1,0);
            }

            // x-z
            if ((fx(i,j,k) == Type::covered || fz(i-1,j,k) == Type::covered) &&
                (fx(i,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
            {
                flg.setDisconnected(-1,0,-1);
            }

            if ((fx(i+1,j,k) == Type::covered || fz(i+1,j,k) == Type::covered) &&
                (fx(i+1,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
            {
                flg.setDisconnected(1,0,-1);
            }

            if ((fx(i,j,k) == Type::covered || fz(i-1,j,k+1) == Type::covered) &&
                (fx(i,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
            {
                flg.setDisconnected(-1,0,1);
            }

            if ((fx(i+1,j,k) == Type::covered || fz(i+1,j,k+1) == Type::covered) &&
                (fx(i+1,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
            {
                flg.setDisconnected(1,0,1);
            }

            // y-z
            if ((fy(i,j,k) == Type::covered || fz(i,j-1,k) == Type::covered) &&
                (fy(i,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
            {
                flg.setDisconnected(0,-1,-1);
            }

            if ((fy(i,j+1,k) == Type::covered || fz(i,j+1,k) == Type::covered) &&
                (fy(i,j+1,k-1) == Type::covered || fz(i,j,k) == Type::covered))
            {
                flg.setDisconnected(0,1,-1);
            }

            if ((fy(i,j,k) == Type::covered || fz(i,j-1,k+1) == Type::covered) &&
                (fy(i,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
            {
                flg.setDisconnected(0,-1,1);
            }

            if ((fy(i,j+1,k) == Type::covered || fz(i,j+1,k+1) == Type::covered) &&
                (fy(i,j+1,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
            {
                flg.setDisconnected(0,1,1);
            }

            cell(i,j,k) = flg;
        }

        ctmp(i,j,k) = cell(i,j,k);
    });

    AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
    {
        if (!cell(i,j,k).isCovered()) {
            auto tmpflg = ctmp(i,j,k);
            auto newflg = tmpflg;

            // -1, -1, -1 corner
            if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0,-1,-1)) &&
                (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected(-1, 0,-1)) &&
                (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected(-1,-1, 0)))
            {
                newflg.setDisconnected(-1,-1,-1);
            }

            // 1, -1, -1 corner
            if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0,-1,-1)) &&
                (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected( 1, 0,-1)) &&
                (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected( 1,-1, 0)))
            {
                newflg.setDisconnected(1,-1,-1);
            }

            // -1, 1, -1 corner
            if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0, 1,-1)) &&
                (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected(-1, 0,-1)) &&
                (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected(-1, 1, 0)))
            {
                newflg.setDisconnected(-1, 1,-1);
            }

            // 1, 1, -1 corner
            if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0, 1,-1)) &&
                (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected( 1, 0,-1)) &&
                (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected( 1, 1, 0)))
            {
                newflg.setDisconnected(1, 1,-1);
            }

            // -1, -1, 1 corner
            if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0,-1, 1)) &&
                (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected(-1, 0, 1)) &&
                (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected(-1,-1, 0)))
            {
                newflg.setDisconnected(-1,-1, 1);
            }

            // 1, -1, 1 corner
            if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0,-1, 1)) &&
                (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected( 1, 0, 1)) &&
                (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected( 1,-1, 0)))
            {
                newflg.setDisconnected(1,-1, 1);
            }

            // -1, 1, 1 corner
            if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0, 1, 1)) &&
                (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected(-1, 0, 1)) &&
                (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected(-1, 1, 0)))
            {
                newflg.setDisconnected(-1,1,1);
            }

            // 1, 1, 1 corner
            if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0, 1, 1)) &&
                (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected( 1, 0, 1)) &&
                (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected( 1, 1, 0)))
            {
                newflg.setDisconnected(1,1,1);
            }

            cell(i,j,k) = newflg;
        }
    });
}

}
