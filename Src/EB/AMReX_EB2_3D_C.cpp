#include <AMReX_EB2_C.H>

namespace amrex { namespace EB2 {

namespace {
AMREX_GPU_HOST_DEVICE
void cut_face_2d (Real& areafrac, Real& centx, Real& centy,
                  Real& Sx2, Real& Sy2, Real& Sxy,
                  Real axm, Real axp, Real aym, Real ayp,
                  Real bcx, Real bcy) noexcept
{
    constexpr Real small = 1.e-14;
    constexpr Real tiny = 1.e-15;

    Real apnorm = std::hypot(axm-axp,aym-ayp);
    Real nx = (axm-axp) * (1./apnorm); // pointing to the wall
    Real ny = (aym-ayp) * (1./apnorm);

    Real nxabs = std::abs(nx);
    Real nyabs = std::abs(ny);

    if (nxabs < tiny or nyabs > 1.0-tiny) {
        areafrac = 0.5*(axm+axp);
        if (areafrac > 1.0-small) {
            areafrac = 1.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = Sy2 = 1./12.;
            Sxy = 0.0;
        } else if (areafrac < small) {
            areafrac = 0.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = 0.0;
            Sy2 = 0.0;
            Sxy = 0.0;
        } else {
            centx = 0.0;
            centy = (0.125*(ayp-aym) + ny*0.5*bcy*bcy)/areafrac;
            Sx2 = (1./24.)*(axm+axp);
            Sy2 = (1./24.)*(ayp+aym) + ny*(1./3.)*(bcy*bcy*bcy);
            Sxy = 0.0;
        }
    } else if (nyabs < tiny or nxabs > 1.0-tiny) {
        areafrac = 0.5*(aym+ayp);
        if (areafrac > 1.0-small) {
            areafrac = 1.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = Sy2 = 1./12.;
            Sxy = 0.0;
        } else if (areafrac < small) {
            areafrac = 0.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = 0.0;
            Sy2 = 0.0;
            Sxy = 0.0;
        } else {
            centx = (0.125*(axp-axm) + nx*0.5*bcx*bcx)/areafrac;
            centy = 0.0;
            Sx2 = (1./24.)*(axp+axm) + nx*(1./3.)*(bcx*bcx*bcx);
            Sy2 = (1./24.)*(ayp+aym);
            Sxy = 0.0;
        }
    } else {
        Real signx = (nx > 0.0) ? 1.0 : -1.0;
        Real x_ym = (-0.5 + aym)*signx;
        Real x_yp = (-0.5 + ayp)*signx;
        Real aa = nxabs/ny;
        Real dx = x_ym - x_yp;
        Real dx2 = dx * (x_ym + x_yp);
        Real dx3 = dx * (x_ym*x_ym + x_ym*x_yp + x_yp*x_yp);
        Real dx4 = dx * (x_ym + x_yp) * (x_ym*x_ym + x_yp*x_yp);
        Real af1 = 0.5*(axm+axp) + aa*0.5*dx2;
        centx = 0.125*(axp-axm) + aa*(1./6.)*dx3;
        Sx2 = (1./24.)*(axm+axp) + aa*(1./12.)*dx4;

        Real signy = (ny > 0.0) ? 1.0 : -1.0;
        Real y_xm = (-0.5 + axm)*signy;
        Real y_xp = (-0.5 + axp)*signy;
        aa = nyabs/nx;
        Real dy = y_xm - y_xp;
        Real dy2 = dy * (y_xm + y_xp);
        Real dy3 = dy * (y_xm*y_xm + y_xm*y_xp + y_xp*y_xp);
        Real dy4 = dy * (y_xm + y_xp) * (y_xm*y_xm + y_xp*y_xp);
        Real af2 = 0.5*(aym+ayp) + aa*0.5*dy2;
        centy = (1./8.)*(ayp-aym) + aa*(1./6.)*dy3;
        Sy2 = (1./24.)*(aym+ayp) + aa*(1./12.)*dy4;

        Real S_b = (nxabs < nyabs)
            ? (Sx2 - (1./24.) - signx*(1./6.)*(x_ym*x_ym*x_ym+x_yp*x_yp*x_yp)) / ny
            : (Sy2 - (1./24.) - signy*(1./6.)*(y_xm*y_xm*y_xm+y_xp*y_xp*y_xp)) / nx;
        Sxy = (nxabs < nyabs)
            ? -signy*(1./16.)*dy2 + 0.5*nx*S_b
            : -signx*(1./16.)*dx2 + 0.5*ny*S_b;

        areafrac = 0.5*(af1+af2);
        if (areafrac > 1.0-small) {
            areafrac = 1.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = Sy2 = 1./12.;
            Sxy = 0.0;
        } else if (areafrac < small) {
            areafrac = 0.0;
            centx = 0.0;
            centy = 0.0;
            Sx2 = 0.0;
            Sy2 = 0.0;
            Sxy = 0.0;
        } else {
            centx *= 1./areafrac;
            centy *= 1./areafrac;
            centx = amrex::min(amrex::max(centx,-0.5),0.5);
            centy = amrex::min(amrex::max(centy,-0.5),0.5);
        }
    }
}
}

void build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Type_t> const& fz, Array4<Type_t> const& ex,
                  Array4<Type_t> const& ey, Array4<Type_t> const& ez,
                  Array4<Real const> const& levset, Array4<Real const> const& interx,
                  Array4<Real const> const& intery, Array4<Real const> const& interz,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& apz, Array4<Real> const& fcx,
                  Array4<Real> const& fcy, Array4<Real> const& fcz,
                  Array4<Real> const& m2x, Array4<Real> const& m2y,
                  Array4<Real> const& m2z,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo)
{
    constexpr Real small = 1.e-14;
    const Real dxinv = 1.0/dx[0];
    const Real dyinv = 1.0/dx[1];
    const Real dzinv = 1.0/dx[2];

    const Box& xbx = amrex::grow(amrex::surroundingNodes(bx,0),1);
    AMREX_HOST_DEVICE_FOR_3D ( xbx, i, j, k,
    {
        if (fx(i,j,k) == Type::regular) {
            apx(i,j,k) = 1.0;
            fcx(i,j,k,0) = 0.0;
            fcx(i,j,k,1) = 0.0;
            m2x(i,j,k,0) = (1./12.);
            m2x(i,j,k,1) = (1./12.);
            m2x(i,j,k,2) = 0.0;
        } else if (fx(i,j,k) == Type::covered) {
            apx(i,j,k) = 0.0;
            fcx(i,j,k,0) = 0.0;
            fcx(i,j,k,1) = 0.0;
            m2x(i,j,k,0) = 0.0;
            m2x(i,j,k,1) = 0.0;
            m2x(i,j,k,2) = 0.0;
        } else {
            int ncuts = 0;
            Real bcy = 0.0;
            Real bcz = 0.0;

            Real lym;
            if (ey(i,j,k) == Type::regular) {
                lym = 1.0;
            } else if (ey(i,j,k) == Type::covered) {
                lym = 0.0;
            } else  {
                ++ncuts;
                Real cut = (intery(i,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcy  += cut;
                lym = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
                lym = amrex::min(amrex::max(0.0,lym),1.0);
            }

            Real lyp;
            if (ey(i,j,k+1) == Type::regular) {
                lyp = 1.0;
            } else if (ey(i,j,k+1) == Type::covered) {
                lyp = 0.0;
            } else {
                ++ncuts;
                Real cut = (intery(i,j,k+1)-(problo[1]+j*dx[1]))*dyinv;
                bcy += cut;
                bcz += 1.0;
                lyp = (levset(i,j,k+1) < 0.0) ? cut : 1.0-cut;
                lyp = amrex::min(amrex::max(0.0,lyp),1.0);
            }

            Real lzm;
            if (ez(i,j,k) == Type::regular) {
                lzm = 1.0;
            } else if (ez(i,j,k) == Type::covered) {
                lzm = 0.0;
            } else {
                ++ncuts;
                Real cut = (interz(i,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcz += cut;
                lzm = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
                lzm = amrex::min(amrex::max(0.0,lzm),1.0);
            }

            Real lzp;
            if (ez(i,j+1,k) == Type::regular) {
                lzp = 1.0;
            } else if (ez(i,j+1,k) == Type::covered) {
                lzp = 0.0;
            } else {
                ++ncuts;
                Real cut = (interz(i,j+1,k)-(problo[2]+k*dx[2]))*dzinv;
                bcy += 1.0;
                bcz += cut;
                lzp = (levset(i,j+1,k) < 0.0) ? cut : 1.0-cut;
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncuts <= 2,
                                             "amrex::EB2::build_faces: more than 2 cuts not supported");

            if (lym <= small and lyp <= small and lzm <= small and lzp <= small) {
                apx(i,j,k) = 0.0;
                fcx(i,j,k,0) = 0.0;
                fcx(i,j,k,1) = 0.0;
                m2x(i,j,k,0) = 0.0;
                m2x(i,j,k,1) = 0.0;
                m2x(i,j,k,2) = 0.0;
            } else if (lym == lyp and lzm == lzp) {
                apx(i,j,k) = 1.0;
                fcx(i,j,k,0) = 0.0;
                fcx(i,j,k,1) = 0.0;
                m2x(i,j,k,0) = (1./12.);
                m2x(i,j,k,1) = (1./12.);
                m2x(i,j,k,2) = 0.0;
            } else {
                bcy = 0.5*bcy - 0.5;
                bcz = 0.5*bcz - 0.5;
                cut_face_2d(apx(i,j,k),fcx(i,j,k,0),fcx(i,j,k,1),
                            m2x(i,j,k,0),m2x(i,j,k,1),m2x(i,j,k,2),
                            lzm,lzp,lym,lyp,bcy,bcz);
            }

            if (apx(i,j,k) == 0.0) {
                fx(i,j,k) = Type::covered;
            } else if (apx(i,j,k) == 1.0) {
                fx(i,j,k) = Type::regular;
            }
        }
    });

    const Box& ybx = amrex::grow(amrex::surroundingNodes(bx,1),1);
    AMREX_HOST_DEVICE_FOR_3D ( ybx, i, j, k,
    {
        if (fy(i,j,k) == Type::regular) {
            apy(i,j,k) = 1.0;
            fcy(i,j,k,0) = 0.0;
            fcy(i,j,k,1) = 0.0;
            m2y(i,j,k,0) = 1./12.;
            m2y(i,j,k,1) = 1./12.;
            m2y(i,j,k,2) = 0.0;
        } else if (fy(i,j,k) == Type::covered) {
            apy(i,j,k) = 0.0;
            fcy(i,j,k,0) = 0.0;
            fcy(i,j,k,1) = 0.0;
            m2y(i,j,k,0) = 0.0;
            m2y(i,j,k,1) = 0.0;
            m2y(i,j,k,2) = 0.0;
        } else {
            int ncuts = 0;
            Real bcx = 0.0;
            Real bcz = 0.0;

            Real lxm;
            if (ex(i,j,k) == Type::regular) {
                lxm = 1.0;
            } else if (ex(i,j,k) == Type::covered) {
                lxm = 0.0;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                lxm = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
                lxm = amrex::min(amrex::max(0.0,lxm),1.0);
            }

            Real lxp;
            if (ex(i,j,k+1) == Type::regular) {
                lxp = 1.0;
            } else if (ex(i,j,k+1) == Type::covered) {
                lxp = 0.0;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k+1)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                bcz += 1.0;
                lxp = (levset(i,j,k+1) < 0.0) ? cut : 1.0-cut;
                lxp = amrex::min(amrex::max(0.0,lxp),1.0);
            }

            Real lzm;
            if (ez(i,j,k) == Type::regular) {
                lzm = 1.0;
            } else if (ez(i,j,k) == Type::covered) {
                lzm = 0.0;
            } else {
                ++ncuts;
                Real cut = (interz(i,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcz += cut;
                lzm = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
            }

            Real lzp;
            if (ez(i+1,j,k) == Type::regular) {
                lzp = 1.0;
            } else if (ez(i+1,j,k) == Type::covered) {
                lzp = 0.0;
            } else {
                ++ncuts;
                Real cut = (interz(i+1,j,k)-(problo[2]+k*dx[2]))*dzinv;
                bcx += 1.0;
                bcz += cut;
                lzp = (levset(i+1,j,k) < 0.0) ? cut : 1.0-cut;
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncuts <= 2,
                                             "amrex::EB2::build_faces: more than 2 cuts not supported");

            if (lxm <= small and lxp <= small and lzm <= small and lzp <= small) {
                apy(i,j,k) = 0.0;
                fcy(i,j,k,0) = 0.0;
                fcy(i,j,k,1) = 0.0;
                m2y(i,j,k,0) = 0.0;
                m2y(i,j,k,1) = 0.0;
                m2y(i,j,k,2) = 0.0;
            } else if (lxm == lxp and lzm == lzp) {
                apy(i,j,k) = 1.0;
                fcy(i,j,k,0) = 0.0;
                fcy(i,j,k,1) = 0.0;
                m2y(i,j,k,0) = 1./12.;
                m2y(i,j,k,1) = 1./12.;
                m2y(i,j,k,2) = 0.0;
            } else {
                bcx = 0.5*bcx - 0.5;
                bcz = 0.5*bcz - 0.5;
                cut_face_2d(apy(i,j,k),fcy(i,j,k,0),fcy(i,j,k,1),
                            m2y(i,j,k,0),m2y(i,j,k,1),m2y(i,j,k,2),
                            lzm,lzp,lxm,lxp,bcx,bcz);
            }

            if (apy(i,j,k) == 0.0) {
                fy(i,j,k) = Type::covered;
            } else if (apy(i,j,k) == 1.0) {
                fy(i,j,k) = Type::regular;
            }
        }
    });

    const Box& zbx = amrex::grow(amrex::surroundingNodes(bx,2),1);
    AMREX_HOST_DEVICE_FOR_3D ( zbx, i, j, k,
    {
        if (fz(i,j,k) == Type::regular) {
            apz(i,j,k) = 1.0;
            fcz(i,j,k,0) = 0.0;
            fcz(i,j,k,1) = 0.0;
            m2z(i,j,k,0) = 1./12.;
            m2z(i,j,k,1) = 1./12.;
            m2z(i,j,k,2) = 0.0;
        } else if (fz(i,j,k) == Type::covered) {
            apz(i,j,k) = 0.0;
            fcz(i,j,k,0) = 0.0;
            fcz(i,j,k,1) = 0.0;
            m2z(i,j,k,0) = 0.0;
            m2z(i,j,k,1) = 0.0;
            m2z(i,j,k,2) = 0.0;
        } else {
            int ncuts = 0;
            Real bcx = 0.0;
            Real bcy = 0.0;

            Real lxm;
            if (ex(i,j,k) == Type::regular) {
                lxm = 1.0;
            } else if (ex(i,j,k) == Type::covered) {
                lxm = 0.0;
            } else {
                ++ncuts;
                Real cut = (interx(i,j,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                lxm = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
                lxm = amrex::min(amrex::max(0.0,lxm),1.0);
            }

            Real lxp;
            if (ex(i,j+1,k) == Type::regular) {
                lxp = 1.0;
            } else if (ex(i,j+1,k) == Type::covered) {
                lxp = 0.0;
            } else {
                ++ncuts;
                Real cut = (interx(i,j+1,k)-(problo[0]+i*dx[0]))*dxinv;
                bcx += cut;
                bcy += 1.0;
                lxp = (levset(i,j+1,k) < 0.0) ? cut : 1.0-cut;
                lxp = amrex::min(amrex::max(0.0,lxp),1.0);
            }

            Real lym;
            if (ey(i,j,k) == Type::regular) {
                lym = 1.0;
            } else if (ey(i,j,k) == Type::covered) {
                lym = 0.0;
            } else {
                ++ncuts;
                Real cut = (intery(i,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcy += cut;
                lym = (levset(i,j,k) < 0.0) ? cut : 1.0-cut;
            }

            Real lyp;
            if (ey(i+1,j,k) == Type::regular) {
                lyp = 1.0;
            } else if (ey(i+1,j,k) == Type::covered) {
                lyp = 0.0;
            } else {
                ++ncuts;
                Real cut = (intery(i+1,j,k)-(problo[1]+j*dx[1]))*dyinv;
                bcx += 1.0;
                bcy += cut;
                lyp = (levset(i+1,j,k) < 0.0) ? cut : 1.0-cut;
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncuts <= 2,
                                             "amrex::EB2::build_faces: more than 2 cuts not supported");

            if (lxm <= small and lxp <= small and lym <= small and lyp <= small) {
                apz(i,j,k) = 0.0;
                fcz(i,j,k,0) = 0.0;
                fcz(i,j,k,1) = 0.0;
                m2z(i,j,k,0) = 0.0;
                m2z(i,j,k,1) = 0.0;
                m2z(i,j,k,2) = 0.0;
            } else if (lxm == lxp and lym == lyp) {
                apz(i,j,k) = 1.0;
                fcz(i,j,k,0) = 0.0;
                fcz(i,j,k,1) = 0.0;
                m2z(i,j,k,0) = 1./12.;
                m2z(i,j,k,1) = 1./12.;
                m2z(i,j,k,2) = 0.0;
            } else {
                bcx = 0.5*bcx - 0.5;
                bcy = 0.5*bcy - 0.5;
                cut_face_2d(apz(i,j,k),fcz(i,j,k,0),fcz(i,j,k,1),
                            m2z(i,j,k,0),m2z(i,j,k,1),m2z(i,j,k,2),
                            lym,lyp,lxm,lxp,bcx,bcy);
            }

            if (apz(i,j,k) == 0.0) {
                fz(i,j,k) = Type::covered;
            } else if (apz(i,j,k) == 1.0) {
                fz(i,j,k) = Type::regular;
            }
        }
    });

    const Box& bxg1 = amrex::grow(bx,1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,k).isSingleValued()) {
            if (fx(i,j,k) == Type::covered and fx(i+1,j,k) == Type::covered and
                fy(i,j,k) == Type::covered and fy(i,j+1,k) == Type::covered and
                fz(i,j,k) == Type::covered and fz(i,j,k+1) == Type::covered)
            {
                cell(i,j,k).setCovered();
            }
            else if (fx(i,j,k) == Type::regular and fx(i+1,j,k) == Type::regular and
                     fy(i,j,k) == Type::regular and fy(i,j+1,k) == Type::regular and
                     fz(i,j,k) == Type::regular and fz(i,j,k+1) == Type::regular)
            {
                cell(i,j,k).setRegular();
            }
        }
    });
}

}}
