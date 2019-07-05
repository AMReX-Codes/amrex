#include <AMReX_EB2_C.H>

namespace amrex { namespace EB2 {

void build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real const> const& levset,
                  Array4<Real const> const& interx, Array4<Real const> const& intery,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& fcx, Array4<Real> const& fcy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo)
{
    constexpr Real small = 1.e-14;
    const Real dxinv = 1.0/dx[0];
    const Real dyinv = 1.0/dx[1];
    const Box& ndbxg1 = amrex::grow(amrex::surroundingNodes(bx),1);
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( ndbxg1, tbx,
    {
        Box lbx = amrex::grow(amrex::surroundingNodes(bx,0),1);
        auto lo = amrex::max_lbound(tbx, lbx);
        auto hi = amrex::min_ubound(tbx, lbx);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            if (fx(i,j,0) == Type::regular) {
                apx(i,j,0) = 1.0;
                fcx(i,j,0) = 0.0;
            } else if (fx(i,j,0) == Type::covered) {
                apx(i,j,0) = 0.0;
                fcx(i,j,0) = 0.0;
            } else {
                if (levset(i,j,0) < 0.0) {
                    apx(i,j,0) = (interx(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5*apx(i,j,0) - 0.5;
                } else {
                    apx(i,j,0) = 1.0 - (interx(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5 - 0.5*apx(i,j,0);
                }

                if (apx(i,j,0) > 1.0-small) {
                    apx(i,j,0) = 1.0;
                    fcx(i,j,0) = 0.0;
                    fx(i,j,0) = Type::regular;
                } else if (apx(i,j,0) < small) {
                    apx(i,j,0) = 0.0;
                    fcx(i,j,0) = 0.0;
                    fx(i,j,0) = Type::covered;
                }
            }
        }}

        lbx = amrex::grow(amrex::surroundingNodes(bx,1),1);
        lo = amrex::max_lbound(tbx, lbx);
        hi = amrex::min_ubound(tbx, lbx);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            if (fy(i,j,0) == Type::regular) {
                apy(i,j,0) = 1.0;
                fcy(i,j,0) = 0.0;
            } else if (fy(i,j,0) == Type::covered) {
                apy(i,j,0) = 0.0;
                fcy(i,j,0) = 0.0;
            } else {
                if (levset(i,j,0) < 0.0) {
                    apy(i,j,0) = (intery(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5*apy(i,j,0) - 0.5;
                } else {
                    apy(i,j,0) = 1.0 - (intery(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5 - 0.5*apy(i,j,0);
                }

                if (apy(i,j,0) > 1.0-small) {
                    apy(i,j,0) = 1.0;
                    fcy(i,j,0) = 0.0;
                    fy(i,j,0) = Type::regular;
                } else if (apy(i,j,0) < small) {
                    apy(i,j,0) = 0.0;
                    fcy(i,j,0) = 0.0;
                    fy(i,j,0) = Type::covered;
                }
            }
        }}
    });

    const Box& bxg1 = amrex::grow(bx, 1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,0).isSingleValued()) {
            if (fx(i,j,0) == Type::regular and fx(i+1,j,0) == Type::regular and
                fy(i,j,0) == Type::regular and fy(i,j+1,0) == Type::regular)
            {
                cell(i,j,0).setRegular();
            }
            else if (fx(i,j,0) == Type::covered and fx(i+1,j,0) == Type::covered and
                     fy(i,j,0) == Type::covered and fy(i,j+1,0) == Type::covered)
            {
                cell(i,j,0).setCovered();
            }
            else
            {
                int ncuts = 0;
                if (fx(i  ,j  ,0) == Type::irregular) ++ncuts;
                if (fx(i+1,j  ,0) == Type::irregular) ++ncuts;
                if (fy(i  ,j  ,0) == Type::irregular) ++ncuts;
                if (fy(i  ,j+1,0) == Type::irregular) ++ncuts;
                if (ncuts > 2) {
                    amrex::Error("amerx::EB2::build_faces: more than 2 cuts not supported");
                }
            }
        }
    });
}

}}
