#ifndef AMREX_EB2_2D_C_H_
#define AMREX_EB2_2D_C_H_

namespace amrex { namespace EB2 {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
amrex_eb2_build_types (Box const& tbx, Box const& bxg2,
                       Array4<Real const> const& s,
                       Array4<EBCellFlag> const& cell,
                       Array4<Type_t> const& fx,
                       Array4<Type_t> const& fy)
{
    auto lo = amrex::max_lbound(tbx, bxg2);
    auto hi = amrex::min_ubound(tbx, bxg2);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (    s(i,j  ,k) < 0.0 and s(i+1,j  ,k) < 0.0
            and s(i,j+1,k) < 0.0 and s(i+1,j+1,k) < 0.0)
        {
            cell(i,j,k).setRegular();
        }
        else if (s(i,j  ,k) >= 0.0 and s(i+1,j  ,k) >= 0.0
            and  s(i,j+1,k) >= 0.0 and s(i+1,j+1,k) >= 0.0)
        {
            cell(i,j,k).setCovered();
        }
        else
        {
            cell(i,j,k).setSingleValued();
        }
    });

    // x-face
    const Box& xbx = amrex::surroundingNodes(bxg2,0);
    lo = amrex::max_lbound(tbx, xbx);
    hi = amrex::min_ubound(tbx, xbx);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (s(i,j,k) < 0.0 and s(i,j+1,k) < 0.0) {
            fx(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i,j+1,k) >= 0.0) {
            fx(i,j,k) = Type::covered;
        } else {
            fx(i,j,k) = Type::irregular;
        }
    });

    // y-face
    const Box& ybx = amrex::surroundingNodes(bxg2,1);
    lo = amrex::max_lbound(tbx, ybx);
    hi = amrex::min_ubound(tbx, ybx);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (s(i,j,k) < 0.0 and s(i+1,j,k) < 0.0) {
            fy(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i+1,j,k) >= 0.0) {
            fy(i,j,k) = Type::covered;
        } else {
            fy(i,j,k) = Type::irregular;
        }
    });
}

}}

#endif
