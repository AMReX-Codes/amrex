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

namespace {
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int num_cuts (Real a, Real b) noexcept {
        return (a >= 0.0 and b < 0.0) or (b >= 0.0 and a < 0.0);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int check_mvmc (int i, int j, int, Array4<Real const> const& fine)
{
    constexpr int k = 0;
    i *= 2;
    j *= 2;
    int ncuts = num_cuts(fine(i  ,j  ,k),fine(i+1,j  ,k))
        +       num_cuts(fine(i+1,j  ,k),fine(i+2,j  ,k))
        +       num_cuts(fine(i  ,j+2,k),fine(i+1,j+2,k))
        +       num_cuts(fine(i+1,j+2,k),fine(i+2,j+2,k))
        +       num_cuts(fine(i  ,j  ,k),fine(i  ,j+1,k))
        +       num_cuts(fine(i  ,j+1,k),fine(i  ,j+2,k))
        +       num_cuts(fine(i+2,j  ,k),fine(i+2,j+1,k))
        +       num_cuts(fine(i+2,j+1,k),fine(i+2,j+2,k));
    return (ncuts != 0 or ncuts != 2);
}

}}

#endif
