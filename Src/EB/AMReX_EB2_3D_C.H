#ifndef AMREX_EB2_3D_C_H_
#define AMREX_EB2_3D_C_H_

namespace amrex { namespace EB2 {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
amrex_eb2_build_types (Box const& tbx, Box const& bxg2,
                       Array4<Real const> const& s,
                       Array4<EBCellFlag> const& cell,
                       Array4<Type_t> const& fx,
                       Array4<Type_t> const& fy,
                       Array4<Type_t> const& fz,
                       Array4<Type_t> const& ex,
                       Array4<Type_t> const& ey,
                       Array4<Type_t> const& ez)
{
    auto lo = amrex::max_lbound(tbx, bxg2);
    auto hi = amrex::min_ubound(tbx, bxg2);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (    s(i,j  ,k  ) < 0.0 and s(i+1,j  ,k  ) < 0.0
            and s(i,j+1,k  ) < 0.0 and s(i+1,j+1,k  ) < 0.0
            and s(i,j  ,k+1) < 0.0 and s(i+1,j  ,k+1) < 0.0
            and s(i,j+1,k+1) < 0.0 and s(i+1,j+1,k+1) < 0.0)
        {
            cell(i,j,k).setRegular();
        }
        else if (s(i,j  ,k  ) >= 0.0 and s(i+1,j  ,k  ) >= 0.0
            and  s(i,j+1,k  ) >= 0.0 and s(i+1,j+1,k  ) >= 0.0
            and  s(i,j  ,k+1) >= 0.0 and s(i+1,j  ,k+1) >= 0.0
            and  s(i,j+1,k+1) >= 0.0 and s(i+1,j+1,k+1) >= 0.0)
        {
            cell(i,j,k).setCovered();
        }
        else
        {
            cell(i,j,k).setSingleValued();
        }
    });

    // x-face
    Box b = amrex::surroundingNodes(bxg2,0);
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (    s(i,j,k  ) < 0.0 and s(i,j+1,k  ) < 0.0
            and s(i,j,k+1) < 0.0 and s(i,j+1,k+1) < 0.0 )
        {
            fx(i,j,k) = Type::regular;
        }
        else if (s(i,j,k  ) >= 0.0 and s(i,j+1,k  ) >= 0.0 
            and  s(i,j,k+1) >= 0.0 and s(i,j+1,k+1) >= 0.0 )
        {
            fx(i,j,k) = Type::covered;
        }
        else
        {
            fx(i,j,k) = Type::irregular;
        }
    });

    // y-face
    b = amrex::surroundingNodes(bxg2,1);
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (    s(i,j,k  ) < 0.0 and s(i+1,j,k  ) < 0.0
            and s(i,j,k+1) < 0.0 and s(i+1,j,k+1) < 0.0 )
        {
            fy(i,j,k) = Type::regular;
        }
        else if (s(i,j,k  ) >= 0.0 and s(i+1,j,k  ) >= 0.0
            and  s(i,j,k+1) >= 0.0 and s(i+1,j,k+1) >= 0.0 )
        {
            fy(i,j,k) = Type::covered;
        }
        else
        {
            fy(i,j,k) = Type::irregular;
        }
    });

    // z-face
    b = amrex::surroundingNodes(bxg2,2);
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (    s(i,j  ,k) < 0.0 and s(i+1,j  ,k) < 0.0
            and s(i,j+1,k) < 0.0 and s(i+1,j+1,k) < 0.0)
        {
            fz(i,j,k) = Type::regular;
        }
        else if (s(i,j  ,k) >= 0.0 and s(i+1,j  ,k) >= 0.0
            and  s(i,j+1,k) >= 0.0 and s(i+1,j+1,k) >= 0.0)
        {
            fz(i,j,k) = Type::covered;
        }
        else
        {
            fz(i,j,k) = Type::irregular;
        }
    });

    // x-edge
    b = amrex::convert(bxg2,IntVect(0,1,1));
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (s(i,j,k) < 0.0 and s(i+1,j,k) < 0.0) {
            ex(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i+1,j,k) >= 0.0) {
            ex(i,j,k) = Type::covered;
        } else {
            ex(i,j,k) = Type::irregular;
        }
    });

    // y-edge
    b = amrex::convert(bxg2,IntVect(1,0,1));
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (s(i,j,k) < 0.0 and s(i,j+1,k) < 0.0) {
            ey(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i,j+1,k) >= 0.0) {
            ey(i,j,k) = Type::covered;
        } else {
            ey(i,j,k) = Type::irregular;
        }
    });

    // z-edge
    b = amrex::convert(bxg2,IntVect(1,1,0));
    lo = amrex::max_lbound(tbx, b);
    hi = amrex::min_ubound(tbx, b);
    amrex::Loop(lo, hi,
    [=] (int i, int j, int k) noexcept
    {
        if (s(i,j,k) < 0.0 and s(i,j,k+1) < 0.0) {
            ez(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i,j,k+1) >= 0.0) {
            ez(i,j,k) = Type::covered;
        } else {
            ez(i,j,k) = Type::irregular;
        }
    });
}

}}

#endif
