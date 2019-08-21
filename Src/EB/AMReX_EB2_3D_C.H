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

namespace {
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int num_cuts (Real a, Real b) noexcept {
        return (a >= 0.0 and b < 0.0) or (b >= 0.0 and a < 0.0);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int check_mvmc (int i, int j, int k, Array4<Real const> const& fine)
{
    int ierr = 0;

    i *= 2;
    j *= 2;
    k *= 2;

    // x-edges
    int nx00 = num_cuts(fine(i,j,k),fine(i+1,j,k)) + num_cuts(fine(i+1,j,k),fine(i+2,j,k));
    int nx10 = num_cuts(fine(i,j+1,k),fine(i+1,j+1,k)) + num_cuts(fine(i+1,j+1,k),fine(i+2,j+1,k));
    int nx01 = num_cuts(fine(i,j,k+1),fine(i+1,j,k+1)) + num_cuts(fine(i+1,j,k+1),fine(i+2,j,k+1));
    int nx11 = num_cuts(fine(i,j+1,k+1),fine(i+1,j+1,k+1)) + num_cuts(fine(i+1,j+1,k+1),fine(i+2,j+1,k+1));

    // y-edges
    int ny00 = num_cuts(fine(i,j,k),fine(i,j+1,k)) + num_cuts(fine(i,j+1,k),fine(i,j+2,k));
    int ny10 = num_cuts(fine(i+1,j,k),fine(i+1,j+1,k)) + num_cuts(fine(i+1,j+1,k),fine(i+1,j+2,k));
    int ny01 = num_cuts(fine(i,j,k+1),fine(i,j+1,k+1)) + num_cuts(fine(i,j+1,k+1),fine(i,j+2,k+1));
    int ny11 = num_cuts(fine(i+1,j,k+1),fine(i+1,j+1,k+1)) + num_cuts(fine(i+1,j+1,k+1),fine(i+1,j+2,k+1));

    // z-edges
    int nz00 = num_cuts(fine(i,j,k),fine(i,j,k+1)) + num_cuts(fine(i,j,k+1),fine(i,j,k+2));
    int nz10 = num_cuts(fine(i+1,j,k),fine(i+1,j,k+1)) + num_cuts(fine(i+1,j,k+1),fine(i+1,j,k+2));
    int nz01 = num_cuts(fine(i,j+1,k),fine(i,j+1,k+1)) + num_cuts(fine(i,j+1,k+1),fine(i,j+1,k+2));
    int nz11 = num_cuts(fine(i+1,j+1,k),fine(i+1,j+1,k+1)) + num_cuts(fine(i+1,j+1,k+1),fine(i+1,j+1,k+2));

    // x-faces
    int nxm;
    int n = ny00 + ny01 + nz00 + nz01;
    if (n == 0) {
        nxm = 0;
    } else if (n == 2) {
        nxm = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on x-face");
    }

    int nxp;
    n = ny10 + ny11 + nz10 + nz11;
    if (n == 0) {
        nxp = 0;
    } else if (n == 2) {
        nxp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on x-face");
    }

    // y-faces
    int nym;
    n = nx00 + nx01 + nz00 + nz10;
    if (n == 0) {
        nym = 0;
    } else if (n == 2) {
        nym = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on y-face");
    }

    int nyp;
    n = nx10 + nx11 + nz01 + nz11;
    if (n == 0) {
        nyp = 0;
    } else if (n == 2) {
        nyp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on y-face");
    }

    // z-faces
    int nzm;
    n = nx00 + nx10 + ny00 + ny10;
    if (n == 0) {
        nzm = 0;
    } else if (n == 2) {
        nzm = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on z-face");
    }

    int nzp;
    n = nx01 + nx11 + ny01 + ny11;
    if (n == 0) {
        nzp = 0;
    } else if (n == 2) {
        nzp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on z-face");
    }

    if (nxm == 1 and nym == 1 and nzm == 1 and nxp == 1 and nyp == 1 and nzp == 1) {
        n = (fine(i  ,j  ,k  ) < 0.0) + (fine(i+2,j  ,k  ) < 0.0) +
            (fine(i  ,j+2,k  ) < 0.0) + (fine(i+2,j+2,k  ) < 0.0) +
            (fine(i  ,j  ,k+2) < 0.0) + (fine(i+2,j  ,k+2) < 0.0) +
            (fine(i  ,j+2,k+2) < 0.0) + (fine(i+2,j+2,k+2) < 0.0);
        if (n == 2 or n == 6) {
            ierr = 1;
        } else if (n != 4) {
            ierr = 1;
            amrex::Abort("amrex::check_mvmc: how did this happen? nopen != 4");
        }
    }

    return ierr;
}

}}

#endif
