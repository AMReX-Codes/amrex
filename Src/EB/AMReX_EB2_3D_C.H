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
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int num_cuts (Real a, Real b) noexcept {
        return (a >= 0.0 and b < 0.0) or (b >= 0.0 and a < 0.0);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int check_mvmc (int i, int j, int k, Array4<Real const> const& fine)
{
    int ierr = 0;

    i *= 2;
    j *= 2;
    k *= 2;

    // x-edges
    int nx00 = num_cuts(fine(i,j,k),fine(i+1,j,k)) + num_cuts(fine(i+1,j,k),fine(i+2,j,k));
    int nx10 = num_cuts(fine(i,j+2,k),fine(i+1,j+2,k)) + num_cuts(fine(i+1,j+2,k),fine(i+2,j+2,k));
    int nx01 = num_cuts(fine(i,j,k+2),fine(i+1,j,k+2)) + num_cuts(fine(i+1,j,k+2),fine(i+2,j,k+2));
    int nx11 = num_cuts(fine(i,j+2,k+2),fine(i+1,j+2,k+2)) + num_cuts(fine(i+1,j+2,k+2),fine(i+2,j+2,k+2));

    // y-edges
    int ny00 = num_cuts(fine(i,j,k),fine(i,j+1,k)) + num_cuts(fine(i,j+1,k),fine(i,j+2,k));
    int ny10 = num_cuts(fine(i+2,j,k),fine(i+2,j+1,k)) + num_cuts(fine(i+2,j+1,k),fine(i+2,j+2,k));
    int ny01 = num_cuts(fine(i,j,k+2),fine(i,j+1,k+2)) + num_cuts(fine(i,j+1,k+2),fine(i,j+2,k+2));
    int ny11 = num_cuts(fine(i+2,j,k+2),fine(i+2,j+1,k+2)) + num_cuts(fine(i+2,j+1,k+2),fine(i+2,j+2,k+2));

    // z-edges
    int nz00 = num_cuts(fine(i,j,k),fine(i,j,k+1)) + num_cuts(fine(i,j,k+1),fine(i,j,k+2));
    int nz10 = num_cuts(fine(i+2,j,k),fine(i+2,j,k+1)) + num_cuts(fine(i+2,j,k+1),fine(i+2,j,k+2));
    int nz01 = num_cuts(fine(i,j+2,k),fine(i,j+2,k+1)) + num_cuts(fine(i,j+2,k+1),fine(i,j+2,k+2));
    int nz11 = num_cuts(fine(i+2,j+2,k),fine(i+2,j+2,k+1)) + num_cuts(fine(i+2,j+2,k+1),fine(i+2,j+2,k+2));

    // x-faces
    int nxm = -1;
    int n = ny00 + ny01 + nz00 + nz01;
    if (n == 0) {
        nxm = 0;
    } else if (n == 2) {
        nxm = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on xlo-face");
    }

    int nxp = -1;
    n = ny10 + ny11 + nz10 + nz11;
    if (n == 0) {
        nxp = 0;
    } else if (n == 2) {
        nxp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on xhi-face");
    }

    // y-faces
    int nym = -1;
    n = nx00 + nx01 + nz00 + nz10;
    if (n == 0) {
        nym = 0;
    } else if (n == 2) {
        nym = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on ylo-face");
    }

    int nyp = -1;
    n = nx10 + nx11 + nz01 + nz11;
    if (n == 0) {
        nyp = 0;
    } else if (n == 2) {
        nyp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on yhi-face");
    }

    // z-faces
    int nzm = -1;
    n = nx00 + nx10 + ny00 + ny10;
    if (n == 0) {
        nzm = 0;
    } else if (n == 2) {
        nzm = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on zlo-face");
    }

    int nzp = -1;
    n = nx01 + nx11 + ny01 + ny11;
    if (n == 0) {
        nzp = 0;
    } else if (n == 2) {
        nzp = 1;
    } else if (n == 4) {
        ierr = 1;
    } else {
        ierr = 1;
        amrex::Abort("amrex::check_mvmc: how did this happen? wrong number of cuts on zhi-face");
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen_from_fine (int i, int j, int k, Box const& bx, Box const& gbx,
                       Box const& xbx, Box const& ybx, Box const& zbx,
                       Box const& xgbx, Box const& ygbx, Box const& zgbx,
                       Array4<Real> const& cvol, Array4<Real> const& ccent,
                       Array4<Real> const& cba, Array4<Real> const& cbc,
                       Array4<Real> const& cbn, Array4<Real> const& capx,
                       Array4<Real> const& capy, Array4<Real> const& capz,
                       Array4<Real> const& cfcx, Array4<Real> const& cfcy,
                       Array4<Real> const& cfcz, Array4<EBCellFlag> const& cflag,
                       Array4<Real const> const& fvol, Array4<Real const> const& fcent,
                       Array4<Real const> const& fba, Array4<Real const> const& fbc,
                       Array4<Real const> const& fbn, Array4<Real const> const& fapx,
                       Array4<Real const> const& fapy, Array4<Real const> const& fapz,
                       Array4<Real const> const& ffcx, Array4<Real const> const& ffcy,
                       Array4<Real const> const& ffcz, Array4<EBCellFlag const> const& fflag)
{
    int ierr = 0;
    IntVect iv(i,j,k);
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;

    if (bx.contains(iv))
    {
        if (fflag(ii,jj  ,kk  ).isRegular() && fflag(ii+1,jj  ,kk  ).isRegular() &&
            fflag(ii,jj+1,kk  ).isRegular() && fflag(ii+1,jj+1,kk  ).isRegular() &&
            fflag(ii,jj  ,kk+1).isRegular() && fflag(ii+1,jj  ,kk+1).isRegular() &&
            fflag(ii,jj+1,kk+1).isRegular() && fflag(ii+1,jj+1,kk+1).isRegular())
        {
            cflag(i,j,k).setRegular();
            cvol(i,j,k) = 1.0;
            ccent(i,j,k,0) = 0.0;
            ccent(i,j,k,1) = 0.0;
            ccent(i,j,k,2) = 0.0;
            cba(i,j,k) = 0.0;
            cbc(i,j,k,0) = -1.0;
            cbc(i,j,k,1) = -1.0;
            cbc(i,j,k,2) = -1.0;
            cbn(i,j,k,0) = 0.0;
            cbn(i,j,k,1) = 0.0;
            cbn(i,j,k,2) = 0.0;
        }
        else if (fflag(ii,jj  ,kk  ).isCovered() && fflag(ii+1,jj  ,kk  ).isCovered() &&
                 fflag(ii,jj+1,kk  ).isCovered() && fflag(ii+1,jj+1,kk  ).isCovered() &&
                 fflag(ii,jj  ,kk+1).isCovered() && fflag(ii+1,jj  ,kk+1).isCovered() &&
                 fflag(ii,jj+1,kk+1).isCovered() && fflag(ii+1,jj+1,kk+1).isCovered())
        {
            cflag(i,j,k).setCovered();
            cvol(i,j,k) = 0.0;
            ccent(i,j,k,0) = 0.0;
            ccent(i,j,k,1) = 0.0;
            ccent(i,j,k,2) = 0.0;
            cba(i,j,k) = 0.0;
            cbc(i,j,k,0) = -1.0;
            cbc(i,j,k,1) = -1.0;
            cbc(i,j,k,2) = -1.0;
            cbn(i,j,k,0) = 0.0;
            cbn(i,j,k,1) = 0.0;
            cbn(i,j,k,2) = 0.0;
        }
        else
        {
            cflag(i,j,k).setSingleValued();

            cvol(i,j,k) = 0.125*(fvol(ii  ,jj  ,kk  ) + fvol(ii+1,jj  ,kk  ) +
                                 fvol(ii  ,jj+1,kk  ) + fvol(ii+1,jj+1,kk  ) +
                                 fvol(ii  ,jj  ,kk+1) + fvol(ii+1,jj  ,kk+1) +
                                 fvol(ii  ,jj+1,kk+1) + fvol(ii+1,jj+1,kk+1));
            Real cvolinv = 1.0 / cvol(i,j,k);

            ccent(i,j,k,0) = 0.125 * cvolinv *
                (fvol(ii  ,jj  ,kk  )*(0.5*fcent(ii  ,jj  ,kk  ,0)-0.25) +
                 fvol(ii+1,jj  ,kk  )*(0.5*fcent(ii+1,jj  ,kk  ,0)+0.25) +
                 fvol(ii  ,jj+1,kk  )*(0.5*fcent(ii  ,jj+1,kk  ,0)-0.25) +
                 fvol(ii+1,jj+1,kk  )*(0.5*fcent(ii+1,jj+1,kk  ,0)+0.25) +
                 fvol(ii  ,jj  ,kk+1)*(0.5*fcent(ii  ,jj  ,kk+1,0)-0.25) +
                 fvol(ii+1,jj  ,kk+1)*(0.5*fcent(ii+1,jj  ,kk+1,0)+0.25) +
                 fvol(ii  ,jj+1,kk+1)*(0.5*fcent(ii  ,jj+1,kk+1,0)-0.25) +
                 fvol(ii+1,jj+1,kk+1)*(0.5*fcent(ii+1,jj+1,kk+1,0)+0.25));
            ccent(i,j,k,1) = 0.125 * cvolinv *
                (fvol(ii  ,jj  ,kk  )*(0.5*fcent(ii  ,jj  ,kk  ,1)-0.25) +
                 fvol(ii+1,jj  ,kk  )*(0.5*fcent(ii+1,jj  ,kk  ,1)-0.25) +
                 fvol(ii  ,jj+1,kk  )*(0.5*fcent(ii  ,jj+1,kk  ,1)+0.25) +
                 fvol(ii+1,jj+1,kk  )*(0.5*fcent(ii+1,jj+1,kk  ,1)+0.25) +
                 fvol(ii  ,jj  ,kk+1)*(0.5*fcent(ii  ,jj  ,kk+1,1)-0.25) +
                 fvol(ii+1,jj  ,kk+1)*(0.5*fcent(ii+1,jj  ,kk+1,1)-0.25) +
                 fvol(ii  ,jj+1,kk+1)*(0.5*fcent(ii  ,jj+1,kk+1,1)+0.25) +
                 fvol(ii+1,jj+1,kk+1)*(0.5*fcent(ii+1,jj+1,kk+1,1)+0.25));
            ccent(i,j,k,2) = 0.125 * cvolinv *
                (fvol(ii  ,jj  ,kk  )*(0.5*fcent(ii  ,jj  ,kk  ,2)-0.25) +
                 fvol(ii+1,jj  ,kk  )*(0.5*fcent(ii+1,jj  ,kk  ,2)-0.25) +
                 fvol(ii  ,jj+1,kk  )*(0.5*fcent(ii  ,jj+1,kk  ,2)-0.25) +
                 fvol(ii+1,jj+1,kk  )*(0.5*fcent(ii+1,jj+1,kk  ,2)-0.25) +
                 fvol(ii  ,jj  ,kk+1)*(0.5*fcent(ii  ,jj  ,kk+1,2)+0.25) +
                 fvol(ii+1,jj  ,kk+1)*(0.5*fcent(ii+1,jj  ,kk+1,2)+0.25) +
                 fvol(ii  ,jj+1,kk+1)*(0.5*fcent(ii  ,jj+1,kk+1,2)+0.25) +
                 fvol(ii+1,jj+1,kk+1)*(0.5*fcent(ii+1,jj+1,kk+1,2)+0.25));

            cba(i,j,k) = 0.25*(fba(ii  ,jj  ,kk  ) + fba(ii+1,jj  ,kk  ) +
                               fba(ii  ,jj+1,kk  ) + fba(ii+1,jj+1,kk  ) +
                               fba(ii  ,jj  ,kk+1) + fba(ii+1,jj  ,kk+1) +
                               fba(ii  ,jj+1,kk+1) + fba(ii+1,jj+1,kk+1));
            Real cbainv = 1.0 / cba(i,j,k);

            cbc(i,j,k,0) = 0.25 * cbainv *
                  ( fba(ii  ,jj  ,kk  )*(0.5*fbc(ii  ,jj  ,kk  ,0)-0.25)
                  + fba(ii+1,jj  ,kk  )*(0.5*fbc(ii+1,jj  ,kk  ,0)+0.25)
                  + fba(ii  ,jj+1,kk  )*(0.5*fbc(ii  ,jj+1,kk  ,0)-0.25)
                  + fba(ii+1,jj+1,kk  )*(0.5*fbc(ii+1,jj+1,kk  ,0)+0.25)
                  + fba(ii  ,jj  ,kk+1)*(0.5*fbc(ii  ,jj  ,kk+1,0)-0.25)
                  + fba(ii+1,jj  ,kk+1)*(0.5*fbc(ii+1,jj  ,kk+1,0)+0.25)
                  + fba(ii  ,jj+1,kk+1)*(0.5*fbc(ii  ,jj+1,kk+1,0)-0.25)
                  + fba(ii+1,jj+1,kk+1)*(0.5*fbc(ii+1,jj+1,kk+1,0)+0.25) );
            cbc(i,j,k,1) = 0.25 * cbainv *
                  ( fba(ii  ,jj  ,kk  )*(0.5*fbc(ii  ,jj  ,kk  ,1)-0.25)
                  + fba(ii+1,jj  ,kk  )*(0.5*fbc(ii+1,jj  ,kk  ,1)-0.25)
                  + fba(ii  ,jj+1,kk  )*(0.5*fbc(ii  ,jj+1,kk  ,1)+0.25)
                  + fba(ii+1,jj+1,kk  )*(0.5*fbc(ii+1,jj+1,kk  ,1)+0.25)
                  + fba(ii  ,jj  ,kk+1)*(0.5*fbc(ii  ,jj  ,kk+1,1)-0.25)
                  + fba(ii+1,jj  ,kk+1)*(0.5*fbc(ii+1,jj  ,kk+1,1)-0.25)
                  + fba(ii  ,jj+1,kk+1)*(0.5*fbc(ii  ,jj+1,kk+1,1)+0.25)
                  + fba(ii+1,jj+1,kk+1)*(0.5*fbc(ii+1,jj+1,kk+1,1)+0.25) );
            cbc(i,j,k,2) = 0.25 * cbainv *
                  ( fba(ii  ,jj  ,kk  )*(0.5*fbc(ii  ,jj  ,kk  ,2)-0.25)
                  + fba(ii+1,jj  ,kk  )*(0.5*fbc(ii+1,jj  ,kk  ,2)-0.25)
                  + fba(ii  ,jj+1,kk  )*(0.5*fbc(ii  ,jj+1,kk  ,2)-0.25)
                  + fba(ii+1,jj+1,kk  )*(0.5*fbc(ii+1,jj+1,kk  ,2)-0.25)
                  + fba(ii  ,jj  ,kk+1)*(0.5*fbc(ii  ,jj  ,kk+1,2)+0.25)
                  + fba(ii+1,jj  ,kk+1)*(0.5*fbc(ii+1,jj  ,kk+1,2)+0.25)
                  + fba(ii  ,jj+1,kk+1)*(0.5*fbc(ii  ,jj+1,kk+1,2)+0.25)
                  + fba(ii+1,jj+1,kk+1)*(0.5*fbc(ii+1,jj+1,kk+1,2)+0.25) );

            Real nx = fbn(ii  ,jj  ,kk  ,0)*fba(ii  ,jj  ,kk  )
                +     fbn(ii+1,jj  ,kk  ,0)*fba(ii+1,jj  ,kk  )
                +     fbn(ii  ,jj+1,kk  ,0)*fba(ii  ,jj+1,kk  )
                +     fbn(ii+1,jj+1,kk  ,0)*fba(ii+1,jj+1,kk  )
                +     fbn(ii  ,jj  ,kk+1,0)*fba(ii  ,jj  ,kk+1)
                +     fbn(ii+1,jj  ,kk+1,0)*fba(ii+1,jj  ,kk+1)
                +     fbn(ii  ,jj+1,kk+1,0)*fba(ii  ,jj+1,kk+1)
                +     fbn(ii+1,jj+1,kk+1,0)*fba(ii+1,jj+1,kk+1);
            Real ny = fbn(ii  ,jj  ,kk  ,1)*fba(ii  ,jj  ,kk  )
                +     fbn(ii+1,jj  ,kk  ,1)*fba(ii+1,jj  ,kk  )
                +     fbn(ii  ,jj+1,kk  ,1)*fba(ii  ,jj+1,kk  )
                +     fbn(ii+1,jj+1,kk  ,1)*fba(ii+1,jj+1,kk  )
                +     fbn(ii  ,jj  ,kk+1,1)*fba(ii  ,jj  ,kk+1)
                +     fbn(ii+1,jj  ,kk+1,1)*fba(ii+1,jj  ,kk+1)
                +     fbn(ii  ,jj+1,kk+1,1)*fba(ii  ,jj+1,kk+1)
                +     fbn(ii+1,jj+1,kk+1,1)*fba(ii+1,jj+1,kk+1);
            Real nz = fbn(ii  ,jj  ,kk  ,2)*fba(ii  ,jj  ,kk  )
                +     fbn(ii+1,jj  ,kk  ,2)*fba(ii+1,jj  ,kk  )
                +     fbn(ii  ,jj+1,kk  ,2)*fba(ii  ,jj+1,kk  )
                +     fbn(ii+1,jj+1,kk  ,2)*fba(ii+1,jj+1,kk  )
                +     fbn(ii  ,jj  ,kk+1,2)*fba(ii  ,jj  ,kk+1)
                +     fbn(ii+1,jj  ,kk+1,2)*fba(ii+1,jj  ,kk+1)
                +     fbn(ii  ,jj+1,kk+1,2)*fba(ii  ,jj+1,kk+1)
                +     fbn(ii+1,jj+1,kk+1,2)*fba(ii+1,jj+1,kk+1);
            Real nfac = 1.0 / std::sqrt(nx*nx+ny*ny+nz*nz+1.e-50);
            cbn(i,j,k,0) = nx*nfac;
            cbn(i,j,k,1) = ny*nfac;
            cbn(i,j,k,2) = nz*nfac;
            ierr = (nx == 0.0 and ny == 0.0 and nz == 0.0)
                // we must check the enclosing surface to make sure the coarse cell does not
                // fully contains the geometry object.
                or ( fapx(ii,jj  ,kk  )==1.0 and fapx(ii+2,jj  ,kk  )==1.0
                 and fapx(ii,jj+1,kk  )==1.0 and fapx(ii+2,jj+1,kk  )==1.0
                 and fapx(ii,jj  ,kk+1)==1.0 and fapx(ii+2,jj  ,kk+1)==1.0
                 and fapx(ii,jj+1,kk+1)==1.0 and fapx(ii+2,jj+1,kk+1)==1.0
                 and fapy(ii,jj  ,kk  )==1.0 and fapy(ii+1,jj  ,kk  )==1.0
                 and fapy(ii,jj+2,kk  )==1.0 and fapy(ii+1,jj+2,kk  )==1.0
                 and fapy(ii,jj  ,kk+1)==1.0 and fapy(ii+1,jj  ,kk+1)==1.0
                 and fapy(ii,jj+2,kk+1)==1.0 and fapy(ii+1,jj+2,kk+1)==1.0
                 and fapz(ii,jj  ,kk  )==1.0 and fapz(ii+1,jj  ,kk  )==1.0
                 and fapz(ii,jj+1,kk  )==1.0 and fapz(ii+1,jj+1,kk  )==1.0
                 and fapz(ii,jj  ,kk+2)==1.0 and fapz(ii+1,jj  ,kk+2)==1.0
                 and fapz(ii,jj+1,kk+2)==1.0 and fapz(ii+1,jj+1,kk+2)==1.0);

        }
    }
    else if (gbx.contains(iv))
    {
        cvol(i,j,k) = 1.0;
        ccent(i,j,k,0) = 0.0;
        ccent(i,j,k,1) = 0.0;
        ccent(i,j,k,2) = 0.0;
        cba(i,j,k) = 0.0;
        cbc(i,j,k,0) = -1.0;
        cbc(i,j,k,1) = -1.0;
        cbc(i,j,k,2) = -1.0;
        cbn(i,j,k,0) = 0.0;
        cbn(i,j,k,1) = 0.0;
        cbn(i,j,k,2) = 0.0;
    }

    if (xbx.contains(iv))
    {
        capx(i,j,k) = 0.25*(fapx(ii,jj  ,kk  ) + fapx(ii,jj+1,kk  ) +
                            fapx(ii,jj  ,kk+1) + fapx(ii,jj+1,kk+1));
        if (capx(i,j,k) != 0.0) {
            Real apinv = 1.0 / capx(i,j,k);
            cfcx(i,j,k,0) = 0.25 * apinv *
                (fapx(ii,jj  ,kk  )*(0.5*ffcx(ii,jj  ,kk  ,0)-0.25) +
                 fapx(ii,jj+1,kk  )*(0.5*ffcx(ii,jj+1,kk  ,0)+0.25) +
                 fapx(ii,jj  ,kk+1)*(0.5*ffcx(ii,jj  ,kk+1,0)-0.25) +
                 fapx(ii,jj+1,kk+1)*(0.5*ffcx(ii,jj+1,kk+1,0)+0.25) );
            cfcx(i,j,k,1) = 0.25 * apinv *
                (fapx(ii,jj  ,kk  )*(0.5*ffcx(ii,jj  ,kk  ,1)-0.25) +
                 fapx(ii,jj+1,kk  )*(0.5*ffcx(ii,jj+1,kk  ,1)-0.25) +
                 fapx(ii,jj  ,kk+1)*(0.5*ffcx(ii,jj  ,kk+1,1)+0.25) +
                 fapx(ii,jj+1,kk+1)*(0.5*ffcx(ii,jj+1,kk+1,1)+0.25) );
        } else {
            cfcx(i,j,k,0) = 0.0;
            cfcx(i,j,k,1) = 0.0;
        }
    }
    else if (xgbx.contains(iv))
    {
        capx(i,j,k) = 1.0;
        cfcx(i,j,k,0) = 0.0;
        cfcx(i,j,k,1) = 0.0;
    }

    if (ybx.contains(iv))
    {
        capy(i,j,k) = 0.25*(fapy(ii  ,jj,kk  ) + fapy(ii+1,jj,kk  ) +
                            fapy(ii  ,jj,kk+1) + fapy(ii+1,jj,kk+1));
        if (capy(i,j,k) != 0.0) {
            Real apinv = 1.0 / capy(i,j,k);
            cfcy(i,j,k,0) = 0.25 * apinv *
                (fapy(ii  ,jj,kk  )*(0.5*ffcy(ii  ,jj,kk  ,0)-0.25) +
                 fapy(ii+1,jj,kk  )*(0.5*ffcy(ii+1,jj,kk  ,0)+0.25) +
                 fapy(ii  ,jj,kk+1)*(0.5*ffcy(ii  ,jj,kk+1,0)-0.25) +
                 fapy(ii+1,jj,kk+1)*(0.5*ffcy(ii+1,jj,kk+1,0)+0.25) );
            cfcy(i,j,k,1) = 0.25 * apinv *
                (fapy(ii  ,jj,kk  )*(0.5*ffcy(ii  ,jj,kk  ,1)-0.25) +
                 fapy(ii+1,jj,kk  )*(0.5*ffcy(ii+1,jj,kk  ,1)-0.25) +
                 fapy(ii  ,jj,kk+1)*(0.5*ffcy(ii  ,jj,kk+1,1)+0.25) +
                 fapy(ii+1,jj,kk+1)*(0.5*ffcy(ii+1,jj,kk+1,1)+0.25) );
        } else {
            cfcy(i,j,k,0) = 0.0;
            cfcy(i,j,k,1) = 0.0;
        }
    }
    else if (ygbx.contains(iv))
    {
        capy(i,j,k) = 1.0;
        cfcy(i,j,k,0) = 0.0;
        cfcy(i,j,k,1) = 0.0;
    }

    if (zbx.contains(iv))
    {
        capz(i,j,k) = 0.25 * (fapz(ii  ,jj  ,kk) + fapz(ii+1,jj  ,kk) +
                              fapz(ii  ,jj+1,kk) + fapz(ii+1,jj+1,kk));
        if (capz(i,j,k) != 0.0) {
            Real apinv = 1.0 / capz(i,j,k);
            cfcz(i,j,k,0) = 0.25 * apinv *
                (fapz(ii  ,jj  ,kk)*(0.5*ffcz(ii  ,jj  ,kk,0)-0.25) +
                 fapz(ii+1,jj  ,kk)*(0.5*ffcz(ii+1,jj  ,kk,0)+0.25) +
                 fapz(ii  ,jj+1,kk)*(0.5*ffcz(ii  ,jj+1,kk,0)-0.25) +
                 fapz(ii+1,jj+1,kk)*(0.5*ffcz(ii+1,jj+1,kk,0)+0.25) );
            cfcz(i,j,k,1) = 0.25 * apinv *
                (fapz(ii  ,jj  ,kk)*(0.5*ffcz(ii  ,jj  ,kk,1)-0.25) +
                 fapz(ii+1,jj  ,kk)*(0.5*ffcz(ii+1,jj  ,kk,1)-0.25) +
                 fapz(ii  ,jj+1,kk)*(0.5*ffcz(ii  ,jj+1,kk,1)+0.25) +
                 fapz(ii+1,jj+1,kk)*(0.5*ffcz(ii+1,jj+1,kk,1)+0.25) );
        } else {
            cfcz(i,j,k,0) = 0.0;
            cfcz(i,j,k,1) = 0.0;
        }
    }
    else if (zgbx.contains(iv))
    {
        capz(i,j,k) = 1.0;
        cfcz(i,j,k,0) = 0.0;
        cfcz(i,j,k,1) = 0.0;
    }

    return ierr;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void build_cellflag_from_ap (int i, int j, int k, Array4<EBCellFlag> const& cflag,
                             Array4<Real const> const& apx, Array4<Real const> const& apy,
                             Array4<Real const> const& apz)
{
    auto flg = cflag(i,j,k);
    flg.setDisconnected();

    if (!flg.isCovered())
    {
        flg.setConnected(0,0,0);

        if (apx(i  ,j,k) != 0.0) flg.setConnected(-1,  0,  0);
        if (apx(i+1,j,k) != 0.0) flg.setConnected( 1,  0,  0);
        if (apy(i,j  ,k) != 0.0) flg.setConnected( 0, -1,  0);
        if (apy(i,j+1,k) != 0.0) flg.setConnected( 0,  1,  0);
        if (apz(i,j,k  ) != 0.0) flg.setConnected( 0,  0, -1);
        if (apz(i,j,k+1) != 0.0) flg.setConnected( 0,  0,  1);

        if ( (apx(i,j,k) != 0.0 and apy(i-1,j,k) != 0.0) or
             (apy(i,j,k) != 0.0 and apx(i,j-1,k) != 0.0) )
        {
            flg.setConnected(-1, -1, 0);
            if (apz(i-1,j-1,k  ) != 0.0) flg.setConnected(-1,-1,-1);
            if (apz(i-1,j-1,k+1) != 0.0) flg.setConnected(-1,-1, 1);
        }

        if ( (apx(i+1,j,k) != 0.0 and apy(i+1,j  ,k) != 0.0) or
             (apy(i  ,j,k) != 0.0 and apx(i+1,j-1,k) != 0.0) )
        {
            flg.setConnected(1, -1, 0);
            if (apz(i+1,j-1,k  ) != 0.0) flg.setConnected(1,-1,-1);
            if (apz(i+1,j-1,k+1) != 0.0) flg.setConnected(1,-1, 1);
        }

        if ( (apx(i,j  ,k) != 0.0 and apy(i-1,j+1,k) != 0.0) or
             (apy(i,j+1,k) != 0.0 and apx(i  ,j+1,k) != 0.0) )
        {
            flg.setConnected(-1, 1, 0);
            if (apz(i-1,j+1,k  ) != 0.0) flg.setConnected(-1, 1,-1);
            if (apz(i-1,j+1,k+1) != 0.0) flg.setConnected(-1, 1, 1);
        }

        if ( (apx(i+1,j  ,k) != 0.0 and apy(i+1,j+1,k) != 0.0) or
             (apy(i  ,j+1,k) != 0.0 and apx(i+1,j+1,k) != 0.0) )
        {
            flg.setConnected(1, 1, 0);
            if (apz(i+1,j+1,k  ) != 0.0) flg.setConnected(1, 1,-1);
            if (apz(i+1,j+1,k+1) != 0.0) flg.setConnected(1, 1, 1);
        }

        if ( (apx(i,j,k) != 0.0 and apz(i-1,j,k  ) != 0.0) or
             (apz(i,j,k) != 0.0 and apx(i  ,j,k-1) != 0.0) )
        {
            flg.setConnected(-1, 0, -1);
            if (apy(i-1,j  ,k-1) != 0.0) flg.setConnected(-1,-1,-1);
            if (apy(i-1,j+1,k-1) != 0.0) flg.setConnected(-1, 1,-1);
        }

        if ( (apx(i+1,j,k) != 0.0 and apz(i+1,j,k  ) != 0.0) or
             (apz(i  ,j,k) != 0.0 and apx(i+1,j,k-1) != 0.0) )
        {
            flg.setConnected(1, 0, -1);
            if (apy(i+1,j  ,k-1) != 0.0) flg.setConnected(1,-1,-1);
            if (apy(i+1,j+1,k-1) != 0.0) flg.setConnected(1, 1,-1);
        }

        if ( (apx(i,j,k  ) != 0.0 and apz(i-1,j,k+1) != 0.0) or
             (apz(i,j,k+1) != 0.0 and apx(i  ,j,k+1) != 0.0) )
        {
            flg.setConnected(-1, 0, 1);
            if (apy(i-1,j  ,k+1) != 0.0) flg.setConnected(-1,-1, 1);
            if (apy(i-1,j+1,k+1) != 0.0) flg.setConnected(-1, 1, 1);
        }

        if ( (apx(i+1,j,k  ) != 0.0 and apz(i+1,j,k+1) != 0.0) or
             (apz(i  ,j,k+1) != 0.0 and apx(i+1,j,k+1) != 0.0) )
        {
            flg.setConnected(1, 0, 1);
            if (apy(i+1,j  ,k+1) != 0.0) flg.setConnected(1,-1, 1);
            if (apy(i+1,j+1,k+1) != 0.0) flg.setConnected(1, 1, 1);
        }

        if ( (apy(i,j,k) != 0.0 and apz(i,j-1,k  ) != 0.0) or
             (apz(i,j,k) != 0.0 and apy(i,j  ,k-1) != 0.0) )
        {
            flg.setConnected(0, -1, -1);
            if (apx(i  ,j-1,k-1) != 0.0) flg.setConnected(-1,-1,-1);
            if (apx(i+1,j-1,k-1) != 0.0) flg.setConnected( 1,-1,-1);
        }

        if ( (apy(i,j+1,k) != 0.0 and apz(i,j+1,k  ) != 0.0) or
             (apz(i,j  ,k) != 0.0 and apy(i,j+1,k-1) != 0.0) )
        {
            flg.setConnected(0, 1, -1);
            if (apx(i  ,j+1,k-1) != 0.0) flg.setConnected(-1, 1,-1);
            if (apx(i+1,j+1,k-1) != 0.0) flg.setConnected( 1, 1,-1);
        }

        if ( (apy(i,j,k  ) != 0.0 and apz(i,j-1,k+1) != 0.0) or
             (apz(i,j,k+1) != 0.0 and apy(i,j  ,k+1) != 0.0) )
        {
            flg.setConnected(0, -1, 1);
            if (apx(i  ,j-1,k+1) != 0.0) flg.setConnected(-1,-1, 1);
            if (apx(i+1,j-1,k+1) != 0.0) flg.setConnected( 1,-1, 1);
        }

        if ( (apy(i,j+1,k  ) != 0.0 and apz(i,j+1,k+1) != 0.0) or
             (apz(i,j  ,k+1) != 0.0 and apy(i,j+1,k+1) != 0.0) )
        {
            flg.setConnected(0, 1, 1);
            if (apx(i  ,j+1,k+1) != 0.0) flg.setConnected(-1, 1, 1);
            if (apx(i+1,j+1,k+1) != 0.0) flg.setConnected( 1, 1, 1);
        }
    }

    cflag(i,j,k) = flg;
}

}}

#endif
