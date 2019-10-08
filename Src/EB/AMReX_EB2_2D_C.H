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
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int num_cuts (Real a, Real b) noexcept {
        return (a >= 0.0 and b < 0.0) or (b >= 0.0 and a < 0.0);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
    return (ncuts != 0 and ncuts != 2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen_from_fine (int i, int j, Box const& bx, Box const& gbx,
                       Box const& xbx, Box const& ybx, Box const& xgbx, Box const& ygbx,
                       Array4<Real> const& cvol, Array4<Real> const& ccent,
                       Array4<Real> const& cba, Array4<Real> const& cbc,
                       Array4<Real> const& cbn, Array4<Real> const& capx,
                       Array4<Real> const& capy, Array4<Real> const& cfcx,
                       Array4<Real> const& cfcy, Array4<EBCellFlag> const& cflag,
                       Array4<Real const> const& fvol, Array4<Real const> const& fcent,
                       Array4<Real const> const& fba, Array4<Real const> const& fbc,
                       Array4<Real const> const& fbn, Array4<Real const> const& fapx,
                       Array4<Real const> const& fapy, Array4<Real const> const& ffcx,
                       Array4<Real const> const& ffcy, Array4<EBCellFlag const> const& fflag)
{
    int ierr = 0;
    IntVect iv(i,j);
    constexpr int k = 0;
    int ii = i*2;
    int jj = j*2;

    if (bx.contains(iv))
    {
        if (fflag(ii,jj  ,k).isRegular() && fflag(ii+1,jj  ,k).isRegular() &&
            fflag(ii,jj+1,k).isRegular() && fflag(ii+1,jj+1,k).isRegular())
        {
            cflag(i,j,k).setRegular();
            cvol(i,j,k) = 1.0;
            ccent(i,j,k,0) = 0.0;
            ccent(i,j,k,1) = 0.0;
            cba(i,j,k) = 0.0;
            cbc(i,j,k,0) = -1.0;
            cbc(i,j,k,1) = -1.0;
            cbn(i,j,k,0) = 0.0;
            cbn(i,j,k,1) = 0.0;
        }
        else if (fflag(ii,jj  ,k).isCovered() && fflag(ii+1,jj  ,k).isCovered() &&
                 fflag(ii,jj+1,k).isCovered() && fflag(ii+1,jj+1,k).isCovered())
        {
            cflag(i,j,k).setCovered();
            cvol(i,j,k) = 0.0;
            ccent(i,j,k,0) = 0.0;
            ccent(i,j,k,1) = 0.0;
            cba(i,j,k) = 0.0;
            cbc(i,j,k,0) = -1.0;
            cbc(i,j,k,1) = -1.0;
            cbn(i,j,k,0) = 0.0;
            cbn(i,j,k,1) = 0.0;
        }
        else
        {
            cflag(i,j,k).setSingleValued();

            cvol(i,j,k) = 0.25*(fvol(ii,jj  ,k) + fvol(ii+1,jj  ,k) +
                                fvol(ii,jj+1,k) + fvol(ii+1,jj+1,k));
            Real cvolinv = 1.0/cvol(i,j,k);

            ccent(i,j,k,0) = 0.25 * cvolinv *
                (fvol(ii  ,jj  ,k)*(0.5*fcent(ii  ,jj  ,k,0)-0.25) +
                 fvol(ii+1,jj  ,k)*(0.5*fcent(ii+1,jj  ,k,0)+0.25) +
                 fvol(ii  ,jj+1,k)*(0.5*fcent(ii  ,jj+1,k,0)-0.25) +
                 fvol(ii+1,jj+1,k)*(0.5*fcent(ii+1,jj+1,k,0)+0.25));
            ccent(i,j,k,1) = 0.25 * cvolinv *
                (fvol(ii  ,jj  ,k)*(0.5*fcent(ii  ,jj  ,k,1)-0.25) +
                 fvol(ii+1,jj  ,k)*(0.5*fcent(ii+1,jj  ,k,1)-0.25) +
                 fvol(ii  ,jj+1,k)*(0.5*fcent(ii  ,jj+1,k,1)+0.25) +
                 fvol(ii+1,jj+1,k)*(0.5*fcent(ii+1,jj+1,k,1)+0.25));

            cba(i,j,k) = 0.5*(fba(ii,jj  ,k) + fba(ii+1,jj  ,k) +
                              fba(ii,jj+1,k) + fba(ii+1,jj+1,k));
            Real cbainv = 1.0/cba(i,j,k);

            cbc(i,j,k,0) = 0.5 * cbainv *
                (fba(ii  ,jj  ,k)*(0.5*fbc(ii  ,jj  ,k,0)-0.25) +
                 fba(ii+1,jj  ,k)*(0.5*fbc(ii+1,jj  ,k,0)+0.25) +
                 fba(ii  ,jj+1,k)*(0.5*fbc(ii  ,jj+1,k,0)-0.25) +
                 fba(ii+1,jj+1,k)*(0.5*fbc(ii+1,jj+1,k,0)+0.25));
            cbc(i,j,k,1) = 0.5 * cbainv *
                (fba(ii  ,jj  ,k)*(0.5*fbc(ii  ,jj  ,k,1)-0.25) +
                 fba(ii+1,jj  ,k)*(0.5*fbc(ii+1,jj  ,k,1)-0.25) +
                 fba(ii  ,jj+1,k)*(0.5*fbc(ii  ,jj+1,k,1)+0.25) +
                 fba(ii+1,jj+1,k)*(0.5*fbc(ii+1,jj+1,k,1)+0.25));

             Real nx = fbn(ii  ,jj  ,k,0)*fba(ii  ,jj  ,k)
                 +     fbn(ii+1,jj  ,k,0)*fba(ii+1,jj  ,k)
                 +     fbn(ii  ,jj+1,k,0)*fba(ii  ,jj+1,k)
                 +     fbn(ii+1,jj+1,k,0)*fba(ii+1,jj+1,k);
             Real ny = fbn(ii  ,jj  ,k,1)*fba(ii  ,jj  ,k)
                 +     fbn(ii+1,jj  ,k,1)*fba(ii+1,jj  ,k)
                 +     fbn(ii  ,jj+1,k,1)*fba(ii  ,jj+1,k)
                 +     fbn(ii+1,jj+1,k,1)*fba(ii+1,jj+1,k);
             Real nfac = 1.0/std::sqrt(nx*nx+ny*ny+1.e-50);
             cbn(i,j,k,0) = nx*nfac;
             cbn(i,j,k,1) = ny*nfac;
             ierr = (nx == 0.0 and ny == 0.0);
        }
    }
    else if (gbx.contains(iv))
    {
        cvol(i,j,k) = 1.0;
        ccent(i,j,k,0) = 0.0;
        ccent(i,j,k,1) = 0.0;
        cba(i,j,k) = 0.0;
        cbc(i,j,k,0) = -1.0;
        cbc(i,j,k,1) = -1.0;
        cbn(i,j,k,0) = 0.0;
        cbn(i,j,k,1) = 0.0;
    }

    if (xbx.contains(iv))
    {
        capx(i,j,k) = 0.5*(fapx(ii,jj,k) + fapx(ii,jj+1,k));
        if (capx(i,j,k) != 0.0) {
            cfcx(i,j,k) = 0.5 / capx(i,j,k) *
                (fapx(ii,jj  ,k)*(0.5*ffcx(ii,jj  ,k)-0.25) +
                 fapx(ii,jj+1,k)*(0.5*ffcx(ii,jj+1,k)+0.25));
        }
        else {
            cfcx(i,j,k) = 0.0;
        }
    }
    else if (xgbx.contains(iv))
    {
        capx(i,j,k) = 1.0;
        cfcx(i,j,k) = 0.0;
    }

    if (ybx.contains(iv))
    {
        capy(i,j,k) = 0.5*(fapy(ii,jj,k) + fapy(ii+1,jj,k));
        if (capy(i,j,k) != 0.0) {
            cfcy(i,j,k) = 0.5 / capy(i,j,k) *
                (fapy(ii  ,jj,k)*(0.5*ffcy(ii  ,jj,k)-0.25) +
                 fapy(ii+1,jj,k)*(0.5*ffcy(ii+1,jj,k)+0.25));
        } else {
            cfcy(i,j,k) = 0.0;
        }
    }
    else if (ygbx.contains(iv))
    {
        capy(i,j,k) = 1.0;
        cfcy(i,j,k) = 0.0;
    }

    return ierr;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void build_cellflag_from_ap (int i, int j, Array4<EBCellFlag> const& cflag,
                             Array4<Real const> const& apx, Array4<Real const> const& apy)
{
    constexpr int k = 0;

    // By default, all neighbors are already set.
    auto flg = cflag(i,j,k);

    if (apx(i  ,j  ,k) == 0.0) flg.setDisconnected(-1, 0, 0);
    if (apx(i+1,j  ,k) == 0.0) flg.setDisconnected( 1, 0, 0);
    if (apy(i  ,j  ,k) == 0.0) flg.setDisconnected( 0,-1, 0);
    if (apy(i  ,j+1,k) == 0.0) flg.setDisconnected( 0, 1, 0);

    if ((apx(i,j  ,k) == 0.0 or apy(i-1,j,k) == 0.0) and
        (apx(i,j-1,k) == 0.0 or apy(i  ,j,k) == 0.0))
    {
        flg.setDisconnected(-1,-1,0);
    }

    if ((apx(i+1,j  ,k) == 0.0 or apy(i+1,j,k) == 0.0) and
        (apx(i+1,j-1,k) == 0.0 or apy(i  ,j,k) == 0.0))
    {
        flg.setDisconnected(1,-1,0);
    }

    if ((apx(i,j  ,k) == 0.0 or apy(i-1,j+1,k) == 0.0) and
        (apx(i,j+1,k) == 0.0 or apy(i  ,j+1,k) == 0.0))
    {
        flg.setDisconnected(-1,1,0);
    }

    if ((apx(i+1,j  ,k) == 0.0 or apy(i+1,j+1,k) == 0.0) and
        (apx(i+1,j+1,k) == 0.0 or apy(i  ,j+1,k) == 0.0))
    {
        flg.setDisconnected(1,1,0);
    }

    cflag(i,j,k) = flg;
}

}}

#endif
