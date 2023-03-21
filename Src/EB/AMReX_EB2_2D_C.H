#ifndef AMREX_EB2_2D_C_H_
#define AMREX_EB2_2D_C_H_
#include <AMReX_Config.H>

namespace amrex::EB2 {

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
        if (    s(i,j  ,k) < 0.0_rt && s(i+1,j  ,k) < 0.0_rt
            && s(i,j+1,k) < 0.0_rt && s(i+1,j+1,k) < 0.0_rt)
        {
            cell(i,j,k).setRegular();
        }
        else if (s(i,j  ,k) >= 0.0_rt && s(i+1,j  ,k) >= 0.0_rt
            &&  s(i,j+1,k) >= 0.0_rt && s(i+1,j+1,k) >= 0.0_rt)
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
        if (s(i,j,k) < 0.0_rt && s(i,j+1,k) < 0.0_rt) {
            fx(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0_rt && s(i,j+1,k) >= 0.0_rt) {
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
        if (s(i,j,k) < 0.0_rt && s(i+1,j,k) < 0.0_rt) {
            fy(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0_rt && s(i+1,j,k) >= 0.0_rt) {
            fy(i,j,k) = Type::covered;
        } else {
            fy(i,j,k) = Type::irregular;
        }
    });
}

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int num_cuts (Real a, Real b) noexcept {
        return (a >= 0.0_rt && b < 0.0_rt) || (b >= 0.0_rt && a < 0.0_rt);
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
    return (ncuts != 0 && ncuts != 2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen_from_fine (int i, int j, Box const& bx, int ngrow,
                       Array4<Real> const& cvol, Array4<Real> const& ccent,
                       Array4<Real> const& cba, Array4<Real> const& cbc,
                       Array4<Real> const& cbn, Array4<Real> const& capx,
                       Array4<Real> const& capy, Array4<Real> const& cfcx,
                       Array4<Real> const& cfcy, Array4<Real> const& cecx,
                       Array4<Real> const& cecy, Array4<EBCellFlag> const& cflag,
                       Array4<Real const> const& fvol, Array4<Real const> const& fcent,
                       Array4<Real const> const& fba, Array4<Real const> const& fbc,
                       Array4<Real const> const& fbn, Array4<Real const> const& fapx,
                       Array4<Real const> const& fapy, Array4<Real const> const& ffcx,
                       Array4<Real const> const& ffcy, Array4<Real const> const& fecx,
                       Array4<Real const> const& fecy, Array4<EBCellFlag const> const& fflag)
{
    const Box& gbx = amrex::grow(bx,ngrow);
    const Box& xbx = amrex::surroundingNodes(bx,0);
    const Box& ybx = amrex::surroundingNodes(bx,1);
    const Box& xgbx = amrex::surroundingNodes(gbx,0);
    const Box& ygbx = amrex::surroundingNodes(gbx,1);

    int ierr = 0;
    constexpr int k = 0;
    int ii = i*2;
    int jj = j*2;

    if (bx.contains(i,j,k))
    {
        if (fflag(ii,jj  ,k).isRegular() && fflag(ii+1,jj  ,k).isRegular() &&
            fflag(ii,jj+1,k).isRegular() && fflag(ii+1,jj+1,k).isRegular())
        {
            cflag(i,j,k).setRegular();
            cvol(i,j,k) = 1.0_rt;
            ccent(i,j,k,0) = 0.0_rt;
            ccent(i,j,k,1) = 0.0_rt;
            cba(i,j,k) = 0.0_rt;
            cbc(i,j,k,0) = -1.0_rt;
            cbc(i,j,k,1) = -1.0_rt;
            cbn(i,j,k,0) = 0.0_rt;
            cbn(i,j,k,1) = 0.0_rt;
        }
        else if (fflag(ii,jj  ,k).isCovered() && fflag(ii+1,jj  ,k).isCovered() &&
                 fflag(ii,jj+1,k).isCovered() && fflag(ii+1,jj+1,k).isCovered())
        {
            cflag(i,j,k).setCovered();
            cvol(i,j,k) = 0.0_rt;
            ccent(i,j,k,0) = 0.0_rt;
            ccent(i,j,k,1) = 0.0_rt;
            cba(i,j,k) = 0.0_rt;
            cbc(i,j,k,0) = -1.0_rt;
            cbc(i,j,k,1) = -1.0_rt;
            cbn(i,j,k,0) = 0.0_rt;
            cbn(i,j,k,1) = 0.0_rt;
        }
        else
        {
            cflag(i,j,k).setSingleValued();

            cvol(i,j,k) = 0.25_rt*(fvol(ii,jj  ,k) + fvol(ii+1,jj  ,k) +
                                fvol(ii,jj+1,k) + fvol(ii+1,jj+1,k));
            Real cvolinv = 1.0_rt/cvol(i,j,k);

            ccent(i,j,k,0) = 0.25_rt * cvolinv *
                (fvol(ii  ,jj  ,k)*(0.5_rt*fcent(ii  ,jj  ,k,0)-0.25_rt) +
                 fvol(ii+1,jj  ,k)*(0.5_rt*fcent(ii+1,jj  ,k,0)+0.25_rt) +
                 fvol(ii  ,jj+1,k)*(0.5_rt*fcent(ii  ,jj+1,k,0)-0.25_rt) +
                 fvol(ii+1,jj+1,k)*(0.5_rt*fcent(ii+1,jj+1,k,0)+0.25_rt));
            ccent(i,j,k,1) = 0.25_rt * cvolinv *
                (fvol(ii  ,jj  ,k)*(0.5_rt*fcent(ii  ,jj  ,k,1)-0.25_rt) +
                 fvol(ii+1,jj  ,k)*(0.5_rt*fcent(ii+1,jj  ,k,1)-0.25_rt) +
                 fvol(ii  ,jj+1,k)*(0.5_rt*fcent(ii  ,jj+1,k,1)+0.25_rt) +
                 fvol(ii+1,jj+1,k)*(0.5_rt*fcent(ii+1,jj+1,k,1)+0.25_rt));

            cba(i,j,k) = 0.5_rt*(fba(ii,jj  ,k) + fba(ii+1,jj  ,k) +
                              fba(ii,jj+1,k) + fba(ii+1,jj+1,k));
            Real cbainv = 1.0_rt/cba(i,j,k);

            cbc(i,j,k,0) = 0.5_rt * cbainv *
                (fba(ii  ,jj  ,k)*(0.5_rt*fbc(ii  ,jj  ,k,0)-0.25_rt) +
                 fba(ii+1,jj  ,k)*(0.5_rt*fbc(ii+1,jj  ,k,0)+0.25_rt) +
                 fba(ii  ,jj+1,k)*(0.5_rt*fbc(ii  ,jj+1,k,0)-0.25_rt) +
                 fba(ii+1,jj+1,k)*(0.5_rt*fbc(ii+1,jj+1,k,0)+0.25_rt));
            cbc(i,j,k,1) = 0.5_rt * cbainv *
                (fba(ii  ,jj  ,k)*(0.5_rt*fbc(ii  ,jj  ,k,1)-0.25_rt) +
                 fba(ii+1,jj  ,k)*(0.5_rt*fbc(ii+1,jj  ,k,1)-0.25_rt) +
                 fba(ii  ,jj+1,k)*(0.5_rt*fbc(ii  ,jj+1,k,1)+0.25_rt) +
                 fba(ii+1,jj+1,k)*(0.5_rt*fbc(ii+1,jj+1,k,1)+0.25_rt));

             Real nx = fbn(ii  ,jj  ,k,0)*fba(ii  ,jj  ,k)
                 +     fbn(ii+1,jj  ,k,0)*fba(ii+1,jj  ,k)
                 +     fbn(ii  ,jj+1,k,0)*fba(ii  ,jj+1,k)
                 +     fbn(ii+1,jj+1,k,0)*fba(ii+1,jj+1,k);
             Real ny = fbn(ii  ,jj  ,k,1)*fba(ii  ,jj  ,k)
                 +     fbn(ii+1,jj  ,k,1)*fba(ii+1,jj  ,k)
                 +     fbn(ii  ,jj+1,k,1)*fba(ii  ,jj+1,k)
                 +     fbn(ii+1,jj+1,k,1)*fba(ii+1,jj+1,k);
             Real nfac = 1.0_rt/std::sqrt(nx*nx+ny*ny+1.e-30_rt);
             cbn(i,j,k,0) = nx*nfac;
             cbn(i,j,k,1) = ny*nfac;
             ierr = (nx == 0.0_rt && ny == 0.0_rt)
                 // we must check the enclosing surface to make sure the coarse cell does not
                 // fully contains the geometry object.
                 || ( fapx(ii,jj  ,k)==1.0_rt && fapx(ii+2,jj  ,k)==1.0_rt
                  && fapx(ii,jj+1,k)==1.0_rt && fapx(ii+2,jj+1,k)==1.0_rt
                  && fapy(ii,jj  ,k)==1.0_rt && fapy(ii+1,jj  ,k)==1.0_rt
                  && fapy(ii,jj+2,k)==1.0_rt && fapy(ii+1,jj+2,k)==1.0_rt);
        }
    }
    else if (gbx.contains(i,j,k))
    {
        cvol(i,j,k) = 1.0_rt;
        ccent(i,j,k,0) = 0.0_rt;
        ccent(i,j,k,1) = 0.0_rt;
        cba(i,j,k) = 0.0_rt;
        cbc(i,j,k,0) = -1.0_rt;
        cbc(i,j,k,1) = -1.0_rt;
        cbn(i,j,k,0) = 0.0_rt;
        cbn(i,j,k,1) = 0.0_rt;
    }

    if (xbx.contains(i,j,k))
    {
        capx(i,j,k) = 0.5_rt*(fapx(ii,jj,k) + fapx(ii,jj+1,k));
        if (capx(i,j,k) != 0.0_rt) {
            cfcx(i,j,k) = 0.5_rt / capx(i,j,k) *
                (fapx(ii,jj  ,k)*(0.5_rt*ffcx(ii,jj  ,k)-0.25_rt) +
                 fapx(ii,jj+1,k)*(0.5_rt*ffcx(ii,jj+1,k)+0.25_rt));
            if (fecy(ii,jj,k) == Real(1.0) && fecy(ii,jj+1,k) == Real(1.0)) {
                cecy(i,j,k) = Real(1.0);
            } else {
                cecy(i,j,k) = cfcx(i,j,k);
            }
        }
        else {
            cfcx(i,j,k) = 0.0_rt;
            cecy(i,j,k) = -1.0_rt;
        }
    }
    else if (xgbx.contains(i,j,k))
    {
        capx(i,j,k) = 1.0_rt;
        cfcx(i,j,k) = 0.0_rt;
        cecy(i,j,k) = 1.0_rt;
    }

    if (ybx.contains(i,j,k))
    {
        capy(i,j,k) = 0.5_rt*(fapy(ii,jj,k) + fapy(ii+1,jj,k));
        if (capy(i,j,k) != 0.0_rt) {
            cfcy(i,j,k) = 0.5_rt / capy(i,j,k) *
                (fapy(ii  ,jj,k)*(0.5_rt*ffcy(ii  ,jj,k)-0.25_rt) +
                 fapy(ii+1,jj,k)*(0.5_rt*ffcy(ii+1,jj,k)+0.25_rt));
            if (fecx(ii,jj,k) == Real(1.0) && fecx(ii+1,jj,k) == Real(1.0)) {
                cecx(i,j,k) = Real(1.0);
            } else {
                cecx(i,j,k) = cfcy(i,j,k);
            }
        } else {
            cfcy(i,j,k) = 0.0_rt;
            cecx(i,j,k) = -1.0_rt;
        }
    }
    else if (ygbx.contains(i,j,k))
    {
        capy(i,j,k) = 1.0_rt;
        cfcy(i,j,k) = 0.0_rt;
        cecx(i,j,k) = 1.0_rt;
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

    if (apx(i  ,j  ,k) == 0.0_rt) flg.setDisconnected(-1, 0, 0);
    if (apx(i+1,j  ,k) == 0.0_rt) flg.setDisconnected( 1, 0, 0);
    if (apy(i  ,j  ,k) == 0.0_rt) flg.setDisconnected( 0,-1, 0);
    if (apy(i  ,j+1,k) == 0.0_rt) flg.setDisconnected( 0, 1, 0);

    if ((apx(i,j  ,k) == 0.0_rt || apy(i-1,j,k) == 0.0_rt) &&
        (apx(i,j-1,k) == 0.0_rt || apy(i  ,j,k) == 0.0_rt))
    {
        flg.setDisconnected(-1,-1,0);
    }

    if ((apx(i+1,j  ,k) == 0.0_rt || apy(i+1,j,k) == 0.0_rt) &&
        (apx(i+1,j-1,k) == 0.0_rt || apy(i  ,j,k) == 0.0_rt))
    {
        flg.setDisconnected(1,-1,0);
    }

    if ((apx(i,j  ,k) == 0.0_rt || apy(i-1,j+1,k) == 0.0_rt) &&
        (apx(i,j+1,k) == 0.0_rt || apy(i  ,j+1,k) == 0.0_rt))
    {
        flg.setDisconnected(-1,1,0);
    }

    if ((apx(i+1,j  ,k) == 0.0_rt || apy(i+1,j+1,k) == 0.0_rt) &&
        (apx(i+1,j+1,k) == 0.0_rt || apy(i  ,j+1,k) == 0.0_rt))
    {
        flg.setDisconnected(1,1,0);
    }

    cflag(i,j,k) = flg;
}

}

#endif
