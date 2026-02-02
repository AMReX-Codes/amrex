#ifndef AMREX_MLEBABECLAP_2D_K_H_
#define AMREX_MLEBABECLAP_2D_K_H_
#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_EBData.H>

#include <AMReX_EB_LeastSquares_2D_K.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx_centroid (int i, int j, int k, int n, Array4<T> const& y,
                        Array4<T const> const& x, Array4<T const> const& a,
                        Array4<T const> const& bX, Array4<T const> const& bY,
                        EBData const& ebdata,
                        Array4<T const> const& beb,
                        Array4<T const> const& phieb,
                        const int& domlo_x,    const int& domlo_y,
                        const int& domhi_x,    const int& domhi_y,
                        const bool& on_x_face, const bool& on_y_face,
                        bool is_eb_dirichlet, bool is_eb_inhomog,
                        GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
                        T alpha, T beta) noexcept
{
    T dhx = beta*dxinv[0]*dxinv[0];
    T dhy = beta*dxinv[1]*dxinv[1];
    T dh  = beta*dxinv[0]*dxinv[0];

    Array4<EBCellFlag const> const& flag = ebdata.get<EBData_t::cellflag>();
    Array4<Real const> const& vfrc = ebdata.get<EBData_t::volfrac>();

    Array4<Real const> const& apx = ebdata.get<EBData_t::apx>();
    Array4<Real const> const& apy = ebdata.get<EBData_t::apy>();

    Array4<Real const> const& fcx = ebdata.get<EBData_t::fcx>();
    Array4<Real const> const& fcy = ebdata.get<EBData_t::fcy>();

    Array4<Real const> const& ccent = ebdata.get<EBData_t::centroid>();
    Array4<Real const> const& bcent = ebdata.get<EBData_t::bndrycent>();
    Array4<Real const> const& ba = ebdata.get<EBData_t::bndryarea>();

    if (flag(i,j,k).isCovered())
    {
        y(i,j,k,n) = T(0.0);
    }
    else if (flag(i,j,k).isRegular() &&
            ((flag(i-1,j  ,k).isRegular() && flag(i+1,j  ,k).isRegular() &&
             flag(i  ,j-1,k).isRegular() && flag(i  ,j+1,k).isRegular()) ))
    {
        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
            - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                   - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
            - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                   - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)));

    }
    else
    {
        T kappa = vfrc(i,j,k);
        T apxm = apx(i,j,k);
        T apxp = apx(i+1,j,k);
        T apym = apy(i,j,k);
        T apyp = apy(i,j+1,k);

        // First get EB-aware slope that doesn't know about extdir
        bool needs_bdry_stencil = (i <= domlo_x) || (i >= domhi_x) ||
                                  (j <= domlo_y) || (j >= domhi_y);

        // if phi_on_centroid -- A second order least squares fit is used
        // to approximate the slope on the high and low faces. Note that if
        // any of the three cells --e.g., (i-1,j), (i,j), or (i-1,j)-- are
        // cut, then the least squares fit is needed. This is a bit more than
        // is actually needed for most cases but it will return the correct
        // value in all cases.

        T fxm = bX(i,j,k,n) * (x(i,j,k,n)-x(i-1,j,k,n));
        if ( (apxm != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i-1,j,k) != T(1.0) || vfrc(i+1,j,k) != T(1.0)) )
        {
            T yloc_on_xface = fcx(i,j,k);

            if(needs_bdry_stencil) {

              fxm = grad_x_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                      yloc_on_xface,is_eb_dirichlet,is_eb_inhomog,
                                                      on_x_face, domlo_x, domhi_x,
                                                      on_y_face, domlo_y, domhi_y);

            } else {
              fxm = grad_x_of_phi_on_centroids(i,j,k,n,x,phieb,flag,ccent,bcent,
                                               yloc_on_xface,is_eb_dirichlet,is_eb_inhomog);
            }
            fxm *= bX(i,j,k,n);
        }

        T fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n)-x(i,j,k,n));
        if ( (apxp != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i+1,j,k) != T(1.0) || vfrc(i-1,j,k) != T(1.0)) ) {
            T yloc_on_xface = fcx(i+1,j,k,0);
            if(needs_bdry_stencil) {
              fxp = grad_x_of_phi_on_centroids_extdir(i+1,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                      yloc_on_xface,is_eb_dirichlet,is_eb_inhomog,
                                                      on_x_face, domlo_x, domhi_x,
                                                      on_y_face, domlo_y, domhi_y);

            } else {
              fxp = grad_x_of_phi_on_centroids(i+1,j,k,n,x,phieb,flag,ccent,bcent,
                                               yloc_on_xface,is_eb_dirichlet,is_eb_inhomog);
            }
            fxp *= bX(i+1,j,k,n);

        }

        T fym = bY(i,j,k,n)*(x(i,j,k,n)-x(i,j-1,k,n));
        if ( (apym != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j-1,k) != T(1.0) || vfrc(i,j+1,k) != T(1.0)) ) {
            T xloc_on_yface = fcy(i,j,k,0);

            if(needs_bdry_stencil) {

              fym = grad_y_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                      xloc_on_yface,is_eb_dirichlet,is_eb_inhomog,
                                                      on_x_face, domlo_x, domhi_x,
                                                      on_y_face, domlo_y, domhi_y);

            } else {
              fym = grad_y_of_phi_on_centroids(i,j,k,n,x,phieb,flag,ccent,bcent,
                                               xloc_on_yface,is_eb_dirichlet,is_eb_inhomog);
            }
            fym *= bY(i,j,k,n);
        }

        T fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n)-x(i,j,k,n));
        if ( (apyp != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j+1,k) != T(1.0) || vfrc(i,j-1,k) != T(1.0)) ) {
            T xloc_on_yface = fcy(i,j+1,k,0);
            if(needs_bdry_stencil) {
              fyp = grad_y_of_phi_on_centroids_extdir(i,j+1,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                      xloc_on_yface,is_eb_dirichlet,is_eb_inhomog,
                                                      on_x_face, domlo_x, domhi_x,
                                                      on_y_face, domlo_y, domhi_y);

            } else {
              fyp = grad_y_of_phi_on_centroids(i,j+1,k,n,x,phieb,flag,ccent,bcent,
                                               xloc_on_yface,is_eb_dirichlet,is_eb_inhomog);
            }
            fyp *= bY(i,j+1,k,n);
        }

        T feb = T(0.0);
        if (is_eb_dirichlet && flag(i,j,k).isSingleValued())
        {
            T dapx = (apxm-apxp) / dxinv[0];
            T dapy = (apym-apyp) / dxinv[1];
            T anorm = std::hypot(dapx,dapy);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;

            feb = grad_eb_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                     anrmx,anrmy,is_eb_inhomog,
                                                     on_x_face, domlo_x, domhi_x,
                                                     on_y_face, domlo_y, domhi_y);

            feb *= ba(i,j,k) * beb(i,j,k,n);
        }

        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (T(1.0)/kappa) *
            (dhx*(apxm*fxm-apxp*fxp) + dhy*(apym*fym-apyp*fyp) - dh*feb);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx (int i, int j, int k, int n, Array4<T> const& y,
                        Array4<T const> const& x, Array4<T const> const& a,
                        Array4<T const> const& bX, Array4<T const> const& bY,
                        Array4<const int> const& ccm,
                        EBData const& ebdata,
                        Array4<T const> const& beb,
                        bool is_dirichlet, Array4<T const> const& phieb,
                        bool is_inhomog, T alpha,T beta,
                        GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                        bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    T dhx = beta*dxinv[0]*dxinv[0];
    T dhy = beta*dxinv[1]*dxinv[1];
    T dh  = beta*dxinv[0]*dxinv[1];

    bool beta_on_center = !(beta_on_centroid);
    bool  phi_on_center = !( phi_on_centroid);
    auto const flag = ebdata.get<EBData_t::cellflag>(i,j,k);

    if (flag.isCovered())
    {
        y(i,j,k,n) = T(0.0);
    }
    else if (flag.isRegular())
    {
        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
            - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                    - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
            - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                    - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)));
    }
    else
    {
        T kappa = ebdata.get<EBData_t::volfrac>(i,j,k);
        T apxm = ebdata.get<EBData_t::apx>(i  ,j  ,k  );
        T apxp = ebdata.get<EBData_t::apx>(i+1,j  ,k  );
        T apym = ebdata.get<EBData_t::apy>(i  ,j  ,k  );
        T apyp = ebdata.get<EBData_t::apy>(i  ,j+1,k  );

        T fxm = bX(i,j,k,n) * (x(i,j,k,n)-x(i-1,j,k,n));
        if (apxm != T(0.0) && apxm != T(1.0)) {
            T fcx = ebdata.get<EBData_t::fcx>(i,j,k);
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx) : T(0.0);
            if (beta_on_center && phi_on_center) {
                fxm = (T(1.0)-fracy)*fxm + fracy*bX(i,jj,k,n)*(x(i,jj,k,n)-x(i-1,jj,k,n));
            } else if (beta_on_centroid && phi_on_center) {
                fxm = bX(i,j,k,n) * ( (T(1.0)-fracy)*(x(i, j,k,n)-x(i-1, j,k,n))
                                        + fracy  *(x(i,jj,k,n)-x(i-1,jj,k,n)) );
            }
        }

        T fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n)-x(i,j,k,n));
        if (apxp != T(0.0) && apxp != T(1.0)) {
            T fcx = ebdata.get<EBData_t::fcx>(i+1,j,k);
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx));
            T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx) : T(0.0);
            if (beta_on_center && phi_on_center) {
                fxp = (T(1.0)-fracy)*fxp + fracy*bX(i+1,jj,k,n)*(x(i+1,jj,k,n)-x(i,jj,k,n));
            } else if (beta_on_centroid && phi_on_center) {
                fxp = bX(i+1,j,k,n) * ( (T(1.0)-fracy)*(x(i+1, j,k,n)-x(i, j,k,n))
                                            + fracy *(x(i+1,jj,k,n)-x(i,jj,k,n)) );
            }
        }

        T fym = bY(i,j,k,n)*(x(i,j,k,n)-x(i,j-1,k,n));
        if (apym != T(0.0) && apym != T(1.0)) {
            T fcy = ebdata.get<EBData_t::fcy>(i,j,k);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy) : T(0.0);
            if (beta_on_center && phi_on_center) {
                fym = (T(1.0)-fracx)*fym + fracx*bY(ii,j,k,n)*(x(ii,j,k,n)-x(ii,j-1,k,n));
            } else if (beta_on_centroid && phi_on_center) {
                fym = bY(i,j,k,n) * ( (T(1.0)-fracx)*(x( i,j,k,n)-x( i,j-1,k,n))
                                            + fracx *(x(ii,j,k,n)-x(ii,j-1,k,n)) );
            }
        }

        T fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n)-x(i,j,k,n));
        if (apyp != T(0.0) && apyp != T(1.0)) {
            T fcy = ebdata.get<EBData_t::fcy>(i,j+1,k);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy));
            T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy) : T(0.0);
            if (beta_on_center && phi_on_center) {
                fyp = (T(1.0)-fracx)*fyp + fracx*bY(ii,j+1,k,n)*(x(ii,j+1,k,n)-x(ii,j,k,n));
            } else if (beta_on_centroid && phi_on_center) {
                fyp = bY(i,j+1,k,n) * ( (T(1.0)-fracx)*(x( i,j+1,k,n)-x( i,j,k,n))
                                            + fracx *(x(ii,j+1,k,n)-x(ii,j,k,n)) );
            }
        }

        T feb = T(0.0);
        if (is_dirichlet) {
            T dapx = (apxm-apxp) / dxinv[0];
            T dapy = (apym-apyp) / dxinv[1];
            T anorm = std::hypot(dapx,dapy);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;

            T phib = is_inhomog ? phieb(i,j,k,n) : T(0.0);

            T bctx = ebdata.get<EBData_t::bndrycent>(i,j,k,0);
            T bcty = ebdata.get<EBData_t::bndrycent>(i,j,k,1);
            T dx_eb = get_dx_eb(kappa);

            T dg, gx, gy, sx, sy;
            if (std::abs(anrmx) > std::abs(anrmy)) {
                dg = dx_eb / std::abs(anrmx);
            } else {
                dg = dx_eb / std::abs(anrmy);
            }
            gx = (bctx - dg*anrmx);
            gy = (bcty - dg*anrmy);
            sx = std::copysign(T(1.0),anrmx);
            sy = std::copysign(T(1.0),anrmy);

            int ii = i - static_cast<int>(sx);
            int jj = j - static_cast<int>(sy);

            T phig = (T(1.0) + gx*sx + gy*sy + gx*gy*sx*sy) * x(i ,j ,k,n)
                +       (    - gx*sx         - gx*gy*sx*sy) * x(ii,j ,k,n)
                +       (            - gy*sy - gx*gy*sx*sy) * x(i ,jj,k,n)
                +       (                    + gx*gy*sx*sy) * x(ii,jj,k,n) ;

            T dphidn = (phib-phig) / dg;

            T ba = ebdata.get<EBData_t::bndryarea>(i,j,k);

            feb = dphidn * ba * beb(i,j,k,n);
        }

        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (T(1.0)/kappa) *
            (dhx*(apxm*fxm-apxp*fxp) +
             dhy*(apym*fym-apyp*fyp) - dh*feb);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_ebflux (int i, int j, int k, int n,
                         Array4<T> const& feb,
                         Array4<T const> const& x,
                         Array4<EBCellFlag const> const& flag,
                         Array4<T const> const& vfrc,
                         Array4<T const> const& apx,
                         Array4<T const> const& apy,
                         Array4<T const> const& bc,
                         Array4<T const> const& beb,
                         Array4<T const> const& phieb,
                         bool is_inhomog,
                         GpuArray<T,AMREX_SPACEDIM> const& dxinv) noexcept
{

    if (!flag(i,j,k).isSingleValued())
    {
        feb(i,j,k,n) = T(0.0);
    }
    else
    {
        T kappa = vfrc(i,j,k);
        T apxm = apx(i,j,k);
        T apxp = apx(i+1,j,k);
        T apym = apy(i,j,k);
        T apyp = apy(i,j+1,k);

        T dapx = (apxm-apxp)/dxinv[1];
        T dapy = (apym-apyp)/dxinv[0];
        T anorm = std::hypot(dapx,dapy);
        T anorminv = T(1.0)/anorm;
        T anrmx = dapx * anorminv;
        T anrmy = dapy * anorminv;
        const T bareascaling = std::sqrt( (anrmx/dxinv[0])*(anrmx/dxinv[0]) +
                (anrmy/dxinv[1])*(anrmy/dxinv[1]) );


        T phib = is_inhomog ? phieb(i,j,k,n) : T(0.0);

        T bctx = bc(i,j,k,0);
        T bcty = bc(i,j,k,1);
        T dx_eb = get_dx_eb(kappa);

        T dg, gx, gy, sx, sy;
        if (std::abs(anrmx) > std::abs(anrmy)) {
            dg = dx_eb / std::abs(anrmx);
        } else {
            dg = dx_eb / std::abs(anrmy);
        }
        gx = bctx - dg*anrmx;
        gy = bcty - dg*anrmy;
        sx = std::copysign(T(1.0),anrmx);
        sy = std::copysign(T(1.0),anrmy);

        int ii = i - static_cast<int>(sx);
        int jj = j - static_cast<int>(sy);

        T phig = (T(1.0) + gx*sx + gy*sy + gx*gy*sx*sy) * x(i ,j ,k,n)
            +       (    - gx*sx         - gx*gy*sx*sy) * x(ii,j ,k,n)
            +       (            - gy*sy - gx*gy*sx*sy) * x(i ,jj,k,n)
            +       (                    + gx*gy*sx*sy) * x(ii,jj,k,n) ;

        T dphidn = (phib-phig)/(dg * bareascaling);
        feb(i,j,k,n) = -beb(i,j,k,n) * dphidn;
    }
}


template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_gsrb (int i, int j, int k, int n,
                       Array4<T> const& phi, Array4<T const> const& rhs,
                       T alpha, Array4<T const> const& a, T beta,
                       GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                       Array4<T const> const& bX, Array4<T const> const& bY,
                       Array4<int const> const& m0, Array4<int const> const& m2,
                       Array4<int const> const& m1, Array4<int const> const& m3,
                       Array4<T const> const& f0, Array4<T const> const& f2,
                       Array4<T const> const& f1, Array4<T const> const& f3,
                       Array4<const int> const& ccm, Array4<T const> const& beb,
                       EBData const& ebdata,
                       bool is_dirichlet, bool beta_on_centroid, bool phi_on_centroid,
                       Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);
    T dhx = beta*dxinv[0]*dxinv[0];
    T dhy = beta*dxinv[1]*dxinv[1];
    T dh  = beta*dxinv[0]*dxinv[1];

    if ((i+j+k+redblack) % 2 == 0)
    {
        auto const flag = ebdata.get<EBData_t::cellflag>(i,j,k);
        if (flag.isCovered())
        {
            phi(i,j,k,n) = T(0.0);
        }
        else
        {
            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                ? f0(vlo.x,j,k,n) : T(0.0);
            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                ? f1(i,vlo.y,k,n) : T(0.0);
            T cf2 = (i == vhi.x && m2(vhi.x+1,j,k) > 0)
                ? f2(vhi.x,j,k,n) : T(0.0);
            T cf3 = (j == vhi.y && m3(i,vhi.y+1,k) > 0)
                ? f3(i,vhi.y,k,n) : T(0.0);

            if (flag.isRegular())
            {
                T gamma = alpha*a(i,j,k)
                    + dhx * (bX(i+1,j,k,n) + bX(i,j,k,n))
                    + dhy * (bY(i,j+1,k,n) + bY(i,j,k,n));

                T rho =  dhx * (bX(i+1,j,k,n)*phi(i+1,j,k,n)
                              + bX(i  ,j,k,n)*phi(i-1,j,k,n))
                       + dhy * (bY(i,j+1,k,n)*phi(i,j+1,k,n)
                              + bY(i,j  ,k,n)*phi(i,j-1,k,n));

                T delta = dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf2)
                 +        dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf3);

                T res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                phi(i,j,k,n) += res/(gamma-delta);
            }
            else
            {
                T kappa = ebdata.get<EBData_t::volfrac>(i,j,k);
                T apxm = ebdata.get<EBData_t::apx>(i  ,j  ,k);
                T apxp = ebdata.get<EBData_t::apx>(i+1,j  ,k);
                T apym = ebdata.get<EBData_t::apy>(i  ,j  ,k);
                T apyp = ebdata.get<EBData_t::apy>(i  ,j+1,k);

                T fxm = -bX(i,j,k,n)*phi(i-1,j,k,n);
                T oxm = -bX(i,j,k,n)*cf0;
                T sxm =  bX(i,j,k,n);
                if (apxm != T(0.0) && apxm != T(1.0)) {
                    T fcx = ebdata.get<EBData_t::fcx>(i,j,k);
                    int jj = j + static_cast<int>(std::copysign(T(1.0),fcx));
                    T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fxm = (T(1.0)-fracy)*fxm +
                                    fracy *bX(i,jj,k,n)*(phi(i,jj,k,n)-phi(i-1,jj,k,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fxm = (T(1.0)-fracy)*(             -phi(i-1,j,k,n)) +
                                      fracy *(phi(i,jj,k,n)-phi(i-1,jj,k,n));
                        fxm *= bX(i,j,k,n);
                    }
                    oxm = T(0.0);
                    sxm = (T(1.0)-fracy)*sxm;
                }

                T fxp =  bX(i+1,j,k,n)*phi(i+1,j,k,n);
                T oxp =  bX(i+1,j,k,n)*cf2;
                T sxp = -bX(i+1,j,k,n);
                if (apxp != T(0.0) && apxp != T(1.0)) {
                    T fcx = ebdata.get<EBData_t::fcx>(i+1,j,k);
                    int jj = j + static_cast<int>(std::copysign(T(1.0),fcx));
                    T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fxp = (T(1.0)-fracy)*fxp +
                                    fracy *bX(i+1,jj,k,n)*(phi(i+1,jj,k,n)-phi(i,jj,k,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fxp = (T(1.0)-fracy)*(phi(i+1,j,k,n)               ) +
                                      fracy *(phi(i+1,jj,k,n)-phi(i,jj,k,n));
                        fxp *= bX(i+1,j,k,n);
                    }
                    oxp = T(0.0);
                    sxp = (T(1.0)-fracy)*sxp;
                }

                T fym = -bY(i,j,k,n)*phi(i,j-1,k,n);
                T oym = -bY(i,j,k,n)*cf1;
                T sym =  bY(i,j,k,n);
                if (apym != T(0.0) && apym != T(1.0)) {
                    T fcy = ebdata.get<EBData_t::fcy>(i,j,k);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcy));
                    T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fym = (T(1.0)-fracx)*fym +
                                    fracx *bY(ii,j,k,n)*(phi(ii,j,k,n)-phi(ii,j-1,k,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fym = (T(1.0)-fracx)*(             -phi( i,j-1,k,n)) +
                                      fracx *(phi(ii,j,k,n)-phi(ii,j-1,k,n));
                        fym *= bY(i,j,k,n);
                    }

                    oym = T(0.0);
                    sym = (T(1.0)-fracx)*sym;
                }

                T fyp =  bY(i,j+1,k,n)*phi(i,j+1,k,n);
                T oyp =  bY(i,j+1,k,n)*cf3;
                T syp = -bY(i,j+1,k,n);
                if (apyp != T(0.0) && apyp != T(1.0)) {
                    T fcy = ebdata.get<EBData_t::fcy>(i,j+1,k);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcy));
                    T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fyp = (T(1.0)-fracx)*fyp +
                                    fracx*bY(ii,j+1,k,n)*(phi(ii,j+1,k,n)-phi(ii,j,k,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fyp = (T(1.0)-fracx)*(phi(i,j+1,k,n)               )+
                                      fracx *(phi(ii,j+1,k,n)-phi(ii,j,k,n));
                        fyp *= bY(i,j+1,k,n);
                    }
                    oyp = T(0.0);
                    syp = (T(1.0)-fracx)*syp;
                }

                T vfrcinv = (T(1.0)/kappa);
                T gamma = alpha*a(i,j,k) + vfrcinv *
                    (dhx*(apxm*sxm-apxp*sxp) +
                        dhy*(apym*sym-apyp*syp));
                T rho = -vfrcinv *
                    (dhx*(apxm*fxm-apxp*fxp) +
                        dhy*(apym*fym-apyp*fyp));

                T delta = -vfrcinv *
                    (dhx*(apxm*oxm-apxp*oxp) +
                        dhy*(apym*oym-apyp*oyp));

                if (is_dirichlet) {
                    T dapx = (apxm-apxp) / dxinv[0];
                    T dapy = (apym-apyp) / dxinv[1];
                    T anorm = std::hypot(dapx,dapy);
                    T anorminv = T(1.0)/anorm;
                    T anrmx = dapx * anorminv;
                    T anrmy = dapy * anorminv;

                    T bctx = ebdata.get<EBData_t::bndrycent>(i,j,k,0);
                    T bcty = ebdata.get<EBData_t::bndrycent>(i,j,k,1);
                    T dx_eb = get_dx_eb(kappa);

                    T dg, gx, gy, sx, sy;
                    if (std::abs(anrmx) > std::abs(anrmy)) {
                        dg = dx_eb / std::abs(anrmx);
                    } else {
                        dg = dx_eb / std::abs(anrmy);
                    }
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = std::copysign(T(1.0),anrmx);
                    sy = std::copysign(T(1.0),anrmy);

                    int ii = i - static_cast<int>(sx);
                    int jj = j - static_cast<int>(sy);

                    T phig_gamma = (T(1.0) + gx*sx + gy*sy + gx*gy*sx*sy);
                    T phig = (    - gx*sx         - gx*gy*sx*sy) * phi(ii,j ,k,n)
                     +       (            - gy*sy - gx*gy*sx*sy) * phi(i ,jj,k,n)
                     +       (                    + gx*gy*sx*sy) * phi(ii,jj,k,n);

                    // In gsrb we are always in residual-correction form so phib = 0
                    T dphidn =  (    -phig)/dg;

                    T ba = ebdata.get<EBData_t::bndryarea>(i,j,k);

                    T feb = dphidn * ba * beb(i,j,k,n);
                    rho += -vfrcinv*(-dh)*feb;

                    T feb_gamma = -phig_gamma/dg * ba * beb(i,j,k,n);
                    gamma += vfrcinv*(-dh)*feb_gamma;
                }

                T res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                phi(i,j,k,n) += res/(gamma-delta);
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_x (Box const& box, Array4<T> const& fx, Array4<T const> const& apx,
                         Array4<T const> const& fcx, Array4<T const> const& sol,
                         Array4<T const> const& bX, Array4<int const> const& ccm,
                         T dhx, int face_only, int ncomp, Box const& xbox,
                         bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == i || hif == i) {
            if (apx(i,j,k) == T(0.0)) {
                fx(i,j,k,n) = T(0.0);
            } else if (apx(i,j,k) == T(1.0)) {
                fx(i,j,k,n) = -dhx*bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            } else {
                T fxm = bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
                int jj = j + static_cast<int>(std::copysign(T(1.0),fcx(i,j,k)));
                T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k)) : T(0.0);
                if (!beta_on_centroid && !phi_on_centroid)
                {
                    fxm = (T(1.0)-fracy)*fxm + fracy*bX(i,jj,k,n)*(sol(i,jj,k,n)-sol(i-1,jj,k,n));
                } else if (beta_on_centroid && !phi_on_centroid) {
                    fxm = bX(i,j,k,n) * ( (T(1.0)-fracy)*(sol(i, j,k,n)-sol(i-1, j,k,n)) +
                                               fracy *(sol(i,jj,k,n)-sol(i-1,jj,k,n)) );
                }
                fx(i,j,k,n) = -fxm*dhx;
            }
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_y (Box const& box, Array4<T> const& fy, Array4<T const> const& apy,
                         Array4<T const> const& fcy, Array4<T const> const& sol,
                         Array4<T const> const& bY, Array4<int const> const& ccm,
                         T dhy, int face_only, int ncomp, Box const& ybox,
                         bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == j || hif == j) {
            if (apy(i,j,k) == T(0.0)) {
                fy(i,j,k,n) = T(0.0);
            } else if (apy(i,j,k) == T(1.0)) {
                fy(i,j,k,n) = -dhy*bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            } else {
                T fym = bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
                int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k)));
                T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k)) : T(0.0);
                if (!beta_on_centroid && !phi_on_centroid)
                {
                    fym = (T(1.0)-fracx)*fym + fracx*bY(ii,j,k,n)*(sol(ii,j,k,n)-sol(ii,j-1,k,n));
                } else if (beta_on_centroid && !phi_on_centroid) {
                    fym = bY(i,j,k,n) * ( (T(1.0)-fracx)*(sol( i,j,k,n)-sol( i,j-1,k,n)) +
                                               fracx *(sol(ii,j,k,n)-sol(ii,j-1,k,n)) );
                }
                fy(i,j,k,n) = -fym*dhy;
            }
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_x_0 (Box const& box, Array4<T> const& fx, Array4<T const> const& apx,
                           Array4<T const> const& sol, Array4<T const> const& bX,
                           T dhx, int face_only, int ncomp, Box const& xbox) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == i || hif == i) {
            if (apx(i,j,k) == T(0.0)) {
                fx(i,j,k,n) = T(0.0);
            } else {
                fx(i,j,k,n) = -dhx*bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            }
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_y_0 (Box const& box, Array4<T> const& fy, Array4<T const> const& apy,
                           Array4<T const> const& sol, Array4<T const> const& bY,
                           T dhy, int face_only, int ncomp, Box const& ybox) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == j || hif == j) {
            if (apy(i,j,k) == T(0.0)) {
                fy(i,j,k,n) = T(0.0);
            } else {
                fy(i,j,k,n) = -dhy*bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_x (Box const& box, Array4<T> const& gx, Array4<T const> const& sol,
                         Array4<T const> const& apx, Array4<T const> const& fcx,
                         Array4<int const> const& ccm,
                         T dxi, int ncomp, bool phi_on_centroid) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) == T(0.0)) {
            gx(i,j,k,n) = T(0.0);
        } else if (apx(i,j,k) == T(1.0)) {
            gx(i,j,k,n) = dxi*(sol(i,j,k,n)-sol(i-1,j,k,n));
        } else {
            T gxm = (sol(i,j,k,n)-sol(i-1,j,k,n));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx(i,j,k)));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k)) : T(0.0);
            if (!phi_on_centroid) {
               gxm = (T(1.0)-fracy)*gxm + fracy*(sol(i,jj,k,n)-sol(i-1,jj,k,n));
            }
            gx(i,j,k,n) = gxm*dxi;
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_y (Box const& box, Array4<T> const& gy, Array4<T const> const& sol,
                         Array4<T const> const& apy, Array4<T const> const& fcy,
                         Array4<int const> const& ccm,
                         T dyi, int ncomp, bool phi_on_centroid) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == T(0.0)) {
            gy(i,j,k,n) = T(0.0);
        } else if (apy(i,j,k) == T(1.0)) {
            gy(i,j,k,n) = dyi*(sol(i,j,k,n)-sol(i,j-1,k,n));
        } else {
            T gym = (sol(i,j,k,n)-sol(i,j-1,k,n));
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k)));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k)) : T(0.0);
            if (!phi_on_centroid) {
                gym = (T(1.0)-fracx)*gym + fracx*(sol(ii,j,k,n)-sol(ii,j-1,k,n));
            }
            gy(i,j,k,n) = gym*dyi;
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_x_0 (Box const& box, Array4<T> const& gx, Array4<T const> const& sol,
                           Array4<T const> const& apx, T dxi, int ncomp) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) == T(0.0)) {
            gx(i,j,k,n) = T(0.0);
        } else {
            gx(i,j,k,n) = dxi*(sol(i,j,k,n)-sol(i-1,j,k,n));
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_y_0 (Box const& box, Array4<T> const& gy, Array4<T const> const& sol,
                           Array4<T const> const& apy, T dyi, int ncomp) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == T(0.0)) {
            gy(i,j,k,n) = T(0.0);
        } else {
            gy(i,j,k,n) = dyi*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_normalize (int i, int j, int k, int n, Array4<T> const& phi,
                            T alpha, Array4<T const> const& a,
                            T dhx, T dhy, T dh,
                            const amrex::GpuArray<T, AMREX_SPACEDIM>& dx,
                            Array4<T const> const& bX, Array4<T const> const& bY,
                            Array4<const int> const& ccm, Array4<EBCellFlag const> const& flag,
                            Array4<T const> const& vfrc,
                            Array4<T const> const& apx, Array4<T const> const& apy,
                            Array4<T const> const& fcx, Array4<T const> const& fcy,
                            Array4<T const> const& ba, Array4<T const> const& bc,
                            Array4<T const> const& beb,
                            bool is_dirichlet, bool beta_on_centroid) noexcept
{
    if (flag(i,j,k).isRegular())
    {
        phi(i,j,k,n) /= alpha*a(i,j,k) + dhx*(bX(i,j,k,n) + bX(i+1,j,k,n))
                                       + dhy*(bY(i,j,k,n) + bY(i,j+1,k,n));
    }
    else if (flag(i,j,k).isSingleValued())
    {
        T kappa = vfrc(i,j,k);
        T apxm = apx(i,j,k);
        T apxp = apx(i+1,j,k);
        T apym = apy(i,j,k);
        T apyp = apy(i,j+1,k);

        T sxm =  bX(i,j,k,n);
        if (apxm != T(0.0) && apxm != T(1.0) && !beta_on_centroid) {
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx(i,j,k)));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                ? std::abs(fcx(i,j,k)) : T(0.0);
            sxm = (T(1.0)-fracy)*sxm;
        }

        T sxp = -bX(i+1,j,k,n);
        if (apxp != T(0.0) && apxp != T(1.0) && !beta_on_centroid) {
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx(i+1,j,k)));
            T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                ? std::abs(fcx(i+1,j,k)) : T(0.0);
            sxp = (T(1.0)-fracy)*sxp;
        }

        T sym =  bY(i,j,k,n);
        if (apym != T(0.0) && apym != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k)));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                ? std::abs(fcy(i,j,k)) : T(0.0);
            sym = (T(1.0)-fracx)*sym;
        }

        T syp = -bY(i,j+1,k,n);
        if (apyp != T(0.0) && apyp != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j+1,k)));
            T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                ? std::abs(fcy(i,j+1,k)) : T(0.0);
            syp = (T(1.0)-fracx)*syp;
        }

        T vfrcinv = (T(1.0)/kappa);
        T gamma = alpha*a(i,j,k) + vfrcinv *
            (dhx*(apxm*sxm-apxp*sxp) +
             dhy*(apym*sym-apyp*syp));

        if (is_dirichlet) {
            T dapx = (apxm-apxp)*dx[1];
            T dapy = (apym-apyp)*dx[0];
            T anorm = std::hypot(dapx,dapy);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;

            T bctx = bc(i,j,k,0);
            T bcty = bc(i,j,k,1);
            T dx_eb = get_dx_eb(vfrc(i,j,k));

            T dg, gx, gy, sx, sy;
            if (std::abs(anrmx) > std::abs(anrmy)) {
                dg = dx_eb / std::abs(anrmx);
            } else {
                dg = dx_eb / std::abs(anrmy);
            }
            gx = bctx - dg*anrmx;
            gy = bcty - dg*anrmy;
            sx = std::copysign(T(1.0),anrmx);
            sy = std::copysign(T(1.0),anrmy);

            T phig_gamma = (T(1.0) + gx*sx + gy*sy + gx*gy*sx*sy);
            T feb_gamma = -phig_gamma/dg * ba(i,j,k) * beb(i,j,k,n);
            gamma += vfrcinv*(-dh)*feb_gamma;
        }
        phi(i,j,k,n) /= gamma;
    }
}

}
#endif
