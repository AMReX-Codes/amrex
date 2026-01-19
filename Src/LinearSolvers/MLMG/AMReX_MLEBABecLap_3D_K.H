#ifndef AMREX_MLEBABECLAP_3D_K_H_
#define AMREX_MLEBABECLAP_3D_K_H_
#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_EBData.H>

#include <AMReX_EB_LeastSquares_3D_K.H>

namespace amrex {

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx_centroid (int i, int j, int k, int n, Array4<T> const& y,
                        Array4<T const> const& x, Array4<T const> const& a,
                        Array4<T const> const& bX, Array4<T const> const& bY,
                        Array4<T const> const& bZ,
                        EBData const& ebdata,
                        Array4<T const> const& beb,
                        Array4<T const> const& phieb,
                        const int& domlo_x, const int& domlo_y, const int& domlo_z,
                        const int& domhi_x, const int& domhi_y, const int& domhi_z,
                        const bool& on_x_face, const bool& on_y_face, const bool& on_z_face,
                        bool is_eb_dirichlet, bool is_eb_inhomog,
                        GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
                        T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];

    Array4<EBCellFlag const> const& flag = ebdata.get<EBData_t::cellflag>();
    Array4<Real const>  const& vfrc = ebdata.get<EBData_t::volfrac>();

    Array4<Real const> const& apx = ebdata.get<EBData_t::apx>();
    Array4<Real const> const& apy = ebdata.get<EBData_t::apy>();
    Array4<Real const> const& apz = ebdata.get<EBData_t::apz>();

    Array4<Real const> const& fcx = ebdata.get<EBData_t::fcx>();
    Array4<Real const> const& fcy = ebdata.get<EBData_t::fcy>();
    Array4<Real const> const& fcz = ebdata.get<EBData_t::fcz>();

    Array4<Real const> const& ccent = ebdata.get<EBData_t::centroid>();
    Array4<Real const> const& bcent = ebdata.get<EBData_t::bndrycent>();
    Array4<Real const> const& ba = ebdata.get<EBData_t::bndryarea>();

    if (flag(i,j,k).isCovered())
    {
        y(i,j,k,n) = T(0.0);
    }
    else if (flag(i,j,k).isRegular() &&
                ((flag(i-1,j  ,k  ).isRegular() && flag(i+1,j  ,k  ).isRegular() &&
                  flag(i  ,j-1,k  ).isRegular() && flag(i  ,j+1,k  ).isRegular() &&
                  flag(i  ,j  ,k-1).isRegular() && flag(i  ,j  ,k+1).isRegular()) ))
    {
        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
            - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                    -bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
            - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                    -bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
            - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                    -bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
    }
    else
    {
        T kappa = vfrc(i,j,k);
        T apxm = apx(i,j,k);
        T apxp = apx(i+1,j,k);
        T apym = apy(i,j,k);
        T apyp = apy(i,j+1,k);
        T apzm = apz(i,j,k);
        T apzp = apz(i,j,k+1);

        // First get EB-aware slope that doesn't know about extdir
        bool needs_bdry_stencil = (i <= domlo_x) || (i >= domhi_x) ||
                                    (j <= domlo_y) || (j >= domhi_y) ||
                                    (k <= domlo_z) || (k >= domhi_z);

        T fxm = bX(i,j,k,n)*(x(i,j,k,n) - x(i-1,j,k,n));
        if ( (apxm != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i-1,j,k) != T(1.0) || vfrc(i+1,j,k) != T(1.0)) ) {
            T yloc_on_xface = fcx(i,j,k,0);
            T zloc_on_xface = fcx(i,j,k,1);

            if(needs_bdry_stencil) {
                fxm = grad_x_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        yloc_on_xface,zloc_on_xface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fxm = grad_x_of_phi_on_centroids(i,j,k,n,x,phieb,flag,ccent,bcent,
                                                yloc_on_xface,zloc_on_xface,is_eb_dirichlet,is_eb_inhomog);
            }

            fxm *= bX(i,j,k,n);
        }

        T fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i,j,k,n));
        if ( (apxp != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i+1,j,k) != T(1.0) || vfrc(i-1,j,k) != T(1.0)) ) {
            T yloc_on_xface = fcx(i+1,j,k,0);
            T zloc_on_xface = fcx(i+1,j,k,1);

            if(needs_bdry_stencil) {
                fxp = grad_x_of_phi_on_centroids_extdir(i+1,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        yloc_on_xface,zloc_on_xface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fxp = grad_x_of_phi_on_centroids(i+1,j,k,n,x,phieb,flag,ccent,bcent,
                                                yloc_on_xface,zloc_on_xface,is_eb_dirichlet,is_eb_inhomog);
            }

            fxp *= bX(i+1,j,k,n);
        }

        T fym = bY(i,j,k,n)*(x(i,j,k,n) - x(i,j-1,k,n));
        if ( (apym != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j-1,k) != T(1.0) || vfrc(i,j+1,k) != T(1.0)) ) {
            T xloc_on_yface = fcy(i,j,k,0);
            T zloc_on_yface = fcy(i,j,k,1);

            if(needs_bdry_stencil) {
                fym = grad_y_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        xloc_on_yface,zloc_on_yface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fym = grad_y_of_phi_on_centroids(i,j,k,n,x,phieb,flag,ccent,bcent,
                                                xloc_on_yface,zloc_on_yface,is_eb_dirichlet,is_eb_inhomog);
            }

            fym *= bY(i,j,k,n);
        }

        T fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j,k,n));
        if ( (apyp != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j+1,k) != T(1.0) || vfrc(i,j-1,k) != T(1.0)) ) {
            T xloc_on_yface = fcy(i,j+1,k,0);
            T zloc_on_yface = fcy(i,j+1,k,1);

            if(needs_bdry_stencil) {
                fyp = grad_y_of_phi_on_centroids_extdir(i,j+1,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        xloc_on_yface,zloc_on_yface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fyp = grad_y_of_phi_on_centroids(i,j+1,k,n,x,phieb,flag,ccent,bcent,
                                                xloc_on_yface,zloc_on_yface,is_eb_dirichlet,is_eb_inhomog);
            }

            fyp *= bY(i,j+1,k,n);
        }

        T fzm = bZ(i,j,k,n)*(x(i,j,k,n) - x(i,j,k-1,n));
        if ( (apzm != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j,k-1) != T(1.0) || vfrc(i,j,k+1) != T(1.0)) ) {
            T xloc_on_zface = fcz(i,j,k,0);
            T yloc_on_zface = fcz(i,j,k,1);

            if(needs_bdry_stencil) {
                fzm = grad_z_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        xloc_on_zface,yloc_on_zface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fzm = grad_z_of_phi_on_centroids(i,j,k,n,x,phieb,flag,ccent,bcent,
                                                xloc_on_zface,yloc_on_zface,is_eb_dirichlet,is_eb_inhomog);
            }

            fzm *= bZ(i,j,k,n);
        }

        T fzp = bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k,n));
        if ( (apzp != T(0.0)) && (vfrc(i,j,k) != T(1.0) || vfrc(i,j,k+1) != T(1.0)  || vfrc(i,j,k-1) != T(1.0)) ) {
            T xloc_on_zface = fcz(i,j,k+1,0);
            T yloc_on_zface = fcz(i,j,k+1,1);

            if(needs_bdry_stencil) {
                fzp = grad_z_of_phi_on_centroids_extdir(i,j,k+1,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        xloc_on_zface,yloc_on_zface,
                                                        is_eb_dirichlet,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            } else {
                fzp = grad_z_of_phi_on_centroids(i,j,k+1,n,x,phieb,flag,ccent,bcent,
                                                xloc_on_zface,yloc_on_zface,is_eb_dirichlet,is_eb_inhomog);
            }

            fzp *= bZ(i,j,k+1,n);
        }

        T feb = T(0.0);
        if (is_eb_dirichlet && flag(i,j,k).isSingleValued()) {
            T dapx = apxm-apxp;
            T dapy = apym-apyp;
            T dapz = apzm-apzp;
            T anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;
            T anrmz = dapz * anorminv;

            feb = grad_eb_of_phi_on_centroids_extdir(i,j,k,n,x,phieb,flag,ccent,bcent,vfrc,
                                                        anrmx,anrmy,anrmz,is_eb_inhomog,
                                                        on_x_face,domlo_x,domhi_x,
                                                        on_y_face,domlo_y,domhi_y,
                                                        on_z_face,domlo_z,domhi_z);
            feb *= ba(i,j,k) * beb(i,j,k,n);
        }

        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (T(1.0)/kappa) *
            (dhx*(apxm*fxm - apxp*fxp) +
             dhy*(apym*fym - apyp*fyp) +
             dhz*(apzm*fzm - apzp*fzp) - dhx*feb);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx (int i, int j, int k, int n, Array4<T> const& y,
                        Array4<T const> const& x, Array4<T const> const& a,
                        Array4<T const> const& bX, Array4<T const> const& bY,
                        Array4<T const> const& bZ, Array4<const int> const& ccm,
                        EBData const& ebdata,
                        Array4<T const> const& beb,
                        bool is_dirichlet, Array4<T const> const& phieb,
                        bool is_inhomog, T alpha, T beta,
                        GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
                        bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];

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
                    -bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
            - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                    -bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
            - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                    -bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
    }
    else
    {
        T kappa = ebdata.get<EBData_t::volfrac>(i,j,k);
        T apxm = ebdata.get<EBData_t::apx>(i  ,j  ,k  );
        T apxp = ebdata.get<EBData_t::apx>(i+1,j  ,k  );
        T apym = ebdata.get<EBData_t::apy>(i  ,j  ,k  );
        T apyp = ebdata.get<EBData_t::apy>(i  ,j+1,k  );
        T apzm = ebdata.get<EBData_t::apz>(i  ,j  ,k  );
        T apzp = ebdata.get<EBData_t::apz>(i  ,j  ,k+1);

        T fxm = bX(i,j,k,n)*(x(i,j,k,n) - x(i-1,j,k,n));
        if (apxm != T(0.0) && apxm != T(1.0)) {
            auto fcx0 = ebdata.get<EBData_t::fcx>(i,j,k,0);
            auto fcx1 = ebdata.get<EBData_t::fcx>(i,j,k,1);
            int jj = j + static_cast<int>(std::copysign(T(1.0), fcx0));
            int kk = k + static_cast<int>(std::copysign(T(1.0), fcx1));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx0) : T(0.0);
            T fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*fxm +
                    fracy*(T(1.0)-fracz)*bX(i,jj,k ,n)*(x(i,jj,k ,n)-x(i-1,jj,k ,n)) +
                    fracz*(T(1.0)-fracy)*bX(i,j ,kk,n)*(x(i,j ,kk,n)-x(i-1,j ,kk,n)) +
                    fracy*        fracz *bX(i,jj,kk,n)*(x(i,jj,kk,n)-x(i-1,jj,kk,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*(x(i, j, k,n)-x(i-1, j, k,n)) +
                           fracy *(T(1.0)-fracz)*(x(i,jj, k,n)-x(i-1,jj, k,n)) +
                           fracz *(T(1.0)-fracy)*(x(i, j,kk,n)-x(i-1, j,kk,n)) +
                           fracy *        fracz *(x(i,jj,kk,n)-x(i-1,jj,kk,n));
                fxm *= bX(i,j,k,n);
            }
        }

        T fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i,j,k,n));
        if (apxp != T(0.0) && apxp != T(1.0)) {
            auto fcx0 = ebdata.get<EBData_t::fcx>(i+1,j,k,0);
            auto fcx1 = ebdata.get<EBData_t::fcx>(i+1,j,k,1);
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx0));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcx1));
            T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx0) : T(0.0);
            T fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? std::abs(fcx1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fxp = (T(1.0)-fracy)*(T(1.0)-fracz)*fxp +
                    fracy*(T(1.0)-fracz)*bX(i+1,jj,k ,n)*(x(i+1,jj,k ,n)-x(i,jj,k ,n)) +
                    fracz*(T(1.0)-fracy)*bX(i+1,j ,kk,n)*(x(i+1,j ,kk,n)-x(i,j ,kk,n)) +
                    fracy*        fracz *bX(i+1,jj,kk,n)*(x(i+1,jj,kk,n)-x(i,jj,kk,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fxp = (T(1.0)-fracy)*(T(1.0)-fracz)*(x(i+1, j, k,n)-x(i, j, k,n)) +
                           fracy *(T(1.0)-fracz)*(x(i+1,jj, k,n)-x(i,jj, k,n)) +
                           fracz *(T(1.0)-fracy)*(x(i+1, j,kk,n)-x(i, j,kk,n)) +
                           fracy *        fracz *(x(i+1,jj,kk,n)-x(i,jj,kk,n));
                fxp *= bX(i+1,j,k,n);

            }
        }

        T fym = bY(i,j,k,n)*(x(i,j,k,n) - x(i,j-1,k,n));
        if (apym != T(0.0) && apym != T(1.0)) {
            auto fcy0 = ebdata.get<EBData_t::fcy>(i,j,k,0);
            auto fcy1 = ebdata.get<EBData_t::fcy>(i,j,k,1);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy0));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcy1));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy0) : T(0.0);
            T fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fym = (T(1.0)-fracx)*(T(1.0)-fracz)*fym +
                    fracx*(T(1.0)-fracz)*bY(ii,j,k ,n)*(x(ii,j,k ,n)-x(ii,j-1,k ,n)) +
                    fracz*(T(1.0)-fracx)*bY(i ,j,kk,n)*(x(i ,j,kk,n)-x(i ,j-1,kk,n)) +
                    fracx*        fracz *bY(ii,j,kk,n)*(x(ii,j,kk,n)-x(ii,j-1,kk,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fym = (T(1.0)-fracx)*(T(1.0)-fracz)*(x( i,j, k,n)-x( i,j-1, k,n)) +
                           fracx *(T(1.0)-fracz)*(x(ii,j, k,n)-x(ii,j-1, k,n)) +
                           fracz *(T(1.0)-fracx)*(x(i ,j,kk,n)-x( i,j-1,kk,n)) +
                           fracx *        fracz *(x(ii,j,kk,n)-x(ii,j-1,kk,n));
                fym *= bY(i,j,k,n);

            }
        }

        T fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j,k,n));
        if (apyp != T(0.0) && apyp != T(1.0)) {
            auto fcy0 = ebdata.get<EBData_t::fcy>(i,j+1,k,0);
            auto fcy1 = ebdata.get<EBData_t::fcy>(i,j+1,k,1);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy0));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcy1));
            T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy0) : T(0.0);
            T fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? std::abs(fcy1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fyp = (T(1.0)-fracx)*(T(1.0)-fracz)*fyp +
                    fracx*(T(1.0)-fracz)*bY(ii,j+1,k ,n)*(x(ii,j+1,k ,n)-x(ii,j,k ,n)) +
                    fracz*(T(1.0)-fracx)*bY(i ,j+1,kk,n)*(x(i ,j+1,kk,n)-x(i ,j,kk,n)) +
                    fracx*        fracz *bY(ii,j+1,kk,n)*(x(ii,j+1,kk,n)-x(ii,j,kk,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fyp = (T(1.0)-fracx)*(T(1.0)-fracz)*(x( i,j+1, k,n)-x( i,j, k,n)) +
                           fracx *(T(1.0)-fracz)*(x(ii,j+1, k,n)-x(ii,j, k,n)) +
                           fracz *(T(1.0)-fracx)*(x( i,j+1,kk,n)-x( i,j,kk,n)) +
                           fracx *        fracz *(x(ii,j+1,kk,n)-x(ii,j,kk,n));
                fyp *= bY(i,j+1,k,n);

            }
        }

        T fzm = bZ(i,j,k,n)*(x(i,j,k,n) - x(i,j,k-1,n));
        if (apzm != T(0.0) && apzm != T(1.0)) {
            auto fcz0 = ebdata.get<EBData_t::fcz>(i,j,k,0);
            auto fcz1 = ebdata.get<EBData_t::fcz>(i,j,k,1);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcz0));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcz1));
            T fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz0) : T(0.0);
            T fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*fzm +
                    fracx*(T(1.0)-fracy)*bZ(ii,j ,k,n)*(x(ii,j ,k,n)-x(ii,j ,k-1,n)) +
                    fracy*(T(1.0)-fracx)*bZ(i ,jj,k,n)*(x(i ,jj,k,n)-x(i ,jj,k-1,n)) +
                    fracx*        fracy *bZ(ii,jj,k,n)*(x(ii,jj,k,n)-x(ii,jj,k-1,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*(x( i, j,k,n)-x( i, j,k-1,n)) +
                           fracx *(T(1.0)-fracy)*(x(ii, j,k,n)-x(ii, j,k-1,n)) +
                           fracy *(T(1.0)-fracx)*(x( i,jj,k,n)-x( i,jj,k-1,n)) +
                           fracx *        fracy *(x(ii,jj,k,n)-x(ii,jj,k-1,n));
                fzm *= bZ(i,j,k,n);

            }
        }

        T fzp = bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k,n));
        if (apzp != T(0.0) && apzp != T(1.0)) {
            auto fcz0 = ebdata.get<EBData_t::fcz>(i,j,k+1,0);
            auto fcz1 = ebdata.get<EBData_t::fcz>(i,j,k+1,1);
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcz0));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcz1));
            T fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? std::abs(fcz0) : T(0.0);
            T fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? std::abs(fcz1) : T(0.0);
            if (beta_on_center && phi_on_center)
            {
                fzp = (T(1.0)-fracx)*(T(1.0)-fracy)*fzp +
                    fracx*(T(1.0)-fracy)*bZ(ii,j ,k+1,n)*(x(ii,j ,k+1,n)-x(ii,j ,k,n)) +
                    fracy*(T(1.0)-fracx)*bZ(i ,jj,k+1,n)*(x(i ,jj,k+1,n)-x(i ,jj,k,n)) +
                    fracx*        fracy *bZ(ii,jj,k+1,n)*(x(ii,jj,k+1,n)-x(ii,jj,k,n));
            }
            else if (beta_on_centroid && phi_on_center)
            {
                fzp = (T(1.0)-fracx)*(T(1.0)-fracy)*(x( i, j,k+1,n)-x( i, j,k,n)) +
                           fracx *(T(1.0)-fracy)*(x(ii, j,k+1,n)-x(ii, j,k,n)) +
                           fracy *(T(1.0)-fracx)*(x( i,jj,k+1,n)-x( i,jj,k,n)) +
                           fracx *        fracy *(x(ii,jj,k+1,n)-x(ii,jj,k,n));
                fzp *= bZ(i,j,k+1,n);

            }
        }

        T feb = T(0.0);
        if (is_dirichlet) {
            T dapx = apxm-apxp;
            T dapy = apym-apyp;
            T dapz = apzm-apzp;
            T anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;
            T anrmz = dapz * anorminv;

            T phib = is_inhomog ? phieb(i,j,k,n) : T(0.0);

            T bctx = ebdata.get<EBData_t::bndrycent>(i,j,k,0);
            T bcty = ebdata.get<EBData_t::bndrycent>(i,j,k,1);
            T bctz = ebdata.get<EBData_t::bndrycent>(i,j,k,2);
            T dx_eb = get_dx_eb(kappa);

            T dg = dx_eb / amrex::max(std::abs(anrmx), std::abs(anrmy),
                                         std::abs(anrmz));
            T gx = bctx - dg*anrmx;
            T gy = bcty - dg*anrmy;
            T gz = bctz - dg*anrmz;
            T sx = std::copysign(T(1.0),anrmx);
            T sy = std::copysign(T(1.0),anrmy);
            T sz = std::copysign(T(1.0),anrmz);
            int ii = i - static_cast<int>(sx);
            int jj = j - static_cast<int>(sy);
            int kk = k - static_cast<int>(sz);

            gx = sx*gx;
            gy = sy*gy;
            gz = sz*gz;
            T gxy = gx*gy;
            T gxz = gx*gz;
            T gyz = gy*gz;
            T gxyz = gx*gy*gz;
            T phig = (T(1.0)+gx+gy+gz+gxy+gxz+gyz+gxyz) * x(i ,j ,k ,n)
                +       (-gz - gxz - gyz - gxyz)        * x(i ,j ,kk,n)
                +       (-gy - gxy - gyz - gxyz)        * x(i ,jj,k ,n)
                +       (gyz + gxyz)                    * x(i ,jj,kk,n)
                +       (-gx - gxy - gxz - gxyz)        * x(ii,j ,k ,n)
                +       (gxz + gxyz)                    * x(ii,j ,kk,n)
                +       (gxy + gxyz)                    * x(ii,jj,k ,n)
                +       (-gxyz)                         * x(ii,jj,kk,n);

            T dphidn = (phib-phig)/dg;

            T ba = ebdata.get<EBData_t::bndryarea>(i,j,k);
            feb = dphidn * ba * beb(i,j,k,n);
        }

        y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (T(1.0)/kappa) *
            (dhx*(apxm*fxm - apxp*fxp) +
             dhy*(apym*fym - apyp*fyp) +
             dhz*(apzm*fzm - apzp*fzp) - dhx*feb);
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
                         Array4<T const> const& apz,
                         Array4<T const> const& bc,
                         Array4<T const> const& beb,
                         Array4<T const> const& phieb,
                         bool is_inhomog,
                         GpuArray<T,AMREX_SPACEDIM> const& dxinv) noexcept
{
    T dhx = dxinv[0];

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
        T apzm = apz(i,j,k);
        T apzp = apz(i,j,k+1);

        T dapx = apxm-apxp;
        T dapy = apym-apyp;
        T dapz = apzm-apzp;
        T anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
        T anorminv = T(1.0)/anorm;
        T anrmx = dapx * anorminv;
        T anrmy = dapy * anorminv;
        T anrmz = dapz * anorminv;

        T phib = is_inhomog ? phieb(i,j,k,n) : T(0.0);

        T bctx = bc(i,j,k,0);
        T bcty = bc(i,j,k,1);
        T bctz = bc(i,j,k,2);
        T dx_eb = get_dx_eb(kappa);

        T dg = dx_eb / amrex::max(std::abs(anrmx), std::abs(anrmy), std::abs(anrmz));
        T gx = bctx - dg*anrmx;
        T gy = bcty - dg*anrmy;
        T gz = bctz - dg*anrmz;
        T sx = std::copysign(T(1.0),anrmx);
        T sy = std::copysign(T(1.0),anrmy);
        T sz = std::copysign(T(1.0),anrmz);
        int ii = i - static_cast<int>(sx);
        int jj = j - static_cast<int>(sy);
        int kk = k - static_cast<int>(sz);

        gx = sx*gx;
        gy = sy*gy;
        gz = sz*gz;
        T gxy = gx*gy;
        T gxz = gx*gz;
        T gyz = gy*gz;
        T gxyz = gx*gy*gz;
        T phig = (T(1.0)+gx+gy+gz+gxy+gxz+gyz+gxyz) * x(i ,j ,k ,n)
            +       (-gz - gxz - gyz - gxyz)        * x(i ,j ,kk,n)
            +       (-gy - gxy - gyz - gxyz)        * x(i ,jj,k ,n)
            +       (gyz + gxyz)                    * x(i ,jj,kk,n)
            +       (-gx - gxy - gxz - gxyz)        * x(ii,j ,k ,n)
            +       (gxz + gxyz)                    * x(ii,j ,kk,n)
            +       (gxy + gxyz)                    * x(ii,jj,k ,n)
            +       (-gxyz)                         * x(ii,jj,kk,n);

        T dphidn = dhx*(phib-phig)/dg;
        feb(i,j,k,n) = -beb(i,j,k,n) * dphidn;
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_gsrb (int i, int j, int k, int n,
                       Array4<T> const& phi, Array4<T const> const& rhs,
                       T alpha, Array4<T const> const& a,
                       T beta, GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                       Array4<T const> const& bX, Array4<T const> const& bY,
                       Array4<T const> const& bZ,
                       Array4<int const> const& m0, Array4<int const> const& m2,
                       Array4<int const> const& m4,
                       Array4<int const> const& m1, Array4<int const> const& m3,
                       Array4<int const> const& m5,
                       Array4<T const> const& f0, Array4<T const> const& f2,
                       Array4<T const> const& f4,
                       Array4<T const> const& f1, Array4<T const> const& f3,
                       Array4<T const> const& f5,
                       Array4<const int> const& ccm, Array4<T const> const& beb,
                       EBData const& ebdata,
                       bool is_dirichlet, bool beta_on_centroid, bool phi_on_centroid,
                       Box const& vbox, int redblack) noexcept
{
    constexpr T omega = 1.15;
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    const T dhx = beta*dxinv[0]*dxinv[0];
    const T dhy = beta*dxinv[1]*dxinv[1];
    const T dhz = beta*dxinv[2]*dxinv[2];

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
            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                ? f2(i,j,vlo.z,n) : T(0.0);
            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                ? f3(vhi.x,j,k,n) : T(0.0);
            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                ? f4(i,vhi.y,k,n) : T(0.0);
            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                ? f5(i,j,vhi.z,n) : T(0.0);

            if (flag.isRegular())
            {
                T gamma = alpha*a(i,j,k)
                    + dhx*(bX(i+1,j,k,n) + bX(i,j,k,n))
                    + dhy*(bY(i,j+1,k,n) + bY(i,j,k,n))
                    + dhz*(bZ(i,j,k+1,n) + bZ(i,j,k,n));

                T rho = dhx*(bX(i+1,j  ,k  ,n)*phi(i+1,j  ,k  ,n) +
                             bX(i  ,j  ,k  ,n)*phi(i-1,j  ,k  ,n))
                    +   dhy*(bY(i  ,j+1,k  ,n)*phi(i  ,j+1,k  ,n) +
                             bY(i  ,j  ,k  ,n)*phi(i  ,j-1,k  ,n))
                    +   dhz*(bZ(i  ,j  ,k+1,n)*phi(i  ,j  ,k+1,n) +
                             bZ(i  ,j  ,k  ,n)*phi(i  ,j  ,k-1,n));

                T delta = dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                    +     dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                    +     dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5);

                T res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                phi(i,j,k,n) += omega*res/(gamma-delta);
            }
            else
            {
                T kappa = ebdata.get<EBData_t::volfrac>(i,j,k);
                T apxm = ebdata.get<EBData_t::apx>(i  ,j  ,k  );
                T apxp = ebdata.get<EBData_t::apx>(i+1,j  ,k  );
                T apym = ebdata.get<EBData_t::apy>(i  ,j  ,k  );
                T apyp = ebdata.get<EBData_t::apy>(i  ,j+1,k  );
                T apzm = ebdata.get<EBData_t::apz>(i  ,j  ,k  );
                T apzp = ebdata.get<EBData_t::apz>(i  ,j  ,k+1);

                T fxm = -bX(i,j,k,n)*phi(i-1,j,k,n);
                T oxm = -bX(i,j,k,n)*cf0;
                T sxm =  bX(i,j,k,n);
                if (apxm != T(0.0) && apxm != T(1.0)) {
                    auto fcx0 = ebdata.get<EBData_t::fcx>(i,j,k,0);
                    auto fcx1 = ebdata.get<EBData_t::fcx>(i,j,k,1);
                    int jj = j + static_cast<int>(std::copysign(T(1.0), fcx0));
                    int kk = k + static_cast<int>(std::copysign(T(1.0), fcx1));
                    T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                        ? std::abs(fcx0) : T(0.0);
                    T fracz = (ccm(i-1,j,kk) || ccm(i,j,kk))
                        ? std::abs(fcx1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*fxm
                                +     fracy *(T(1.0)-fracz)*bX(i,jj,k ,n)*(phi(i,jj,k ,n)-phi(i-1,jj,k ,n))
                                +(T(1.0)-fracy)*     fracz *bX(i,j ,kk,n)*(phi(i,j ,kk,n)-phi(i-1,j ,kk,n))
                                +        fracy *     fracz *bX(i,jj,kk,n)*(phi(i,jj,kk,n)-phi(i-1,jj,kk,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*(              -phi(i-1, j, k,n))
                                +     fracy *(T(1.0)-fracz)*(phi(i,jj,k ,n)-phi(i-1,jj, k,n))
                                +(T(1.0)-fracy)*     fracz *(phi(i,j ,kk,n)-phi(i-1, j,kk,n))
                                +        fracy *     fracz *(phi(i,jj,kk,n)-phi(i-1,jj,kk,n));
                        fxm *= bX(i,j,k,n);

                    }
                    oxm = T(0.0);
                    sxm = (T(1.0)-fracy)*(T(1.0)-fracz)*sxm;
                }

                T fxp =  bX(i+1,j,k,n)*phi(i+1,j,k,n);
                T oxp =  bX(i+1,j,k,n)*cf3;
                T sxp = -bX(i+1,j,k,n);
                if (apxp != T(0.0) && apxp != T(1.0)) {
                    auto fcx0 = ebdata.get<EBData_t::fcx>(i+1,j,k,0);
                    auto fcx1 = ebdata.get<EBData_t::fcx>(i+1,j,k,1);
                    int jj = j + static_cast<int>(std::copysign(T(1.0),fcx0));
                    int kk = k + static_cast<int>(std::copysign(T(1.0),fcx1));
                    T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                        ? std::abs(fcx0) : T(0.0);
                    T fracz = (ccm(i,j,kk) || ccm(i+1,j,kk))
                        ? std::abs(fcx1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fxp = (T(1.0)-fracy)*(T(1.0)-fracz)*fxp
                                +    fracy *(T(1.0)-fracz)*bX(i+1,jj,k ,n)*(phi(i+1,jj,k ,n)-phi(i,jj,k ,n))
                                +(T(1.0)-fracy)*     fracz *bX(i+1,j ,kk,n)*(phi(i+1,j ,kk,n)-phi(i,j ,kk,n))
                                +        fracy *     fracz *bX(i+1,jj,kk,n)*(phi(i+1,jj,kk,n)-phi(i,jj,kk,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fxp = (T(1.0)-fracy)*(T(1.0)-fracz)*(phi(i+1, j, k,n)               ) +
                                    fracy *(T(1.0)-fracz)*(phi(i+1,jj, k,n)-phi(i,jj, k,n)) +
                                    fracz *(T(1.0)-fracy)*(phi(i+1, j,kk,n)-phi(i, j,kk,n)) +
                                    fracy *        fracz *(phi(i+1,jj,kk,n)-phi(i,jj,kk,n));
                        fxp *= bX(i+1,j,k,n);
                    }
                    oxp = T(0.0);
                    sxp = (T(1.0)-fracy)*(T(1.0)-fracz)*sxp;
                }

                T fym = -bY(i,j,k,n)*phi(i,j-1,k,n);
                T oym = -bY(i,j,k,n)*cf1;
                T sym =  bY(i,j,k,n);
                if (apym != T(0.0) && apym != T(1.0)) {
                    auto fcy0 = ebdata.get<EBData_t::fcy>(i,j,k,0);
                    auto fcy1 = ebdata.get<EBData_t::fcy>(i,j,k,1);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcy0));
                    int kk = k + static_cast<int>(std::copysign(T(1.0),fcy1));
                    T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                        ? std::abs(fcy0) : T(0.0);
                    T fracz = (ccm(i,j-1,kk) || ccm(i,j,kk))
                        ? std::abs(fcy1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fym = (T(1.0)-fracx)*(T(1.0)-fracz)*fym
                            +      fracx *(T(1.0)-fracz)*bY(ii,j,k ,n)*(phi(ii,j,k ,n)-phi(ii,j-1,k ,n))
                            + (T(1.0)-fracx)*     fracz *bY(i ,j,kk,n)*(phi(i ,j,kk,n)-phi(i ,j-1,kk,n))
                            +      fracx *     fracz *bY(ii,j,kk,n)*(phi(ii,j,kk,n)-phi(ii,j-1,kk,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fym = (T(1.0)-fracx)*(T(1.0)-fracz)*(              -phi( i,j-1, k,n))
                            +      fracx *(T(1.0)-fracz)*(phi(ii,j,k ,n)-phi(ii,j-1, k,n))
                            + (T(1.0)-fracx)*     fracz *(phi(i ,j,kk,n)-phi( i,j-1,kk,n))
                            +      fracx *        fracz *(phi(ii,j,kk,n)-phi(ii,j-1,kk,n));
                        fym *= bY(i,j,k,n);
                    }
                    oym = T(0.0);
                    sym = (T(1.0)-fracx)*(T(1.0)-fracz)*sym;
                }

                T fyp =  bY(i,j+1,k,n)*phi(i,j+1,k,n);
                T oyp =  bY(i,j+1,k,n)*cf4;
                T syp = -bY(i,j+1,k,n);
                if (apyp != T(0.0) && apyp != T(1.0)) {
                    auto fcy0 = ebdata.get<EBData_t::fcy>(i,j+1,k,0);
                    auto fcy1 = ebdata.get<EBData_t::fcy>(i,j+1,k,1);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcy0));
                    int kk = k + static_cast<int>(std::copysign(T(1.0),fcy1));
                    T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                        ? std::abs(fcy0) : T(0.0);
                    T fracz = (ccm(i,j,kk) || ccm(i,j+1,kk))
                        ? std::abs(fcy1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fyp = (T(1.0)-fracx)*(T(1.0)-fracz)*fyp
                            +      fracx *(T(1.0)-fracz)*bY(ii,j+1,k ,n)*(phi(ii,j+1,k ,n)-phi(ii,j,k ,n))
                            + (T(1.0)-fracx)*     fracz *bY(i ,j+1,kk,n)*(phi(i ,j+1,kk,n)-phi(i ,j,kk,n))
                            +      fracx *        fracz *bY(ii,j+1,kk,n)*(phi(ii,j+1,kk,n)-phi(ii,j,kk,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fyp = (T(1.0)-fracx)*(T(1.0)-fracz)*(phi( i,j+1, k,n)               )
                            +      fracx *(T(1.0)-fracz)*(phi(ii,j+1, k,n)-phi(ii,j, k,n))
                            + (T(1.0)-fracx)*     fracz *(phi( i,j+1,kk,n)-phi( i,j,kk,n))
                            +      fracx *        fracz *(phi(ii,j+1,kk,n)-phi(ii,j,kk,n));
                        fyp *= bY(i,j+1,k,n);
                    }
                    oyp = T(0.0);
                    syp = (T(1.0)-fracx)*(T(1.0)-fracz)*syp;
                }

                T fzm = -bZ(i,j,k,n)*phi(i,j,k-1,n);
                T ozm = -bZ(i,j,k,n)*cf2;
                T szm =  bZ(i,j,k,n);
                if (apzm != T(0.0) && apzm != T(1.0)) {
                    auto fcz0 = ebdata.get<EBData_t::fcz>(i,j,k,0);
                    auto fcz1 = ebdata.get<EBData_t::fcz>(i,j,k,1);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcz0));
                    int jj = j + static_cast<int>(std::copysign(T(1.0),fcz1));
                    T fracx = (ccm(ii,j,k-1) || ccm(ii,j,k))
                        ? std::abs(fcz0) : T(0.0);
                    T fracy = (ccm(i,jj,k-1) || ccm(i,jj,k))
                        ? std::abs(fcz1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*fzm
                                +     fracx *(T(1.0)-fracy)*bZ(ii, j,k,n)*(phi(ii, j,k,n)-phi(ii, j,k-1,n))
                                +(T(1.0)-fracx)*     fracy *bZ( i,jj,k,n)*(phi( i,jj,k,n)-phi( i,jj,k-1,n))
                                +     fracx *        fracy *bZ(ii,jj,k,n)*(phi(ii,jj,k,n)-phi(ii,jj,k-1,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*(              -phi( i, j,k-1,n))
                            +      fracx *(T(1.0)-fracy)*(phi(ii, j,k,n)-phi(ii, j,k-1,n))
                            + (T(1.0)-fracx)*     fracy *(phi( i,jj,k,n)-phi(i ,jj,k-1,n))
                            +      fracx *        fracy *(phi(ii,jj,k,n)-phi(ii,jj,k-1,n));
                        fzm *= bZ(i,j,k,n);

                    }
                    ozm = T(0.0);
                    szm = (T(1.0)-fracx)*(T(1.0)-fracy)*szm;
                }

                T fzp =  bZ(i,j,k+1,n)*phi(i,j,k+1,n);
                T ozp =  bZ(i,j,k+1,n)*cf5;
                T szp = -bZ(i,j,k+1,n);
                if (apzp != T(0.0) && apzp != T(1.0)) {
                    auto fcz0 = ebdata.get<EBData_t::fcz>(i,j,k+1,0);
                    auto fcz1 = ebdata.get<EBData_t::fcz>(i,j,k+1,1);
                    int ii = i + static_cast<int>(std::copysign(T(1.0),fcz0));
                    int jj = j + static_cast<int>(std::copysign(T(1.0),fcz1));
                    T fracx = (ccm(ii,j,k) || ccm(ii,j,k+1))
                        ? std::abs(fcz0) : T(0.0);
                    T fracy = (ccm(i,jj,k) || ccm(i,jj,k+1))
                        ? std::abs(fcz1) : T(0.0);
                    if (!beta_on_centroid && !phi_on_centroid)
                    {
                        fzp = (T(1.0)-fracx)*(T(1.0)-fracy)*fzp
                            +      fracx *(T(1.0)-fracy)*bZ(ii,j ,k+1,n)*(phi(ii,j ,k+1,n)-phi(ii,j ,k,n))
                            + (T(1.0)-fracx)*     fracy *bZ(i ,jj,k+1,n)*(phi(i ,jj,k+1,n)-phi(i ,jj,k,n))
                            +      fracx *        fracy *bZ(ii,jj,k+1,n)*(phi(ii,jj,k+1,n)-phi(ii,jj,k,n));
                    }
                    else if (beta_on_centroid && !phi_on_centroid)
                    {
                        fzp = (T(1.0)-fracx)*(T(1.0)-fracy)*(phi( i, j,k+1,n)               )
                            +      fracx *(T(1.0)-fracy)*(phi(ii, j,k+1,n)-phi(ii, j,k,n))
                            + (T(1.0)-fracx)*     fracy *(phi( i,jj,k+1,n)-phi( i,jj,k,n))
                            +      fracx *        fracy *(phi(ii,jj,k+1,n)-phi(ii,jj,k,n));
                        fzp *= bZ(i,j,k+1,n);

                    }
                    ozp = T(0.0);
                    szp = (T(1.0)-fracx)*(T(1.0)-fracy)*szp;
                }

                T vfrcinv = T(1.0)/kappa;
                T gamma = alpha*a(i,j,k) + vfrcinv *
                    (dhx*(apxm*sxm-apxp*sxp) +
                     dhy*(apym*sym-apyp*syp) +
                     dhz*(apzm*szm-apzp*szp));

                T rho = -vfrcinv *
                    (dhx*(apxm*fxm-apxp*fxp) +
                     dhy*(apym*fym-apyp*fyp) +
                     dhz*(apzm*fzm-apzp*fzp));

                T delta = -vfrcinv *
                    (dhx*(apxm*oxm-apxp*oxp) +
                     dhy*(apym*oym-apyp*oyp) +
                     dhz*(apzm*ozm-apzp*ozp));

                if (is_dirichlet) {
                    T dapx = apxm-apxp;
                    T dapy = apym-apyp;
                    T dapz = apzm-apzp;
                    T anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
                    T anorminv = T(1.0)/anorm;
                    T anrmx = dapx * anorminv;
                    T anrmy = dapy * anorminv;
                    T anrmz = dapz * anorminv;
                    T bctx = ebdata.get<EBData_t::bndrycent>(i,j,k,0);
                    T bcty = ebdata.get<EBData_t::bndrycent>(i,j,k,1);
                    T bctz = ebdata.get<EBData_t::bndrycent>(i,j,k,2);
                    T dx_eb = get_dx_eb(kappa);

                    T dg = dx_eb / amrex::max(std::abs(anrmx),std::abs(anrmy),
                                                 std::abs(anrmz));

                    T gx = bctx - dg*anrmx;
                    T gy = bcty - dg*anrmy;
                    T gz = bctz - dg*anrmz;
                    T sx = std::copysign(T(1.0),anrmx);
                    T sy = std::copysign(T(1.0),anrmy);
                    T sz = std::copysign(T(1.0),anrmz);
                    int ii = i - static_cast<int>(sx);
                    int jj = j - static_cast<int>(sy);
                    int kk = k - static_cast<int>(sz);

                    gx *= sx;
                    gy *= sy;
                    gz *= sz;
                    T gxy = gx*gy;
                    T gxz = gx*gz;
                    T gyz = gy*gz;
                    T gxyz = gx*gy*gz;
                    T phig_gamma = (T(1.0)+gx+gy+gz+gxy+gxz+gyz+gxyz);
                    T phig =    (-gz - gxz - gyz - gxyz)        * phi(i ,j ,kk,n)
                        +       (-gy - gxy - gyz - gxyz)        * phi(i ,jj,k ,n)
                        +       (gyz + gxyz)                    * phi(i ,jj,kk,n)
                        +       (-gx - gxy - gxz - gxyz)        * phi(ii,j ,k ,n)
                        +       (gxz + gxyz)                    * phi(ii,j ,kk,n)
                        +       (gxy + gxyz)                    * phi(ii,jj,k ,n)
                        +       (-gxyz)                         * phi(ii,jj,kk,n);

                    T ba = ebdata.get<EBData_t::bndryarea>(i,j,k);

                    T dphidn    = (    -phig)/dg;
                    T feb_gamma = -phig_gamma/dg * ba * beb(i,j,k,n);
                    gamma += vfrcinv*(-dhx)*feb_gamma;
                    T feb = dphidn * ba * beb(i,j,k,n);
                    rho += -vfrcinv*(-dhx)*feb;
                }

                T res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                phi(i,j,k,n) += omega*res/(gamma-delta);
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
                T fxm = bX(i,j,k,n)*(sol(i,j,k,n) - sol(i-1,j,k,n));
                int jj = j + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,1)));
                T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : T(0.0);
                T fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : T(0.0);
                if (!beta_on_centroid && !phi_on_centroid)
                {
                    fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*fxm +
                               fracy *(T(1.0)-fracz)*bX(i,jj,k ,n)*(sol(i,jj,k ,n)-sol(i-1,jj,k ,n)) +
                               fracz *(T(1.0)-fracy)*bX(i,j ,kk,n)*(sol(i,j ,kk,n)-sol(i-1,j ,kk,n)) +
                               fracy*      fracz *bX(i,jj,kk,n)*(sol(i,jj,kk,n)-sol(i-1,jj,kk,n));
                }
                else if (beta_on_centroid && !phi_on_centroid)
                {
                    fxm = (T(1.0)-fracy)*(T(1.0)-fracz)*(sol(i, j, k,n)-sol(i-1, j, k,n)) +
                               fracy *(T(1.0)-fracz)*(sol(i,jj, k,n)-sol(i-1,jj, k,n)) +
                               fracz *(T(1.0)-fracy)*(sol(i, j,kk,n)-sol(i-1, j,kk,n)) +
                               fracy*      fracz *(sol(i,jj,kk,n)-sol(i-1,jj,kk,n));
                    fxm *= bX(i,j,k,n);

                }
                fx(i,j,k,n) = -dhx*fxm;
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
                T fym = bY(i,j,k,n)*(sol(i,j,k,n) - sol(i,j-1,k,n));
                int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,1)));
                T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : T(0.0);
                T fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : T(0.0);
                if (!beta_on_centroid && !phi_on_centroid)
                {
                    fym = (T(1.0)-fracx)*(T(1.0)-fracz)*fym +
                               fracx *(T(1.0)-fracz)*bY(ii,j,k ,n)*(sol(ii,j,k ,n)-sol(ii,j-1,k ,n)) +
                               fracz *(T(1.0)-fracx)*bY(i ,j,kk,n)*(sol(i ,j,kk,n)-sol(i ,j-1,kk,n)) +
                               fracx *     fracz *bY(ii,j,kk,n)*(sol(ii,j,kk,n)-sol(ii,j-1,kk,n));
                }
                else if (beta_on_centroid && !phi_on_centroid)
                {
                    fym = (T(1.0)-fracx)*(T(1.0)-fracz)*(sol( i,j, k,n)-sol( i,j-1, k,n)) +
                               fracx *(T(1.0)-fracz)*(sol(ii,j, k,n)-sol(ii,j-1, k,n)) +
                               fracz *(T(1.0)-fracx)*(sol( i,j,kk,n)-sol( i,j-1,kk,n)) +
                               fracx *     fracz *(sol(ii,j,kk,n)-sol(ii,j-1,kk,n));
                    fym *= bY(i,j,k,n);
                }

                fy(i,j,k,n) = -dhy*fym;
            }
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_z (Box const& box, Array4<T> const& fz, Array4<T const> const& apz,
                         Array4<T const> const& fcz, Array4<T const> const& sol,
                         Array4<T const> const& bZ, Array4<int const> const& ccm,
                         T dhz, int face_only, int ncomp, Box const& zbox,
                         bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    int lof = zbox.smallEnd(2);
    int hif = zbox.bigEnd(2);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == k || hif == k) {
            if (apz(i,j,k) == T(0.0)) {
                fz(i,j,k,n) = T(0.0);
            } else if (apz(i,j,k) == T(1.0)) {
                fz(i,j,k,n) = -dhz*bZ(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
            } else {
                T fzm = bZ(i,j,k,n)*(sol(i,j,k,n) - sol(i,j,k-1,n));
                int ii = i + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,0)));
                int jj = j + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,1)));
                T fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : T(0.0);
                T fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : T(0.0);
                if (!beta_on_centroid && !phi_on_centroid)
                {
                    fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*fzm +
                        fracx*(T(1.0)-fracy)*bZ(ii,j ,k,n)*(sol(ii,j ,k,n)-sol(ii,j ,k-1,n)) +
                        fracy*(T(1.0)-fracx)*bZ(i ,jj,k,n)*(sol(i ,jj,k,n)-sol(i ,jj,k-1,n)) +
                        fracx*        fracy *bZ(ii,jj,k,n)*(sol(ii,jj,k,n)-sol(ii,jj,k-1,n));
                }
                else if (beta_on_centroid && !phi_on_centroid)
                {
                    fzm = (T(1.0)-fracx)*(T(1.0)-fracy)*(sol( i, j,k,n)-sol( i, j,k-1,n)) +
                                  fracx *(T(1.0)-fracy)*(sol(ii, j,k,n)-sol(ii, j,k-1,n)) +
                                  fracy *(T(1.0)-fracx)*(sol( i,jj,k,n)-sol( i,jj,k-1,n)) +
                                  fracx *        fracy *(sol(ii,jj,k,n)-sol(ii,jj,k-1,n));
                    fzm *= bZ(i,j,k,n);
                }

                fz(i,j,k,n) = -dhz*fzm;
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
void mlebabeclap_flux_z_0 (Box const& box, Array4<T> const& fz, Array4<T const> const& apz,
                           Array4<T const> const& sol, Array4<T const> const& bZ,
                           T dhz, int face_only, int ncomp, Box const& zbox) noexcept
{
    int lof = zbox.smallEnd(2);
    int hif = zbox.bigEnd(2);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only || lof == k || hif == k) {
            if (apz(i,j,k) == T(0.0)) {
                fz(i,j,k,n) = T(0.0);
            } else {
                fz(i,j,k,n) = -dhz*bZ(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
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
            T gxm = (sol(i,j,k,n) - sol(i-1,j,k,n));
            int jj = j + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,1)));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : T(0.0);
            T fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : T(0.0);
            if (!phi_on_centroid)
            {
                gxm = (T(1.0)-fracy)*(T(1.0)-fracz)*gxm +
                    fracy*(T(1.0)-fracz)*(sol(i,jj,k ,n)-sol(i-1,jj,k ,n)) +
                    fracz*(T(1.0)-fracy)*(sol(i,j ,kk,n)-sol(i-1,j ,kk,n)) +
                    fracy*        fracz *(sol(i,jj,kk,n)-sol(i-1,jj,kk,n));
            }
            gx(i,j,k,n) = dxi*gxm;
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
            T gym = (sol(i,j,k,n) - sol(i,j-1,k,n));
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,1)));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : T(0.0);
            T fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : T(0.0);
            if (!phi_on_centroid)
            {
                gym = (T(1.0)-fracx)*(T(1.0)-fracz)*gym +
                    fracx*(T(1.0)-fracz)*(sol(ii,j,k ,n)-sol(ii,j-1,k ,n)) +
                    fracz*(T(1.0)-fracx)*(sol(i ,j,kk,n)-sol(i ,j-1,kk,n)) +
                    fracx*        fracz *(sol(ii,j,kk,n)-sol(ii,j-1,kk,n));
            }
            gy(i,j,k,n) = dyi*gym;
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_z (Box const& box, Array4<T> const& gz, Array4<T const> const& sol,
                         Array4<T const> const& apz, Array4<T const> const& fcz,
                         Array4<int const> const& ccm,
                         T dzi, int ncomp, bool phi_on_centroid) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) == T(0.0)) {
            gz(i,j,k,n) = T(0.0);
        } else if (apz(i,j,k) == T(1.0)) {
            gz(i,j,k,n) = dzi*(sol(i,j,k,n)-sol(i,j,k-1,n));
        } else {
            T gzm = (sol(i,j,k,n) - sol(i,j,k-1,n));
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,0)));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,1)));
            T fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : T(0.0);
            T fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : T(0.0);
            if (!phi_on_centroid)
            {
                gzm = (T(1.0)-fracx)*(T(1.0)-fracy)*gzm +
                    fracx*(T(1.0)-fracy)*(sol(ii,j ,k,n)-sol(ii,j ,k-1,n)) +
                    fracy*(T(1.0)-fracx)*(sol(i ,jj,k,n)-sol(i ,jj,k-1,n)) +
                    fracx*     fracy *(sol(ii,jj,k,n)-sol(ii,jj,k-1,n));
            }
            gz(i,j,k,n) = dzi*gzm;
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
void mlebabeclap_grad_z_0 (Box const& box, Array4<T> const& gz, Array4<T const> const& sol,
                           Array4<T const> const& apz, T dzi, int ncomp) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) == T(0.0)) {
            gz(i,j,k,n) = T(0.0);
        } else {
            gz(i,j,k,n) = dzi*(sol(i,j,k,n)-sol(i,j,k-1,n));
        }
    });
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_normalize (int i, int j, int k, int n, Array4<T> const& phi,
                            T alpha, Array4<T const> const& a,
                            T dhx, T dhy, T dhz,
                            Array4<T const> const& bX, Array4<T const> const& bY,
                            Array4<T const> const& bZ,
                            Array4<const int> const& ccm, Array4<EBCellFlag const> const& flag,
                            Array4<T const> const& vfrc,
                            Array4<T const> const& apx, Array4<T const> const& apy,
                            Array4<T const> const& apz,
                            Array4<T const> const& fcx, Array4<T const> const& fcy,
                            Array4<T const> const& fcz,
                            Array4<T const> const& ba, Array4<T const> const& bc,
                            Array4<T const> const& beb,
                            bool is_dirichlet, bool beta_on_centroid) noexcept
{
    if (flag(i,j,k).isRegular())
    {
        phi(i,j,k,n) /= alpha*a(i,j,k) + dhx*(bX(i,j,k,n) + bX(i+1,j,k,n))
                                       + dhy*(bY(i,j,k,n) + bY(i,j+1,k,n))
                                       + dhz*(bZ(i,j,k,n) + bZ(i,j,k+1,n));
    }
    else if (flag(i,j,k).isSingleValued())
    {
        T kappa = vfrc(i,j,k);
        T apxm = apx(i,j,k);
        T apxp = apx(i+1,j,k);
        T apym = apy(i,j,k);
        T apyp = apy(i,j+1,k);
        T apzm = apz(i,j,k);
        T apzp = apz(i,j,k+1);

        T sxm =  bX(i,j,k,n);
        if (apxm != T(0.0) && apxm != T(1.0) && !beta_on_centroid) {
            int jj = j + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0), fcx(i,j,k,1)));
            T fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                ? std::abs(fcx(i,j,k,0)) : T(0.0);
            T fracz = (ccm(i-1,j,kk) || ccm(i,j,kk))
                ? std::abs(fcx(i,j,k,1)) : T(0.0);
            sxm = (T(1.0)-fracy)*(T(1.0)-fracz)*sxm;
        }

        T sxp = -bX(i+1,j,k,n);
        if (apxp != T(0.0) && apxp != T(1.0) && !beta_on_centroid) {
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcx(i+1,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcx(i+1,j,k,1)));
            T fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                ? std::abs(fcx(i+1,j,k,0)) : T(0.0);
            T fracz = (ccm(i,j,kk) || ccm(i+1,j,kk))
                ? std::abs(fcx(i+1,j,k,1)) : T(0.0);
            sxp = (T(1.0)-fracy)*(T(1.0)-fracz)*sxp;
        }

        T sym =  bY(i,j,k,n);
        if (apym != T(0.0) && apym != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcy(i,j,k,1)));
            T fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                ? std::abs(fcy(i,j,k,0)) : T(0.0);
            T fracz = (ccm(i,j-1,kk) || ccm(i,j,kk))
                ? std::abs(fcy(i,j,k,1)) : T(0.0);
            sym = (T(1.0)-fracx)*(T(1.0)-fracz)*sym;
        }

        T syp = -bY(i,j+1,k,n);
        if (apyp != T(0.0) && apyp != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcy(i,j+1,k,0)));
            int kk = k + static_cast<int>(std::copysign(T(1.0),fcy(i,j+1,k,1)));
            T fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                ? std::abs(fcy(i,j+1,k,0)) : T(0.0);
            T fracz = (ccm(i,j,kk) || ccm(i,j+1,kk))
                ? std::abs(fcy(i,j+1,k,1)) : T(0.0);
            syp = (T(1.0)-fracx)*(T(1.0)-fracz)*syp;
        }

        T szm =  bZ(i,j,k,n);
        if (apzm != T(0.0) && apzm != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,0)));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k,1)));
            T fracx = (ccm(ii,j,k-1) || ccm(ii,j,k))
                ? std::abs(fcz(i,j,k,0)) : T(0.0);
            T fracy = (ccm(i,jj,k-1) || ccm(i,jj,k))
                ? std::abs(fcz(i,j,k,1)) : T(0.0);
            szm = (T(1.0)-fracx)*(T(1.0)-fracy)*szm;
        }

        T szp = -bZ(i,j,k+1,n);
        if (apzp != T(0.0) && apzp != T(1.0) && !beta_on_centroid) {
            int ii = i + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k+1,0)));
            int jj = j + static_cast<int>(std::copysign(T(1.0),fcz(i,j,k+1,1)));
            T fracx = (ccm(ii,j,k) || ccm(ii,j,k+1))
                ? std::abs(fcz(i,j,k+1,0)) : T(0.0);
            T fracy = (ccm(i,jj,k) || ccm(i,jj,k+1))
                ? std::abs(fcz(i,j,k+1,1)) : T(0.0);
            szp = (T(1.0)-fracx)*(T(1.0)-fracy)*szp;
        }

        T vfrcinv = T(1.0)/kappa;
        T gamma = alpha*a(i,j,k) + vfrcinv *
            (dhx*(apxm*sxm-apxp*sxp) +
             dhy*(apym*sym-apyp*syp) +
             dhz*(apzm*szm-apzp*szp));

        if (is_dirichlet) {
            T dapx = apxm-apxp;
            T dapy = apym-apyp;
            T dapz = apzm-apzp;
            T anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
            T anorminv = T(1.0)/anorm;
            T anrmx = dapx * anorminv;
            T anrmy = dapy * anorminv;
            T anrmz = dapz * anorminv;
            T bctx = bc(i,j,k,0);
            T bcty = bc(i,j,k,1);
            T bctz = bc(i,j,k,2);
            T dx_eb = get_dx_eb(vfrc(i,j,k));

            T dg = dx_eb / amrex::max(std::abs(anrmx),std::abs(anrmy),
                                         std::abs(anrmz));

            T gx = bctx - dg*anrmx;
            T gy = bcty - dg*anrmy;
            T gz = bctz - dg*anrmz;
            T sx = std::copysign(T(1.0),anrmx);
            T sy = std::copysign(T(1.0),anrmy);
            T sz = std::copysign(T(1.0),anrmz);

            gx *= sx;
            gy *= sy;
            gz *= sz;
            T gxy = gx*gy;
            T gxz = gx*gz;
            T gyz = gy*gz;
            T gxyz = gx*gy*gz;
            T phig_gamma = (T(1.0)+gx+gy+gz+gxy+gxz+gyz+gxyz);
            T feb_gamma = -phig_gamma/dg * ba(i,j,k) * beb(i,j,k,n);
            gamma += vfrcinv*(-dhx)*feb_gamma;
        }

        phi(i,j,k,n) /= gamma;
    }
}

}

#endif
