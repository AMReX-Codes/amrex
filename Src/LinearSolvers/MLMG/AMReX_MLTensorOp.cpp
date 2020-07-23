#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLTensor_K.H>
#include <AMReX_MLABecLap_K.H>

namespace amrex {

namespace {
    constexpr int kappa_num_mglevs = 1;
}

MLTensorOp::MLTensorOp ()
{
    MLABecLaplacian::setScalars(1.0,1.0);
}

MLTensorOp::MLTensorOp (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLABecLaplacian::setScalars(1.0,1.0);
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLTensorOp::MLTensorOp (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const Vector<iMultiFab const*>& a_overset_mask,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLABecLaplacian::setScalars(1.0,1.0);
    define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);
}

MLTensorOp::~MLTensorOp ()
{}

void
MLTensorOp::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLTensorOp::define()");

    MLABecLaplacian::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    m_kappa.clear();
    m_kappa.resize(NAMRLevels());
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        m_kappa[amrlev].resize(std::min(kappa_num_mglevs,NMGLevels(amrlev)));
        for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_kappa[amrlev][mglev][idim].define
                    (amrex::convert(m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     m_dmap[amrlev][mglev], 1, 0,
                     MFInfo(), *m_factory[amrlev][mglev]);
            }
        }
    }
}

void
MLTensorOp::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const Vector<iMultiFab const*>& a_overset_mask,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLTensorOp::define(oveset)");

    MLABecLaplacian::define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);

    m_kappa.clear();
    m_kappa.resize(NAMRLevels());
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        m_kappa[amrlev].resize(std::min(kappa_num_mglevs,NMGLevels(amrlev)));
        for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_kappa[amrlev][mglev][idim].define
                    (amrex::convert(m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     m_dmap[amrlev][mglev], 1, 0,
                     MFInfo(), *m_factory[amrlev][mglev]);
            }
        }
    }
}

void
MLTensorOp::setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta)
{
    MLABecLaplacian::setBCoeffs(amrlev, eta);
}

void
MLTensorOp::setShearViscosity (int amrlev, Real eta)
{
    MLABecLaplacian::setBCoeffs(amrlev, eta);
}

void
MLTensorOp::setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_kappa[amrlev][0][idim], *kappa[idim], 0, 0, 1, 0);
    }
    m_has_kappa = true;
}

void
MLTensorOp::setBulkViscosity (int amrlev, Real kappa)
{
    if (kappa != 0.0) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_kappa[amrlev][0][idim].setVal(kappa);
        }
        m_has_kappa = true;
    }
}

void
MLTensorOp::prepareForSolve ()
{
    if (m_has_kappa) {
        for (int amrlev = NAMRLevels()-1; amrlev >= 0; --amrlev) {
            for (int mglev = 1; mglev < m_kappa[amrlev].size(); ++mglev) {
                amrex::average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev][mglev-1]),
                                          GetArrOfPtrs     (m_kappa[amrlev][mglev  ]),
                                          IntVect(mg_coarsen_ratio), 0);
            }
            if (amrlev > 0) {
                amrex::average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev  ].back()),
                                          GetArrOfPtrs     (m_kappa[amrlev-1].front()),
                                          IntVect(mg_coarsen_ratio), 0);
            }
        }
    } else {
        for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
            for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    m_kappa[amrlev][mglev][idim].setVal(0.0);
                }
            }
        }
    }

    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int icomp = idim;
            MultiFab::Xpay(m_b_coeffs[amrlev][0][idim], 4./3.,
                           m_kappa[amrlev][0][idim], 0, icomp, 1, 0);
        }
    }

    MLABecLaplacian::prepareForSolve();

    for (int amrlev = NAMRLevels()-1; amrlev >= 0; --amrlev) {
        for (int mglev = 1; mglev < m_kappa[amrlev].size(); ++mglev) {
            if (m_has_kappa and m_overset_mask[amrlev][mglev]) {
                const Real fac = static_cast<Real>(1 << mglev); // 2**mglev
                const Real osfac = 2.0*fac/(fac+1.0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(m_kappa[amrlev][mglev][0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    AMREX_D_TERM(Box const& xbx = mfi.nodaltilebox(0);,
                                 Box const& ybx = mfi.nodaltilebox(1);,
                                 Box const& zbx = mfi.nodaltilebox(2));
                    AMREX_D_TERM(Array4<Real> const& bx = m_kappa[amrlev][mglev][0].array(mfi);,
                                 Array4<Real> const& by = m_kappa[amrlev][mglev][1].array(mfi);,
                                 Array4<Real> const& bz = m_kappa[amrlev][mglev][2].array(mfi));
                    Array4<int const> const& osm = m_overset_mask[amrlev][mglev]->const_array(mfi);
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                        (xbx, t_xbx,
                         {
                             overset_rescale_bcoef_x(t_xbx, bx, osm, 1, osfac);
                         },
                         ybx, t_ybx,
                         {
                             overset_rescale_bcoef_y(t_ybx, by, osm, 1, osfac);
                         },
                         zbx, t_zbx,
                         {
                             overset_rescale_bcoef_z(t_zbx, bz, osm, 1, osfac);
                         });
                }
            }
        }
    }
}

void
MLTensorOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                   StateMode s_mode, const MLMGBndry* bndry) const
{
#if (AMREX_SPACEDIM > 1)
    BL_PROFILE("MLTensorOp::apply()");

    MLABecLaplacian::apply(amrlev, mglev, out, in, bc_mode, s_mode, bndry);

    if (mglev >= m_kappa[amrlev].size()) return;

    applyBCTensor(amrlev, mglev, in, bc_mode, s_mode, bndry );

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();

    Array<MultiFab,AMREX_SPACEDIM> const& etamf = m_b_coeffs[amrlev][mglev];
    Array<MultiFab,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    Real bscalar = m_b_scalar;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox fluxfab_tmp[AMREX_SPACEDIM];
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const axfab = out.array(mfi);
            Array4<Real const> const vfab = in.const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const etaxfab = etamf[0].const_array(mfi);,
                         Array4<Real const> const etayfab = etamf[1].const_array(mfi);,
                         Array4<Real const> const etazfab = etamf[2].const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const kapxfab = kapmf[0].const_array(mfi);,
                         Array4<Real const> const kapyfab = kapmf[1].const_array(mfi);,
                         Array4<Real const> const kapzfab = kapmf[2].const_array(mfi););
            AMREX_D_TERM(Box const xbx = amrex::surroundingNodes(bx,0);,
                         Box const ybx = amrex::surroundingNodes(bx,1);,
                         Box const zbx = amrex::surroundingNodes(bx,2););
            AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,AMREX_SPACEDIM);,
                         fluxfab_tmp[1].resize(ybx,AMREX_SPACEDIM);,
                         fluxfab_tmp[2].resize(zbx,AMREX_SPACEDIM););
            AMREX_D_TERM(Elixir fxeli = fluxfab_tmp[0].elixir();,
                         Elixir fyeli = fluxfab_tmp[1].elixir();,
                         Elixir fzeli = fluxfab_tmp[2].elixir(););
            AMREX_D_TERM(Array4<Real> const fxfab = fluxfab_tmp[0].array();,
                         Array4<Real> const fyfab = fluxfab_tmp[1].array();,
                         Array4<Real> const fzfab = fluxfab_tmp[2].array(););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv);
              }
            , ybx, tybx,
              {
                  mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv);
              }
            , zbx, tzbx,
              {
                  mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv);
              }
            );

            if (m_overset_mask[amrlev][mglev]) {
                const auto& osm = m_overset_mask[amrlev][mglev]->array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    mltensor_cross_terms_os(tbx, axfab, AMREX_D_DECL(fxfab,fyfab,fzfab),
                                            osm, dxinv, bscalar);
                });
            } else {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    mltensor_cross_terms(tbx, axfab, AMREX_D_DECL(fxfab,fyfab,fzfab),
                                         dxinv, bscalar);
                });
            }
        }
    }
#endif
}

void
MLTensorOp::applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                           BCMode bc_mode, StateMode, const MLMGBndry* bndry) const
{
#if (AMREX_SPACEDIM > 1)
 
    const int inhomog = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = maxorder;
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];
    const auto& maskvals = m_maskvals[amrlev][mglev];

    FArrayBox foofab(Box::TheUnitBox(),3);
    const auto& foo = foofab.array();

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);

    // Domain and coarse-fine boundaries are handled below.

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();

        const auto& velfab = vel.array(mfi);

        const auto & bdlv = bcondloc.bndryLocs(mfi);
        const auto & bdcv = bcondloc.bndryConds(mfi);

        GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> bct;
        GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> bcl;
        for (OrientationIter face; face; ++face) {
            Orientation ori = face();
            const int iface = ori;
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                bct[iface*AMREX_SPACEDIM+icomp] = bdcv[icomp][ori];
                bcl[iface*AMREX_SPACEDIM+icomp] = bdlv[icomp][ori];
            }
        }

        const auto& mxlo = maskvals[Orientation(0,Orientation::low )].array(mfi);
        const auto& mylo = maskvals[Orientation(1,Orientation::low )].array(mfi);
        const auto& mxhi = maskvals[Orientation(0,Orientation::high)].array(mfi);
        const auto& myhi = maskvals[Orientation(1,Orientation::high)].array(mfi);

        const auto& bvxlo = (bndry != nullptr) ?
	  (*bndry)[Orientation(0,Orientation::low )].array(mfi) : foo;
        const auto& bvylo = (bndry != nullptr) ?
	  (*bndry)[Orientation(1,Orientation::low )].array(mfi) : foo;
        const auto& bvxhi = (bndry != nullptr) ?
	  (*bndry)[Orientation(0,Orientation::high)].array(mfi) : foo;
        const auto& bvyhi = (bndry != nullptr) ?
	  (*bndry)[Orientation(1,Orientation::high)].array(mfi) : foo;

#if (AMREX_SPACEDIM == 2)

        AMREX_HOST_DEVICE_FOR_1D ( 4, icorner,
        {
            mltensor_fill_corners(icorner, vbx, velfab,
                                  mxlo, mylo, mxhi, myhi,
                                  bvxlo, bvylo, bvxhi, bvyhi,
                                  bct, bcl, inhomog, imaxorder,
                                  dxinv, domain);
        });
#else
        const auto& mzlo = maskvals[Orientation(2,Orientation::low )].array(mfi);
        const auto& mzhi = maskvals[Orientation(2,Orientation::high)].array(mfi);

        const auto& bvzlo = (bndry != nullptr) ?
	  (*bndry)[Orientation(2,Orientation::low )].array(mfi) : foo;
        const auto& bvzhi = (bndry != nullptr) ?
	  (*bndry)[Orientation(2,Orientation::high)].array(mfi) : foo;

	// only edge vals used in 3D stencil
#ifdef AMREX_USE_DPCPP
        // xxxxx DPCPP todo: kernel size
        Vector<Array4<int const> > htmp = {mxlo,mylo,mzlo,mxhi,myhi,mzhi};
        Gpu::AsyncArray<Array4<int const> > dtmp(htmp.data(), 6);
        auto dp = dtmp.data();
        AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
        {
            mltensor_fill_edges(iedge, vbx, velfab,
                                dp[0],dp[1],dp[2],dp[3],dp[4],dp[5],
                                bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                bct, bcl, inhomog, imaxorder,
                                dxinv, domain);
        });
#else
        AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
        {
            mltensor_fill_edges(iedge, vbx, velfab,
                                mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                bct, bcl, inhomog, imaxorder,
                                dxinv, domain);
        });
#endif

#endif
    }

    vel.EnforcePeriodicity(0, AMREX_SPACEDIM, IntVect(1),
                           m_geom[amrlev][mglev].periodicity());
#endif
}

void
MLTensorOp::compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location loc) const
{
#if (AMREX_SPACEDIM > 1)
    BL_PROFILE("MLTensorOp::compFlux()");

    const int mglev = 0;
    const int ncomp = getNComp();
    MLABecLaplacian::compFlux(amrlev, fluxes, sol, loc);

    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, m_bndry_sol[amrlev].get());

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();

    Array<MultiFab,AMREX_SPACEDIM> const& etamf = m_b_coeffs[amrlev][mglev];
    Array<MultiFab,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    Real bscalar = m_b_scalar;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FArrayBox,AMREX_SPACEDIM> fluxfab_tmp;

        for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real const> const vfab = sol.const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const etaxfab = etamf[0].const_array(mfi);,
                         Array4<Real const> const etayfab = etamf[1].const_array(mfi);,
                         Array4<Real const> const etazfab = etamf[2].const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const kapxfab = kapmf[0].const_array(mfi);,
                         Array4<Real const> const kapyfab = kapmf[1].const_array(mfi);,
                         Array4<Real const> const kapzfab = kapmf[2].const_array(mfi););
            AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                         Box const ybx = mfi.nodaltilebox(1);,
                         Box const zbx = mfi.nodaltilebox(2););
	    AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,AMREX_SPACEDIM);,
                         fluxfab_tmp[1].resize(ybx,AMREX_SPACEDIM);,
                         fluxfab_tmp[2].resize(zbx,AMREX_SPACEDIM););
            AMREX_D_TERM(Elixir fxeli = fluxfab_tmp[0].elixir();,
                         Elixir fyeli = fluxfab_tmp[1].elixir();,
                         Elixir fzeli = fluxfab_tmp[2].elixir(););
            AMREX_D_TERM(Array4<Real> const fxfab = fluxfab_tmp[0].array();,
                         Array4<Real> const fyfab = fluxfab_tmp[1].array();,
                         Array4<Real> const fzfab = fluxfab_tmp[2].array(););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv);
              }
            , ybx, tybx,
              {
                  mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv);
              }
            , zbx, tzbx,
              {
                  mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv);
              }
            );

	    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	        const Box& nbx = mfi.nodaltilebox(idim);
                Array4<Real      > dst = fluxes[idim]->array(mfi);
		Array4<Real const> src = fluxfab_tmp[idim].const_array();
		AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, ncomp, i, j, k, n,
                {
                    dst(i,j,k,n) += bscalar*src(i,j,k,n);
                });
            }

        }
    }
#endif
}



void
MLTensorOp::compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location /*loc*/) const
{
#if (AMREX_SPACEDIM > 1)
    BL_PROFILE("MLTensorOp::compVelGrad()");

    const int mglev = 0;

    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, m_bndry_sol[amrlev].get());

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const int dim_fluxes = AMREX_SPACEDIM*AMREX_SPACEDIM;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FArrayBox,AMREX_SPACEDIM> fluxfab_tmp;

        for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real const> const vfab = sol.const_array(mfi);
            AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                         Box const ybx = mfi.nodaltilebox(1);,
                         Box const zbx = mfi.nodaltilebox(2););
            AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,dim_fluxes);,
                         fluxfab_tmp[1].resize(ybx,dim_fluxes);,
                         fluxfab_tmp[2].resize(zbx,dim_fluxes););
            AMREX_D_TERM(Elixir fxeli = fluxfab_tmp[0].elixir();,
                         Elixir fyeli = fluxfab_tmp[1].elixir();,
                         Elixir fzeli = fluxfab_tmp[2].elixir(););
            AMREX_D_TERM(Array4<Real> const fxfab = fluxfab_tmp[0].array();,
                         Array4<Real> const fyfab = fluxfab_tmp[1].array();,
                         Array4<Real> const fzfab = fluxfab_tmp[2].array(););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_vel_grads_fx(txbx,fxfab,vfab,dxinv);
              }
            , ybx, tybx,
              {
                  mltensor_vel_grads_fy(tybx,fyfab,vfab,dxinv);
              }
            , zbx, tzbx,
              {
                  mltensor_vel_grads_fz(tzbx,fzfab,vfab,dxinv);
              }
            );

// The derivatives are put in the array with the following order:
// component: 0    ,  1    ,  2    ,  3    ,  4    , 5    ,  6    ,  7    ,  8   
// in 2D:     dU/dx,  dV/dx,  dU/dy,  dV/dy 
// in 3D:     dU/dx,  dV/dx,  dW/dx,  dU/dy,  dV/dy, dW/dy,  dU/dz,  dV/dz,  dW/dz      
            
            
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                Array4<Real      > dst = fluxes[idim]->array(mfi);
                Array4<Real const> src = fluxfab_tmp[idim].const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, dim_fluxes, i, j, k, n,
                {
                    dst(i,j,k,n) = src(i,j,k,n);
                });
            }

        }
    }
#endif
}


}
