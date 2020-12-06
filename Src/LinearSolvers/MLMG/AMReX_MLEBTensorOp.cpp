#include <AMReX_MLEBTensorOp.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLTensor_K.H>
#include <AMReX_MLEBTensor_K.H>
#include <AMReX_MLEBABecLap.H>

namespace amrex {

namespace {
    constexpr int kappa_num_mglevs = 1;
}

MLEBTensorOp::MLEBTensorOp ()
{
    MLEBABecLap::setScalars(1.0,1.0);
}

MLEBTensorOp::MLEBTensorOp (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const LPInfo& a_info,
                            const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    MLEBABecLap::setScalars(1.0,1.0);
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLEBTensorOp::~MLEBTensorOp ()
{}

void
MLEBTensorOp::define (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info,
                      const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    BL_PROFILE("MLEBTensorOp::define()");

    MLEBABecLap::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    m_kappa.clear();
    m_kappa.resize(NAMRLevels());
    m_eb_kappa.resize(NAMRLevels());
    m_tauflux.resize(NAMRLevels());
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        m_kappa[amrlev].resize(std::min(kappa_num_mglevs,NMGLevels(amrlev)));
        m_eb_kappa[amrlev].resize(m_kappa[amrlev].size());
        m_tauflux[amrlev].resize(m_kappa[amrlev].size());
        for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_kappa[amrlev][mglev][idim].define
                    (amrex::convert(m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     m_dmap[amrlev][mglev], 1, 0,
                     MFInfo(), *m_factory[amrlev][mglev]);
                m_tauflux[amrlev][mglev][idim].define
                    (amrex::convert(m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     m_dmap[amrlev][mglev],
                     AMREX_SPACEDIM, IntVect(1)-IntVect::TheDimensionVector(idim),
                     MFInfo(), *m_factory[amrlev][mglev]);
                m_tauflux[amrlev][mglev][idim].setVal(0.0);
            }
            m_eb_kappa[amrlev][mglev].define(m_grids[amrlev][mglev],
                                             m_dmap[amrlev][mglev],
                                             1, 0, MFInfo(),
                                             *m_factory[amrlev][mglev]);
        }
    }
}

void
MLEBTensorOp::setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta,
                                 Location a_beta_loc)
{
    MLEBABecLap::setBCoeffs(amrlev, eta, a_beta_loc);
}

void
MLEBTensorOp::setShearViscosity (int amrlev, Real eta)
{
    MLEBABecLap::setBCoeffs(amrlev, eta);
}

void
MLEBTensorOp::setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_kappa[amrlev][0][idim], *kappa[idim], 0, 0, 1, 0);
    }
    m_has_kappa = true;
}

void
MLEBTensorOp::setBulkViscosity (int amrlev, Real kappa)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_kappa[amrlev][0][idim].setVal(kappa);
    }
    m_has_kappa = true;
}

void
MLEBTensorOp::setEBShearViscosity (int amrlev, MultiFab const& eta)
{
    MLEBABecLap::setEBHomogDirichlet(amrlev, eta);
}

void
MLEBTensorOp::setEBShearViscosity (int amrlev, Real eta)
{
    MLEBABecLap::setEBHomogDirichlet(amrlev, eta);
}

void
MLEBTensorOp::setEBBulkViscosity (int amrlev, MultiFab const& kappa)
{
    MultiFab::Copy(m_eb_kappa[amrlev][0], kappa, 0, 0, 1, 0);
    m_has_eb_kappa = true;
}

void
MLEBTensorOp::setEBBulkViscosity (int amrlev, Real kappa)
{
    if (kappa != 0.0) {
        m_eb_kappa[amrlev][0].setVal(kappa);
        m_has_eb_kappa = true;
    }
}

void
MLEBTensorOp::prepareForSolve ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_has_kappa == m_has_eb_kappa,
        "MLEBTensorOp: must call both setBulkViscosity and setEBBulkViscosity or none.");

    if (m_has_kappa) {
        for (int amrlev = NAMRLevels()-1; amrlev >= 0; --amrlev) {
            for (int mglev = 1; mglev < m_kappa[amrlev].size(); ++mglev) {
                amrex::EB_average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev][mglev-1]),
                                             GetArrOfPtrs     (m_kappa[amrlev][mglev  ]),
                                             IntVect(mg_coarsen_ratio), 0);
            }
            if (amrlev > 0) {
                amrex::EB_average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev  ].back()),
                                             GetArrOfPtrs     (m_kappa[amrlev-1].front()),
                                             IntVect(mg_coarsen_ratio), m_geom[amrlev-1][0]);
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

    if (m_has_eb_kappa) {
        for (int amrlev = NAMRLevels()-1; amrlev >= 0; --amrlev) {
            for (int mglev = 1; mglev < m_eb_kappa[amrlev].size(); ++mglev) {
                amrex::EB_average_down_boundaries(m_eb_kappa[amrlev][mglev-1],
                                                  m_eb_kappa[amrlev][mglev  ],
                                                  IntVect(mg_coarsen_ratio), 0);
            }
            if (amrlev > 0) {
                amrex::EB_average_down_boundaries(m_eb_kappa[amrlev  ].back(),
                                                  m_eb_kappa[amrlev-1].front(),
                                                  IntVect(mg_coarsen_ratio), 0);
            }
        }
    } else {
        for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
            for (int mglev = 0; mglev < m_eb_kappa[amrlev].size(); ++mglev) {
                m_eb_kappa[amrlev][mglev].setVal(0.0);
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

    MLEBABecLap::prepareForSolve();
}

void
MLEBTensorOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                     StateMode s_mode, const MLMGBndry* bndry) const
{
    BL_PROFILE("MLEBTensorOp::apply()");
    MLEBABecLap::apply(amrlev, mglev, out, in, bc_mode, s_mode, bndry);

    if (mglev >= m_kappa[amrlev].size()) return;

    applyBCTensor(amrlev, mglev, in, bc_mode, s_mode, bndry);

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* bcent = (factory) ? &(factory->getBndryCent()) : nullptr;

    const Geometry& geom = m_geom[amrlev][mglev];
    const auto dxinv = geom.InvCellSizeArray();

    Array<MultiFab,AMREX_SPACEDIM>& fluxmf = m_tauflux[amrlev][mglev];
    iMultiFab const& mask = m_cc_mask[amrlev][mglev];
    MultiFab const& etaebmf = *m_eb_b_coeffs[amrlev][mglev];
    MultiFab const& kapebmf = m_eb_kappa[amrlev][mglev];
    Real bscalar = m_b_scalar;

    compCrossTerms(amrlev, mglev, in);

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered) continue;

        Array4<Real> const axfab = out.array(mfi);
        AMREX_D_TERM(Array4<Real const> const fxfab = fluxmf[0].const_array(mfi);,
                     Array4<Real const> const fyfab = fluxmf[1].const_array(mfi);,
                     Array4<Real const> const fzfab = fluxmf[2].const_array(mfi););

        if (fabtyp == FabType::regular)
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mltensor_cross_terms(tbx, axfab, AMREX_D_DECL(fxfab,fyfab,fzfab), dxinv, bscalar);
            });
        }
        else
        {
            Array4<Real const> const& vfab = in.const_array(mfi);
            Array4<Real const> const& etab = etaebmf.const_array(mfi);
            Array4<Real const> const& kapb = kapebmf.const_array(mfi);
            Array4<int const> const& ccm = mask.const_array(mfi);
            Array4<EBCellFlag const> const& flag = flags->const_array(mfi);
            Array4<Real const> const& vol = vfrac->const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
                         Array4<Real const> const& apy = area[1]->const_array(mfi);,
                         Array4<Real const> const& apz = area[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););
            Array4<Real const> const& bc = bcent->const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlebtensor_cross_terms(tbx, axfab,
                                       AMREX_D_DECL(fxfab,fyfab,fzfab),
                                       vfab, etab, kapb, ccm, flag, vol,
                                       AMREX_D_DECL(apx,apy,apz),
                                       AMREX_D_DECL(fcx,fcy,fcz),
                                       bc, dxinv, bscalar);
            });
        }
    }
}

void
MLEBTensorOp::applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                             BCMode bc_mode, StateMode, const MLMGBndry* bndry) const
{
    const int inhomog = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = maxorder;
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];
    const auto& maskvals = m_maskvals[amrlev][mglev];

    FArrayBox foofab(Box::TheUnitBox(),3);
    const auto& foo = foofab.array();

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    
    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();

        auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;
        if (fabtyp != FabType::covered)
        {
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

            AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
            {
                mltensor_fill_edges(iedge, vbx, velfab,
                                    mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                    bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                    bct, bcl, inhomog, imaxorder,
                                    dxinv, domain);
            });

            AMREX_HOST_DEVICE_FOR_1D ( 8, icorner,
            {
                mltensor_fill_corners(icorner, vbx, velfab,
                                      mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                      bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                      bct, bcl, inhomog, imaxorder,
                                      dxinv, domain);
            });
#endif
        }
    }

    // Notet that it is incorrect to call EnforcePeriodicity on vel.
}

void
MLEBTensorOp::compCrossTerms(int amrlev, int mglev, MultiFab const& mf) const
{
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};

    const Geometry& geom = m_geom[amrlev][mglev];
    const auto dxinv = geom.InvCellSizeArray();

    Array<MultiFab,AMREX_SPACEDIM> const& etamf = m_b_coeffs[amrlev][mglev];
    Array<MultiFab,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    Array<MultiFab,AMREX_SPACEDIM>& fluxmf = m_tauflux[amrlev][mglev];
    
    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
	AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
		     Box const ybx = mfi.nodaltilebox(1);,
		     Box const zbx = mfi.nodaltilebox(2););

	// grow by 1 because of corners
	auto fabtyp = (flags) ? (*flags)[mfi].getType(amrex::grow(bx,1)) : FabType::regular;

	if (fabtyp == FabType::covered) {
	  AMREX_D_TERM(Array4<Real> const& fxfab = fluxmf[0].array(mfi);,
		       Array4<Real> const& fyfab = fluxmf[1].array(mfi);,
		       Array4<Real> const& fzfab = fluxmf[2].array(mfi););
	  AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
	  ( xbx, txbx,
	    {
                AMREX_LOOP_4D(txbx, AMREX_SPACEDIM, i, j, k, n,
                {
                    fxfab(i,j,k,n) = 0.0;
                });
	    }
	    , ybx, tybx,
	    {
                AMREX_LOOP_4D(tybx, AMREX_SPACEDIM, i, j, k, n,
                {
                    fyfab(i,j,k,n) = 0.0;
                });
	    }
	    , zbx, tzbx,
	    {
                AMREX_LOOP_4D(tzbx, AMREX_SPACEDIM, i, j, k, n,
                {
                    fzfab(i,j,k,n) = 0.0;
                });
	    }
	  );
	} else {
	  AMREX_D_TERM(Array4<Real> const fxfab = fluxmf[0].array(mfi);,
		       Array4<Real> const fyfab = fluxmf[1].array(mfi);,
		       Array4<Real> const fzfab = fluxmf[2].array(mfi););
	  Array4<Real const> const vfab = mf.const_array(mfi);
	  AMREX_D_TERM(Array4<Real const> const etaxfab = etamf[0].const_array(mfi);,
		       Array4<Real const> const etayfab = etamf[1].const_array(mfi);,
		       Array4<Real const> const etazfab = etamf[2].const_array(mfi););
	  AMREX_D_TERM(Array4<Real const> const kapxfab = kapmf[0].const_array(mfi);,
		       Array4<Real const> const kapyfab = kapmf[1].const_array(mfi);,
		       Array4<Real const> const kapzfab = kapmf[2].const_array(mfi););
	  
	  if (fabtyp == FabType::regular)
	  {
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
	  }
	  else
	  {
	    AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
			 Array4<Real const> const& apy = area[1]->const_array(mfi);,
			 Array4<Real const> const& apz = area[2]->const_array(mfi););
	    Array4<EBCellFlag const> const& flag = flags->const_array(mfi);
	    
	    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
	    ( xbx, txbx,
	      {
		mlebtensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,apx,flag,dxinv);
	      }
	      , ybx, tybx,
	      {
		mlebtensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,apy,flag,dxinv);
	      }
	      , zbx, tzbx,
	      {
		mlebtensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,apz,flag,dxinv);
	      }
	      );
	  }
	}
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxmf[idim].FillBoundary(0, AMREX_SPACEDIM, geom.periodicity());
    }
}
  
void
MLEBTensorOp::compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location loc) const
{
    BL_PROFILE("MLEBTensorOp::compFlux()");

    if ( !(loc==Location::FaceCenter || loc==Location::FaceCentroid) )
      amrex::Abort("MLEBTensorOp::compFlux() unknown location for fluxes.");

    const int mglev = 0;
    const int ncomp = getNComp();
    MLEBABecLap::compFlux(amrlev, fluxes, sol, loc);

    if (mglev >= m_kappa[amrlev].size()) return;

    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, m_bndry_sol[amrlev].get());

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};

    Array<MultiFab,AMREX_SPACEDIM>& fluxmf = m_tauflux[amrlev][mglev];
    Real bscalar = m_b_scalar;

    compCrossTerms(amrlev, mglev, sol);

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered) continue;

        if (fabtyp == FabType::regular)
        {
	    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	      const Box& nbx = mfi.nodaltilebox(idim);
	      Array4<Real      > dst = fluxes[idim]->array(mfi);
	      Array4<Real const> src = fluxmf[idim].array(mfi);
	      AMREX_HOST_DEVICE_FOR_4D (nbx, ncomp, i, j, k, n,
	      {
		  dst(i,j,k,n) += bscalar*src(i,j,k,n);
	      });
	    }
        }
        else if ( loc==Location::FaceCenter )
        {
	    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	      const Box& nbx = mfi.nodaltilebox(idim);
	      Array4<Real      > dst = fluxes[idim]->array(mfi);
	      Array4<Real const> src = fluxmf[idim].array(mfi);
	      Array4<Real const> const& ap = area[idim]->array(mfi);

	      AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( nbx, tbx,
	      {
	        mlebtensor_flux_0(tbx, dst, src, ap, bscalar);
	      });
	    }
        }
	else // loc==Location::FaceCentroid
	{
	    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];

	    AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                         Box const ybx = mfi.nodaltilebox(1);,
                         Box const zbx = mfi.nodaltilebox(2););
	    AMREX_D_TERM(Array4<Real const> fx = fluxmf[0].const_array(mfi);,
			 Array4<Real const> fy = fluxmf[1].const_array(mfi);,
			 Array4<Real const> fz = fluxmf[2].const_array(mfi););
	    AMREX_D_TERM(Array4<Real      > Ax = fluxes[0]->array(mfi);,
			 Array4<Real      > Ay = fluxes[1]->array(mfi);,
			 Array4<Real      > Az = fluxes[2]->array(mfi););

	    const auto& fcent = factory->getFaceCent();
	    AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
			 Array4<Real const> const& apy = area[1]->const_array(mfi);,
			 Array4<Real const> const& apz = area[2]->const_array(mfi););
	    AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
			 Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
			 Array4<Real const> const& fcz = fcent[2]->const_array(mfi););
	    Array4<int const> const& msk = ccmask.const_array(mfi);

	    int face_only = 0;

	    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM (
                xbx, txbx,
		{
		  mlebtensor_flux_x(txbx, Ax, fx, apx, fcx, bscalar, msk, face_only, xbx);
		}
		, ybx, tybx,
		{
		  mlebtensor_flux_y(tybx, Ay, fy, apy, fcy, bscalar, msk, face_only, ybx);
		}
		, zbx, tzbx,
		{
		  mlebtensor_flux_z(tzbx, Az, fz, apz, fcz, bscalar, msk, face_only, zbx);
		}
	    );

	}

    }
}

void
MLEBTensorOp::compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location loc) const
{
    BL_PROFILE("MLEBTensorOp::compVelGrad()");

    if ( !(loc==Location::FaceCenter || loc==Location::FaceCentroid) )
      amrex::Abort("MLEBTensorOp::compVelGrad() unknown location for VelGradients.");

    const int mglev = 0;

    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, m_bndry_sol[amrlev].get());

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;

    const Geometry& geom = m_geom[amrlev][mglev];
    const auto dxinv = geom.InvCellSizeArray();

    const int dim_fluxes = AMREX_SPACEDIM*AMREX_SPACEDIM;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    Array<FArrayBox,AMREX_SPACEDIM> fluxfab_tmp;
    for (MFIter mfi(sol, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered) continue;

        if (fabtyp == FabType::regular)
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
        else if ( loc==Location::FaceCenter )
        {

          amrex::Abort("compVelGrad not yet implemented for cut-cells  ");

        }
        else // loc==Location::FaceCentroid
        {
           
          amrex::Abort("compVelGrad not yet implemented for cut-cells  ");

        }

    }
  }
}

}
