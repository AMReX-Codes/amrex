#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLTensor_K.H>
#include <AMReX_MLABecLap_K.H>

namespace amrex {

namespace {
    constexpr int kappa_num_mglevs = 1;
}

template <typename MF>
MLTensorOpT<MF>::MLTensorOpT ()
{
    MLABecLaplacianT<MF>::setScalars(1.0,1.0);
}

template <typename MF>
MLTensorOpT<MF>::MLTensorOpT (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FAB> const*>& a_factory)
{
    MLABecLaplacianT<MF>::setScalars(1.0,1.0);
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
MLTensorOpT<MF>::MLTensorOpT (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const Vector<iMultiFab const*>& a_overset_mask,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FAB> const*>& a_factory)
{
    MLABecLaplacianT<MF>::setScalars(1.0,1.0);
    define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);
}

template <typename MF>
void
MLTensorOpT<MF>::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLTensorOp::define()");

    MLABecLaplacianT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory,
                            AMREX_SPACEDIM);

    m_kappa.clear();
    m_mapfac.clear();
    m_use_mapped = false;
    m_kappa.resize(this->NAMRLevels());
    for (int amrlev = 0; amrlev < this->NAMRLevels(); ++amrlev) {
        m_kappa[amrlev].resize(std::min(kappa_num_mglevs,this->NMGLevels(amrlev)));
        for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_kappa[amrlev][mglev][idim].define
                    (amrex::convert(this->m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     this->m_dmap[amrlev][mglev], 1, 0,
                     MFInfo(), *this->m_factory[amrlev][mglev]);
            }
        }
    }

    this->makeSPLinOp(a_info, [&] (auto const& sp_info) {
        return std::make_unique<MLTensorOpT<fMultiFab>>(this->spGeom(a_geom), a_grids, a_dmap, sp_info);
    });
}

template <typename MF>
void
MLTensorOpT<MF>::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const Vector<iMultiFab const*>& a_overset_mask,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLTensorOp::define(oveset)");

    MLABecLaplacianT<MF>::define(a_geom, a_grids, a_dmap, a_overset_mask, a_info,
                            a_factory, AMREX_SPACEDIM);

    m_kappa.clear();
    m_mapfac.clear();
    m_use_mapped = false;
    m_kappa.resize(this->NAMRLevels());
    for (int amrlev = 0; amrlev < this->NAMRLevels(); ++amrlev) {
        m_kappa[amrlev].resize(std::min(kappa_num_mglevs,this->NMGLevels(amrlev)));
        for (int mglev = 0; mglev < m_kappa[amrlev].size(); ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_kappa[amrlev][mglev][idim].define
                    (amrex::convert(this->m_grids[amrlev][mglev],
                                    IntVect::TheDimensionVector(idim)),
                     this->m_dmap[amrlev][mglev], 1, 0,
                     MFInfo(), *this->m_factory[amrlev][mglev]);
            }
        }
    }

    this->makeSPLinOp(a_info, [&] (auto const& sp_info) {
        return std::make_unique<MLTensorOpT<fMultiFab>>(this->spGeom(a_geom), a_grids, a_dmap,
                                                        a_overset_mask, sp_info);
    });
}

template <typename MF>
void
MLTensorOpT<MF>::setShearViscosity (int amrlev, const Array<MF const*,AMREX_SPACEDIM>& eta)
{
    MLABecLaplacianT<MF>::setBCoeffs(amrlev, eta);
}

template <typename MF>
void
MLTensorOpT<MF>::setShearViscosity (int amrlev, RT eta)
{
    MLABecLaplacianT<MF>::setBCoeffs(amrlev, eta);
}


template <typename MF>
void
MLTensorOpT<MF>::setBulkViscosity (int amrlev, RT kappa)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_kappa[amrlev][0][idim].setVal(kappa);
    }
    m_has_kappa = true;
    m_needs_update = true;
    this->template spForward<MLTensorOpT>([&] (auto* sp) {
        sp->setBulkViscosity(amrlev, static_cast<float>(kappa)); });
}


template <typename MF>
void
MLTensorOpT<MF>::prepareForSolve ()
{
    if (m_use_mapped) {
        // The factors are set per AMR level; a level without them would
        // silently use the unmapped cross-term kernels.
        for (int amrlev = 0; amrlev < this->NAMRLevels(); ++amrlev) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hasMappingFactors(amrlev),
                "MLTensorOp: setMappingFactors must be called for every AMR level");
        }
    }
    if (m_has_kappa) {
        for (int amrlev = this->NAMRLevels()-1; amrlev >= 0; --amrlev) {
            for (int mglev = 1; mglev < m_kappa[amrlev].size(); ++mglev) {
                amrex::average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev][mglev-1]),
                                          GetArrOfPtrs     (m_kappa[amrlev][mglev  ]),
                                          IntVect(this->mg_coarsen_ratio), 0);
            }
            if (amrlev > 0) {
                amrex::average_down_faces(GetArrOfConstPtrs(m_kappa[amrlev  ].back()),
                                          GetArrOfPtrs     (m_kappa[amrlev-1].front()),
                                          this->AMRRefRatioVect(amrlev-1), this->m_geom[amrlev-1][0]);
            }
        }
    } else {
        for (int amrlev = 0; amrlev < this->NAMRLevels(); ++amrlev) {
            for (auto & mglev : m_kappa[amrlev]) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    mglev[idim].setVal(0.0);
                }
            }
        }
    }

    for (int amrlev = 0; amrlev < this->NAMRLevels(); ++amrlev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int icomp = idim;
            MF::Xpay(this->m_b_coeffs[amrlev][0][idim], RT(4./3.),
                           m_kappa[amrlev][0][idim], 0, icomp, 1, IntVect(0));
        }
    }

    MLABecLaplacianT<MF>::prepareForSolve();

    for (int amrlev = this->NAMRLevels()-1; amrlev >= 0; --amrlev) {
        for (int mglev = 1; mglev < m_kappa[amrlev].size(); ++mglev) {
            if (m_has_kappa && this->m_overset_mask[amrlev][mglev]) {
                const RT fac = static_cast<RT>(1 << mglev); // 2**mglev
                const RT osfac = RT(2.0)*fac/(fac+RT(1.0));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(m_kappa[amrlev][mglev][0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    AMREX_D_TERM(Box const& xbx = mfi.nodaltilebox(0);,
                                 Box const& ybx = mfi.nodaltilebox(1);,
                                 Box const& zbx = mfi.nodaltilebox(2));
                    AMREX_D_TERM(Array4<RT> const& bx = m_kappa[amrlev][mglev][0].array(mfi);,
                                 Array4<RT> const& by = m_kappa[amrlev][mglev][1].array(mfi);,
                                 Array4<RT> const& bz = m_kappa[amrlev][mglev][2].array(mfi));
                    Array4<int const> const& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
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

    m_needs_update = false;
}

template <typename MF>
void
MLTensorOpT<MF>::apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                   StateMode s_mode, const MLMGBndryT<MF>* bndry) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev,mglev,out,in,bc_mode,s_mode,bndry);
#else
    BL_PROFILE("MLTensorOp::apply()");

    MLABecLaplacianT<MF>::apply(amrlev, mglev, out, in, bc_mode, s_mode, bndry);

    if (mglev >= m_kappa[amrlev].size()) { return; }

    applyBCTensor(amrlev, mglev, in, bc_mode, s_mode, bndry);

    const auto& bcondloc = *this->m_bcondloc[amrlev][mglev];

    Array4<RT const> foo;

    const auto dxinv_r = this->m_geom[amrlev][mglev].InvCellSizeArray();
    const GpuArray<RT,AMREX_SPACEDIM> dxinv{AMREX_D_DECL(RT(dxinv_r[0]),RT(dxinv_r[1]),RT(dxinv_r[2]))};
    const Box& domain = this->m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Array<MF,AMREX_SPACEDIM> const& etamf = this->m_b_coeffs[amrlev][mglev];
    Array<MF,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    RT bscalar = this->m_b_scalar;
    const bool mapped = hasMappingFactors(amrlev) && (mglev == 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FAB fluxfab_tmp[AMREX_SPACEDIM];
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> const axfab = out.array(mfi);
            Array4<RT const> const vfab = in.const_array(mfi);
            AMREX_D_TERM(Array4<RT const> const etaxfab = etamf[0].const_array(mfi);,
                         Array4<RT const> const etayfab = etamf[1].const_array(mfi);,
                         Array4<RT const> const etazfab = etamf[2].const_array(mfi););
            AMREX_D_TERM(Array4<RT const> const kapxfab = kapmf[0].const_array(mfi);,
                         Array4<RT const> const kapyfab = kapmf[1].const_array(mfi);,
                         Array4<RT const> const kapzfab = kapmf[2].const_array(mfi););
            // Mesh-mapping face factors (empty Array4s when not mapped)
            AMREX_D_TERM(Array4<RT const> const mapxfab = mapped ? m_mapfac[amrlev][0].const_array(mfi) : foo;,
                         Array4<RT const> const mapyfab = mapped ? m_mapfac[amrlev][1].const_array(mfi) : foo;,
                         Array4<RT const> const mapzfab = mapped ? m_mapfac[amrlev][2].const_array(mfi) : foo;);
            AMREX_D_TERM(Box const xbx = amrex::surroundingNodes(bx,0);,
                         Box const ybx = amrex::surroundingNodes(bx,1);,
                         Box const zbx = amrex::surroundingNodes(bx,2););
            AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,AMREX_SPACEDIM, The_Async_Arena());,
                         fluxfab_tmp[1].resize(ybx,AMREX_SPACEDIM, The_Async_Arena());,
                         fluxfab_tmp[2].resize(zbx,AMREX_SPACEDIM, The_Async_Arena()););
            AMREX_D_TERM(Array4<RT> const fxfab = fluxfab_tmp[0].array();,
                         Array4<RT> const fyfab = fluxfab_tmp[1].array();,
                         Array4<RT> const fzfab = fluxfab_tmp[2].array(););

            if (domain.strictly_contains(bx)) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv,mapxfab);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv,mapyfab);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv,mapzfab);
                  }
                );
            } else {
                const auto & bdcv = bcondloc.bndryConds(mfi);

                Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    for (OrientationIter face; face; ++face) {
                        Orientation ori = face();
                        bct(ori,icomp) = bdcv[icomp][ori];
                    }
                }

                const auto& bvxlo = (bndry != nullptr) ?
                    (*bndry)[Orientation(0,Orientation::low )].array(mfi) : foo;
                const auto& bvylo = (bndry != nullptr) ?
                    (*bndry)[Orientation(1,Orientation::low )].array(mfi) : foo;
                const auto& bvxhi = (bndry != nullptr) ?
                    (*bndry)[Orientation(0,Orientation::high)].array(mfi) : foo;
                const auto& bvyhi = (bndry != nullptr) ?
                    (*bndry)[Orientation(1,Orientation::high)].array(mfi) : foo;
#if (AMREX_SPACEDIM == 3)
                const auto& bvzlo = (bndry != nullptr) ?
                    (*bndry)[Orientation(2,Orientation::low )].array(mfi) : foo;
                const auto& bvzhi = (bndry != nullptr) ?
                    (*bndry)[Orientation(2,Orientation::high)].array(mfi) : foo;
#endif

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv,
                                              bvxlo, bvxhi, bct, dlo, dhi, mapxfab);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv,
                                              bvylo, bvyhi, bct, dlo, dhi, mapyfab);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv,
                                              bvzlo, bvzhi, bct, dlo, dhi, mapzfab);
                  }
                );
            }

            if (this->m_overset_mask[amrlev][mglev]) {
                const auto& osm = this->m_overset_mask[amrlev][mglev]->array(mfi);
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

template <typename MF>
void
MLTensorOpT<MF>::applyBCTensor (int amrlev, int mglev, MF& vel, // NOLINT(readability-convert-member-functions-to-static)
                           BCMode bc_mode, StateMode, const MLMGBndryT<MF>* bndry) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev,mglev,vel,bc_mode,bndry);
#else

    const int inhomog = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = this->maxorder;
    const auto& bcondloc = *this->m_bcondloc[amrlev][mglev];
    const auto& maskvals = this->m_maskvals[amrlev][mglev];

    Array4<RT const> foo;

    const auto dxinv_r = this->m_geom[amrlev][mglev].InvCellSizeArray();
    const GpuArray<RT,AMREX_SPACEDIM> dxinv{AMREX_D_DECL(RT(dxinv_r[0]),RT(dxinv_r[1]),RT(dxinv_r[2]))};
    const Box& domain = this->m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();

        const auto& velfab = vel.array(mfi);

        const auto & bdlv = bcondloc.bndryLocs(mfi);
        const auto & bdcv = bcondloc.bndryConds(mfi);

        Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
        Array2D<RT,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bcl;
        for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
            for (OrientationIter face; face; ++face) {
                Orientation ori = face();
                bct(ori,icomp) = bdcv[icomp][ori];
                bcl(ori,icomp) = bdlv[icomp][ori];
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
                                  dxinv, dlo, dhi);
        });
#else
        const auto& mzlo = maskvals[Orientation(2,Orientation::low )].array(mfi);
        const auto& mzhi = maskvals[Orientation(2,Orientation::high)].array(mfi);

        const auto& bvzlo = (bndry != nullptr) ?
          (*bndry)[Orientation(2,Orientation::low )].array(mfi) : foo;
        const auto& bvzhi = (bndry != nullptr) ?
          (*bndry)[Orientation(2,Orientation::high)].array(mfi) : foo;

        // only edge vals used in 3D stencil
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            amrex::launch<64>(12, Gpu::gpuStream(),
#ifdef AMREX_USE_SYCL
            [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item)
            {
                int bid = item.get_group_linear_id();
                int tid = item.get_local_linear_id();
                int bdim = item.get_local_range(0);
#else
            [=] AMREX_GPU_DEVICE ()
            {
                int bid = blockIdx.x;
                int tid = threadIdx.x;
                int bdim = blockDim.x;
#endif
                mltensor_fill_edges(bid, tid, bdim, vbx, velfab,
                                    mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                    bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                    bct, bcl, inhomog, imaxorder,
                                    dxinv, dlo, dhi);
            });
        } else
#endif
        {
            mltensor_fill_edges(vbx, velfab,
                                mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                bct, bcl, inhomog, imaxorder,
                                dxinv, dlo, dhi);
        }
#endif
    }

#endif
}



template <typename MF>
void
MLTensorOpT<MF>::compFlux (int amrlev, const Array<MF*,AMREX_SPACEDIM>& fluxes,
                       MF& sol, Location loc) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev, fluxes, sol, loc);
#else
    BL_PROFILE("MLTensorOp::compFlux()");

    const int mglev = 0;
    const int ncomp = getNComp();
    MLABecLaplacianT<MF>::compFlux(amrlev, fluxes, sol, loc);

    MLMGBndryT<MF> const* bndry = this->m_bndry_sol[amrlev].get();
    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);

    const auto& bcondloc = *this->m_bcondloc[amrlev][mglev];

    const auto dxinv_r = this->m_geom[amrlev][mglev].InvCellSizeArray();
    const GpuArray<RT,AMREX_SPACEDIM> dxinv{AMREX_D_DECL(RT(dxinv_r[0]),RT(dxinv_r[1]),RT(dxinv_r[2]))};
    const Box& domain = this->m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Array<MF,AMREX_SPACEDIM> const& etamf = this->m_b_coeffs[amrlev][mglev];
    Array<MF,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    RT bscalar = this->m_b_scalar;
    const bool mapped = hasMappingFactors(amrlev);
    Array4<RT const> foo;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FAB,AMREX_SPACEDIM> fluxfab_tmp;

        for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<RT const> const vfab = sol.const_array(mfi);
            AMREX_D_TERM(Array4<RT const> const etaxfab = etamf[0].const_array(mfi);,
                         Array4<RT const> const etayfab = etamf[1].const_array(mfi);,
                         Array4<RT const> const etazfab = etamf[2].const_array(mfi););
            AMREX_D_TERM(Array4<RT const> const kapxfab = kapmf[0].const_array(mfi);,
                         Array4<RT const> const kapyfab = kapmf[1].const_array(mfi);,
                         Array4<RT const> const kapzfab = kapmf[2].const_array(mfi););
            AMREX_D_TERM(Array4<RT const> const mapxfab = mapped ? m_mapfac[amrlev][0].const_array(mfi) : foo;,
                         Array4<RT const> const mapyfab = mapped ? m_mapfac[amrlev][1].const_array(mfi) : foo;,
                         Array4<RT const> const mapzfab = mapped ? m_mapfac[amrlev][2].const_array(mfi) : foo;);
            AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                         Box const ybx = mfi.nodaltilebox(1);,
                         Box const zbx = mfi.nodaltilebox(2););
            AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,AMREX_SPACEDIM, The_Async_Arena());,
                         fluxfab_tmp[1].resize(ybx,AMREX_SPACEDIM, The_Async_Arena());,
                         fluxfab_tmp[2].resize(zbx,AMREX_SPACEDIM, The_Async_Arena()););
            AMREX_D_TERM(Array4<RT> const fxfab = fluxfab_tmp[0].array();,
                         Array4<RT> const fyfab = fluxfab_tmp[1].array();,
                         Array4<RT> const fzfab = fluxfab_tmp[2].array(););

            if (domain.strictly_contains(mfi.tilebox())) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv,mapxfab);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv,mapyfab);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv,mapzfab);
                  }
                );
            } else {
                const auto & bdcv = bcondloc.bndryConds(mfi);

                Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    for (OrientationIter face; face; ++face) {
                        Orientation ori = face();
                        bct(ori,icomp) = bdcv[icomp][ori];
                    }
                }

                const auto& bvxlo = (*bndry)[Orientation(0,Orientation::low )].array(mfi);
                const auto& bvylo = (*bndry)[Orientation(1,Orientation::low )].array(mfi);
                const auto& bvxhi = (*bndry)[Orientation(0,Orientation::high)].array(mfi);
                const auto& bvyhi = (*bndry)[Orientation(1,Orientation::high)].array(mfi);
#if (AMREX_SPACEDIM == 3)
                const auto& bvzlo = (*bndry)[Orientation(2,Orientation::low )].array(mfi);
                const auto& bvzhi = (*bndry)[Orientation(2,Orientation::high)].array(mfi);
#endif

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv,
                                              bvxlo, bvxhi, bct, dlo, dhi, mapxfab);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv,
                                              bvylo, bvyhi, bct, dlo, dhi, mapyfab);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv,
                                              bvzlo, bvzhi, bct, dlo, dhi, mapzfab);
                  }
                );
            }

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                Array4<RT      > dst = fluxes[idim]->array(mfi);
                Array4<RT const> src = fluxfab_tmp[idim].const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, ncomp, i, j, k, n,
                {
                    dst(i,j,k,n) += bscalar*src(i,j,k,n);
                });
            }

        }
    }
#endif
}



template <typename MF>
void
MLTensorOpT<MF>::compVelGrad (int amrlev, const Array<MF*,AMREX_SPACEDIM>& fluxes, // NOLINT(readability-convert-member-functions-to-static)
                       MF& sol, Location /*loc*/) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev,fluxes,sol);
#else
    BL_PROFILE("MLTensorOp::compVelGrad()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!isMapped(),
        "MLTensorOp::compVelGrad is not mapping-aware; apply the factors in the caller");

    const int mglev = 0;

    MLMGBndryT<MF> const* bndry = this->m_bndry_sol[amrlev].get();
    this->applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);
    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);

    const auto& bcondloc = *this->m_bcondloc[amrlev][mglev];

    const auto dxinv_r = this->m_geom[amrlev][mglev].InvCellSizeArray();
    const GpuArray<RT,AMREX_SPACEDIM> dxinv{AMREX_D_DECL(RT(dxinv_r[0]),RT(dxinv_r[1]),RT(dxinv_r[2]))};
    const Box& domain = this->m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<RT const> const vfab = sol.const_array(mfi);
        AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                     Box const ybx = mfi.nodaltilebox(1);,
                     Box const zbx = mfi.nodaltilebox(2);)
        AMREX_D_TERM(Array4<RT> const fxfab = fluxes[0]->array(mfi);,
                     Array4<RT> const fyfab = fluxes[1]->array(mfi);,
                     Array4<RT> const fzfab = fluxes[2]->array(mfi);)

// The derivatives are put in the array with the following order:
// component: 0    ,  1    ,  2    ,  3    ,  4    , 5    ,  6    ,  7    ,  8
// in 2D:     dU/dx,  dV/dx,  dU/dy,  dV/dy
// in 3D:     dU/dx,  dV/dx,  dW/dx,  dU/dy,  dV/dy, dW/dy,  dU/dz,  dV/dz,  dW/dz

        if (domain.strictly_contains(mfi.tilebox())) {
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
        } else {
            const auto & bdcv = bcondloc.bndryConds(mfi);

            Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                for (OrientationIter face; face; ++face) {
                    Orientation ori = face();
                    bct(ori,icomp) = bdcv[icomp][ori];
                }
            }

            const auto& bvxlo = (*bndry)[Orientation(0,Orientation::low )].array(mfi);
            const auto& bvylo = (*bndry)[Orientation(1,Orientation::low )].array(mfi);
            const auto& bvxhi = (*bndry)[Orientation(0,Orientation::high)].array(mfi);
            const auto& bvyhi = (*bndry)[Orientation(1,Orientation::high)].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const auto& bvzlo = (*bndry)[Orientation(2,Orientation::low )].array(mfi);
            const auto& bvzhi = (*bndry)[Orientation(2,Orientation::high)].array(mfi);
#endif
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_vel_grads_fx(txbx,fxfab,vfab,dxinv,bvxlo,bvxhi,bct,dlo,dhi);
              }
            , ybx, tybx,
              {
                  mltensor_vel_grads_fy(tybx,fyfab,vfab,dxinv,bvylo,bvyhi,bct,dlo,dhi);
              }
            , zbx, tzbx,
              {
                  mltensor_vel_grads_fz(tzbx,fzfab,vfab,dxinv,bvzlo,bvzhi,bct,dlo,dhi);
              }
            );
        }
    }
#endif
}

template class MLTensorOpT<MultiFab>;
#ifdef AMREX_MLMG_HAS_FLOAT
template class MLTensorOpT<fMultiFab>; // single precision copy for mixed precision
#endif

}
