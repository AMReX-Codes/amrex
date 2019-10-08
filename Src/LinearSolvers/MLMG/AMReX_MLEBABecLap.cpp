
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>

#include <AMReX_MLABecLap_K.H>
#include <AMReX_MLEBABecLap_K.H>
#include <AMReX_MLLinOp_K.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreABecLap3.H>
#endif

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

namespace amrex {

MLEBABecLap::MLEBABecLap (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info,
                          const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

std::unique_ptr<FabFactory<FArrayBox> >
MLEBABecLap::makeFactory (int amrlev, int mglev) const
{
    return makeEBFabFactory(m_geom[amrlev][mglev],
                            m_grids[amrlev][mglev],
                            m_dmap[amrlev][mglev],
                            {1,1,1}, EBSupport::full);
}

void
MLEBABecLap::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    BL_PROFILE("MLEBABecLap::define()");

    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    MLCellABecLap::define(a_geom, a_grids, a_dmap, a_info, _factory);

    const int ncomp = getNComp();

    m_a_coeffs.resize(m_num_amr_levels);
    m_b_coeffs.resize(m_num_amr_levels);
    m_cc_mask.resize(m_num_amr_levels);
    m_eb_phi.resize(m_num_amr_levels);
    m_eb_b_coeffs.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        m_b_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        m_cc_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        m_eb_b_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],
                                             m_dmap[amrlev][mglev],
                                             1, 0, MFInfo(), *m_factory[amrlev][mglev]);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev],
                                                    IntVect::TheDimensionVector(idim));
                const int ng = 1;
                m_b_coeffs[amrlev][mglev][idim].define(ba,
                                                       m_dmap[amrlev][mglev],
                                                       ncomp, ng, MFInfo(), *m_factory[amrlev][mglev]);
                m_b_coeffs[amrlev][mglev][idim].setVal(0.0);
            }

            m_cc_mask[amrlev][mglev].define(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1);
            m_cc_mask[amrlev][mglev].setVal(0);
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                const std::vector<IntVect>& pshifts = m_geom[amrlev][mglev].periodicity().shiftIntVect();
                iMultiFab& mask = m_cc_mask[amrlev][mglev];
                const BoxArray& ba = mask.boxArray();
                for (MFIter mfi(mask); mfi.isValid(); ++mfi)
                {
                    Array4<int> const& fab = mask.array(mfi);
                    const Box& bx = mfi.fabbox();
                    for (const auto& iv : pshifts)
                    {
                        ba.intersections(bx+iv, isects);
                        for (const auto& is : isects)
                        {
                            const Box& b = is.second-iv;
                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( b, i, j, k,
                            {
                                fab(i,j,k) = 1;
                            });
                        }
                    }
                }
            }
        }
    }
}

MLEBABecLap::~MLEBABecLap ()
{}

void
MLEBABecLap::setScalars (Real a, Real b)
{
    m_a_scalar = a;
    m_b_scalar = b;
    if (a == 0.0)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            m_a_coeffs[amrlev][0].setVal(0.0);
        }
    }
}

void
MLEBABecLap::setACoeffs (int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
    m_needs_update = true;
}

void
MLEBABecLap::setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
    const int ncomp = getNComp();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (int icomp = 0; icomp < ncomp; ++icomp) {
            MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, icomp, 1, 0);
        }
    }
    m_needs_update = true;
}

void
MLEBABecLap::setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta)
{
    const int ncomp = getNComp();
    if (m_eb_phi[amrlev] == nullptr) {
        const int mglev = 0;
        m_eb_phi[amrlev].reset(new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev],
                                            ncomp, 0, MFInfo(), *m_factory[amrlev][mglev]));
    }
    if (m_eb_b_coeffs[amrlev][0] == nullptr) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            m_eb_b_coeffs[amrlev][mglev].reset(new MultiFab(m_grids[amrlev][mglev],
                                                            m_dmap[amrlev][mglev],
                                                            ncomp, 0, MFInfo(),
                                                            *m_factory[amrlev][mglev]));
        }
    }

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& phiout = m_eb_phi[amrlev]->array(mfi);
        Array4<Real> const& betaout = m_eb_b_coeffs[amrlev][0]->array(mfi);
        FabType t = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (FabType::regular == t or FabType::covered == t) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                phiout(i,j,k,n) = 0.0;
                betaout(i,j,k,n) = 0.0;
            });
        } else {
            Array4<Real const> const& phiin = phi.const_array(mfi);
            Array4<Real const> const& betain = beta.const_array(mfi);
            const auto& flag = flags->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                if (flag(i,j,k).isSingleValued()) {
                    phiout(i,j,k,n) = phiin(i,j,k,n);
                    betaout(i,j,k,n) = betain(i,j,k,0);
                } else {
                    phiout(i,j,k,n) = 0.0;
                    betaout(i,j,k,n) = 0.0;
                }
            });
        }
    }
}

void
MLEBABecLap::setEBHomogDirichlet (int amrlev, const MultiFab& beta)
{
    const int ncomp = getNComp();
    if (m_eb_phi[amrlev] == nullptr) {
        const int mglev = 0;
        m_eb_phi[amrlev].reset(new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev],
                                            ncomp, 0, MFInfo(), *m_factory[amrlev][mglev]));
    }
    if (m_eb_b_coeffs[amrlev][0] == nullptr) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            m_eb_b_coeffs[amrlev][mglev].reset(new MultiFab(m_grids[amrlev][mglev],
                                                            m_dmap[amrlev][mglev],
                                                            ncomp, 0, MFInfo(),
                                                            *m_factory[amrlev][mglev]));
        }
    }

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_eb_phi[amrlev], mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& phifab = m_eb_phi[amrlev]->array(mfi);
        Array4<Real> const& betaout = m_eb_b_coeffs[amrlev][0]->array(mfi);
        FabType t = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            phifab(i,j,k,n) = 0.0;
        });
        if (FabType::regular == t or FabType::covered == t) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                betaout(i,j,k,n) = 0.0;
            });
        } else {
            Array4<Real const> const& betain = beta.const_array(mfi);
            const auto& flag = flags->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                if (flag(i,j,k).isSingleValued()) {
                    betaout(i,j,k,n) = betain(i,j,k,0);
                } else {
                    betaout(i,j,k,n) = 0.0;
                }
            });
        }
    }
}

void
MLEBABecLap::averageDownCoeffs ()
{
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];
        auto& fine_b_coeffs = m_b_coeffs[amrlev];
        
        averageDownCoeffsSameAmrLevel(fine_a_coeffs, fine_b_coeffs,
                                      amrex::GetVecOfPtrs(m_eb_b_coeffs[0]));
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(m_a_coeffs[0], m_b_coeffs[0],
                                  amrex::GetVecOfPtrs(m_eb_b_coeffs[0]));

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_b_coeffs[amrlev][mglev][idim].FillBoundary(m_geom[amrlev][mglev].periodicity());
            }
        }
    }
}

void
MLEBABecLap::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                            Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                            const Vector<MultiFab*>& b_eb)
{
    int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        if (m_a_scalar == 0.0)
        {
            a[mglev].setVal(0.0);
        }
        else
        {
            amrex::EB_average_down(a[mglev-1], a[mglev], 0, 1, mg_coarsen_ratio);
        }

        amrex::EB_average_down_faces(amrex::GetArrOfConstPtrs(b[mglev-1]),
                                     amrex::GetArrOfPtrs(b[mglev]),
                                     mg_coarsen_ratio, 0);

        if (b_eb[mglev])
        {
            amrex::EB_average_down_boundaries(*b_eb[mglev-1], *b_eb[mglev],
                                              mg_coarsen_ratio, 0);
        }
    }
}

void
MLEBABecLap::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    auto const& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto const& fine_b_coeffs = m_b_coeffs[flev  ].back();
    auto      & crse_a_coeffs = m_a_coeffs[flev-1].front();
    auto      & crse_b_coeffs = m_b_coeffs[flev-1].front();
    auto const& fine_eb_b_coeffs = m_eb_b_coeffs[flev  ].back();
    auto      & crse_eb_b_coeffs = m_eb_b_coeffs[flev-1].front();

    if (m_a_scalar != 0.0) {
        amrex::EB_average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, mg_coarsen_ratio);
    }

    amrex::EB_average_down_faces(amrex::GetArrOfConstPtrs(fine_b_coeffs),
                                 amrex::GetArrOfPtrs(crse_b_coeffs),
                                 mg_coarsen_ratio, 0);

    if (fine_eb_b_coeffs) {
        amrex::EB_average_down_boundaries(*fine_eb_b_coeffs, *crse_eb_b_coeffs, mg_coarsen_ratio, 0);
    }
}

void
MLEBABecLap::prepareForSolve ()
{
    BL_PROFILE("MLABecLaplacian::prepareForSolve()");

    MLCellABecLap::prepareForSolve();
    
    averageDownCoeffs();

    if (m_eb_phi[0]) {
        for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev) {
            amrex::EB_average_down_boundaries(*m_eb_phi[amrlev], *m_eb_phi[amrlev-1],
                                              mg_coarsen_ratio, 0);
        }
    }

    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, false);
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == m_lobc[0].end() && ithi == m_hibc[0].end() && !isEBDirichlet())
    {  // No Dirichlet
        for (int alev = 0; alev < m_num_amr_levels; ++alev)
        {
            if (m_domain_covered[alev])
            {
                if (m_a_scalar == 0.0)
                {
                    m_is_singular[alev] = true;
                }
                else
                {
                    Real asum = m_a_coeffs[alev].back().sum();
                    Real amax = m_a_coeffs[alev].back().norm0();
                    m_is_singular[alev] = (asum <= amax * 1.e-12);
                }
            }
        }
    }

    m_needs_update = false;
}

void
MLEBABecLap::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLEBABecLap::Fapply()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];
    
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* barea = (factory) ? &(factory->getBndryArea()) : nullptr;
    const MultiCutFab* bcent = (factory) ? &(factory->getBndryCent()) : nullptr;

    const bool is_eb_dirichlet =  isEBDirichlet();
    const bool is_eb_inhomog = m_is_eb_inhomog;

    const int ncomp = getNComp();

    Array4<Real const> foo;

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& xfab = in.const_array(mfi);
        Array4<Real> const& yfab = out.array(mfi);
        Array4<Real const> const& afab = acoef.const_array(mfi);
        AMREX_D_TERM(Array4<Real const> const& bxfab = bxcoef.const_array(mfi);,
                     Array4<Real const> const& byfab = bycoef.const_array(mfi);,
                     Array4<Real const> const& bzfab = bzcoef.const_array(mfi););

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;

        if (fabtyp == FabType::covered) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
            {
                yfab(i,j,k,n) = 0.0;
            });
        } else if (fabtyp == FabType::regular) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlabeclap_adotx(tbx, yfab, xfab, afab,
                                AMREX_D_DECL(bxfab,byfab,bzfab),
                                dxinvarr, ascalar, bscalar, ncomp);
            });
        } else {
            Array4<int const> const& ccmfab = ccmask.const_array(mfi);
            Array4<EBCellFlag const> const& flagfab = flags->const_array(mfi);
            Array4<Real const> const& vfracfab = vfrac->const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const& apxfab = area[0]->const_array(mfi);,
                         Array4<Real const> const& apyfab = area[1]->const_array(mfi);,
                         Array4<Real const> const& apzfab = area[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const& fcxfab = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcyfab = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fczfab = fcent[2]->const_array(mfi););
            Array4<Real const> const& bafab = barea->const_array(mfi);
            Array4<Real const> const& bcfab = bcent->const_array(mfi);
            Array4<Real const> const& bebfab = (is_eb_dirichlet)
                ? m_eb_b_coeffs[amrlev][mglev]->const_array(mfi) : foo;
            Array4<Real const> const& phiebfab = (is_eb_dirichlet && m_is_eb_inhomog)
                ? m_eb_phi[amrlev]->const_array(mfi) : foo;

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlebabeclap_adotx(tbx, yfab, xfab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                                  ccmfab, flagfab, vfracfab,
                                  AMREX_D_DECL(apxfab,apyfab,apzfab),
                                  AMREX_D_DECL(fcxfab,fcyfab,fczfab), bafab, bcfab, bebfab,
                                  is_eb_dirichlet, phiebfab, is_eb_inhomog, dxinvarr,
                                  ascalar, bscalar, ncomp);
            });
        }
    }
}

void
MLEBABecLap::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLEBABecLap::Fsmooth()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];
    const auto& undrrelxr = m_undrrelxr[amrlev][mglev];
    const auto& maskvals  = m_maskvals [amrlev][mglev];

    OrientationIter oitr;

    const FabSet& f0 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f1 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 1)
    const FabSet& f2 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f3 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 2)
    const FabSet& f4 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f5 = undrrelxr[oitr()]; ++oitr;
#endif
#endif

    const MultiMask& mm0 = maskvals[0];
    const MultiMask& mm1 = maskvals[1];
#if (AMREX_SPACEDIM > 1)
    const MultiMask& mm2 = maskvals[2];
    const MultiMask& mm3 = maskvals[3];
#if (AMREX_SPACEDIM > 2)
    const MultiMask& mm4 = maskvals[4];
    const MultiMask& mm5 = maskvals[5];
#endif
#endif

    const int nc = getNComp();
    const Real* h = m_geom[amrlev][mglev].CellSize();
    AMREX_D_TERM(const Real dhx = m_b_scalar/(h[0]*h[0]);,
                 const Real dhy = m_b_scalar/(h[1]*h[1]);,
                 const Real dhz = m_b_scalar/(h[2]*h[2]));
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Real alpha = m_a_scalar;

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* barea = (factory) ? &(factory->getBndryArea()) : nullptr;
    const MultiCutFab* bcent = (factory) ? &(factory->getBndryCent()) : nullptr;

    bool is_eb_dirichlet =  isEBDirichlet();

    Array4<Real const> foo;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, mfi_info);  mfi.isValid(); ++mfi)
    {
	const auto& m0 = mm0.array(mfi);
        const auto& m1 = mm1.array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& m2 = mm2.array(mfi);
        const auto& m3 = mm3.array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& m4 = mm4.array(mfi);
        const auto& m5 = mm5.array(mfi);
#endif
#endif

        const Box& vbx = mfi.validbox();
        const auto& solnfab = sol.array(mfi);
        const auto& rhsfab  = rhs.const_array(mfi);
        const auto& afab    = acoef.const_array(mfi);

        AMREX_D_TERM(const auto& bxfab = bxcoef.const_array(mfi);,
                     const auto& byfab = bycoef.const_array(mfi);,
                     const auto& bzfab = bzcoef.const_array(mfi););

        const auto& f0fab = f0.const_array(mfi);
        const auto& f1fab = f1.const_array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& f2fab = f2.const_array(mfi);
        const auto& f3fab = f3.const_array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& f4fab = f4.const_array(mfi);
        const auto& f5fab = f5.const_array(mfi);
#endif
#endif

        auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;

        if (fabtyp == FabType::covered)
        {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( vbx, nc, i, j, k, n,
            {
                solnfab(i,j,k,n) = 0.0;
            });
        }
        else if (fabtyp == FabType::regular)
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( vbx, thread_box,
            {
                abec_gsrb(thread_box, solnfab, rhsfab, alpha, afab,
                          AMREX_D_DECL(dhx, dhy, dhz),
                          AMREX_D_DECL(bxfab, byfab, bzfab),
                          AMREX_D_DECL(m0,m2,m4),
                          AMREX_D_DECL(m1,m3,m5),
                          AMREX_D_DECL(f0fab,f2fab,f4fab),
                          AMREX_D_DECL(f1fab,f3fab,f5fab),
                          vbx, redblack, nc);
            });
        }
        else
        {
            Array4<int const> const& ccmfab = ccmask.const_array(mfi);
            Array4<EBCellFlag const> const& flagfab = flags->const_array(mfi);
            Array4<Real const> const& vfracfab = vfrac->const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const& apxfab = area[0]->const_array(mfi);,
                         Array4<Real const> const& apyfab = area[1]->const_array(mfi);,
                         Array4<Real const> const& apzfab = area[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const& fcxfab = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcyfab = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fczfab = fcent[2]->const_array(mfi););
            Array4<Real const> const& bafab = barea->const_array(mfi);
            Array4<Real const> const& bcfab = bcent->const_array(mfi);
            Array4<Real const> const& bebfab = (is_eb_dirichlet)
                ? m_eb_b_coeffs[amrlev][mglev]->const_array(mfi) : foo;

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( vbx, thread_box,
            {
                mlebabeclap_gsrb(thread_box, solnfab, rhsfab, alpha, afab,
                                 AMREX_D_DECL(dhx, dhy, dhz),
                                 AMREX_D_DECL(bxfab,byfab,bzfab),
                                 AMREX_D_DECL(m0,m2,m4),
                                 AMREX_D_DECL(m1,m3,m5),
                                 AMREX_D_DECL(f0fab,f2fab,f4fab),
                                 AMREX_D_DECL(f1fab,f3fab,f5fab),
                                 ccmfab, flagfab, vfracfab,
                                 AMREX_D_DECL(apxfab,apyfab,apzfab),
                                 AMREX_D_DECL(fcxfab,fcyfab,fczfab),
                                 bafab, bcfab, bebfab,
                                 is_eb_dirichlet, vbx, redblack, nc);
            });
        }
    }
}

void
MLEBABecLap::FFlux (int amrlev, const MFIter& mfi, const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                    const FArrayBox& sol, Location loc, const int face_only) const
{
    BL_PROFILE("MLEBABecLap::FFlux()");
    const int at_centroid = (Location::FaceCentroid == loc) ? 1 : 0;
    const int mglev = 0; 
    const Box& box = mfi.tilebox();
    const int ncomp = getNComp();

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get()); 
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr; 

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
    AMREX_D_TERM(const auto& bx = m_b_coeffs[amrlev][mglev][0][mfi];,
                 const auto& by = m_b_coeffs[amrlev][mglev][1][mfi];,
                 const auto& bz = m_b_coeffs[amrlev][mglev][2][mfi];);
    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];

    AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(box,0);,
                 Box const& ybx = amrex::surroundingNodes(box,1);,
                 Box const& zbx = amrex::surroundingNodes(box,2););
    AMREX_D_TERM(Array4<Real> const& fx = flux[0]->array();,
                 Array4<Real> const& fy = flux[1]->array();,
                 Array4<Real> const& fz = flux[2]->array(););

    const auto fabtyp = (flags) ? (*flags)[mfi].getType(box) : FabType::regular; 
    if (fabtyp == FabType::covered) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(xbx, ncomp, i, j, k, n,
        {
            fx(i,j,k,n) = 0.0;
        });
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(ybx, ncomp, i, j, k, n,
        {
            fy(i,j,k,n) = 0.0;
        });
#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(zbx, ncomp, i, j, k, n,
        {
            fz(i,j,k,n) = 0.0;
        });
#endif
    } else if (fabtyp == FabType::regular) {
        MLABecLaplacian::FFlux(box, dxinv, m_b_scalar,
                               Array<FArrayBox const*,AMREX_SPACEDIM>{AMREX_D_DECL(&bx,&by,&bz)},
                               flux, sol, face_only, ncomp);
    } else if (at_centroid) {
        const auto& area = factory->getAreaFrac();
        const auto& fcent = factory->getFaceCent();
        AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
                     Array4<Real const> const& apy = area[1]->const_array(mfi);,
                     Array4<Real const> const& apz = area[2]->const_array(mfi););
        AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                     Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                     Array4<Real const> const& fcz = fcent[2]->const_array(mfi););
        Array4<Real const> const& phi = sol.const_array();
        AMREX_D_TERM(Array4<Real const> const& bxcoef = bx.const_array();,
                     Array4<Real const> const& bycoef = by.const_array();,
                     Array4<Real const> const& bzcoef = bz.const_array(););
        Array4<int const> const& msk = ccmask.const_array(mfi);
        Array4<EBCellFlag const> flg = flags->const_array(mfi);
        AMREX_D_TERM(Real dhx = m_b_scalar*dxinv[0];,
                     Real dhy = m_b_scalar*dxinv[1];,
                     Real dhz = m_b_scalar*dxinv[2];);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
            xbx, txbx,
            {
                mlebabeclap_flux_x(txbx, fx, apx, fcx, phi, bxcoef, msk, flg, dhx, face_only, ncomp, xbx);
            }
            , ybx, tybx,
            {
                mlebabeclap_flux_y(tybx, fy, apy, fcy, phi, bycoef, msk, flg, dhy, face_only, ncomp, ybx);
            }
#if (AMREX_SPACEDIM == 3)
            , zbx, tzbx,
            {
                mlebabeclap_flux_z(tzbx, fz, apz, fcz, phi, bzcoef, msk, flg, dhz, face_only, ncomp, zbx);
            }
#endif
        );
    } else {
        const auto& area = factory->getAreaFrac();
        AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
                     Array4<Real const> const& apy = area[1]->const_array(mfi);,
                     Array4<Real const> const& apz = area[2]->const_array(mfi););
        Array4<Real const> const& phi = sol.const_array();
        AMREX_D_TERM(Array4<Real const> const& bxcoef = bx.const_array();,
                     Array4<Real const> const& bycoef = by.const_array();,
                     Array4<Real const> const& bzcoef = bz.const_array(););
        AMREX_D_TERM(Real dhx = m_b_scalar*dxinv[0];,
                     Real dhy = m_b_scalar*dxinv[1];,
                     Real dhz = m_b_scalar*dxinv[2];);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
            xbx, txbx,
            {
                mlebabeclap_flux_x_0(txbx, fx, apx, phi, bxcoef, dhx, face_only, ncomp, xbx);
            }
            , ybx, tybx,
            {
                mlebabeclap_flux_y_0(tybx, fy, apy, phi, bycoef, dhy, face_only, ncomp, ybx);
            }
#if (AMREX_SPACEDIM == 3)
            , zbx, tzbx,
            {
                mlebabeclap_flux_z_0(tzbx, fz, apz, phi, bzcoef, dhz, face_only, ncomp, zbx);
            }
#endif
        );
    }
}

void
MLEBABecLap::compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                       MultiFab& sol, Location loc) const
{
    BL_PROFILE("MLEBABecLap::compGrad()");

    const int ncomp = getNComp();
    const int at_centroid = (Location::FaceCentroid == loc) ? 1 : 0;
    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[amrlev].get());

    AMREX_D_TERM(const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);,
                 const Real dyi = m_geom[amrlev][mglev].InvCellSize(1);,
                 const Real dzi = m_geom[amrlev][mglev].InvCellSize(2););
    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get()); 
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr; 
    auto area = (factory) ? factory->getAreaFrac() : 
        Array<const MultiCutFab*, AMREX_SPACEDIM>{AMREX_D_DECL(nullptr, nullptr, nullptr)}; 
    auto fcent = (factory) ? factory->getFaceCent():
        Array<const MultiCutFab*, AMREX_SPACEDIM>{AMREX_D_DECL(nullptr, nullptr, nullptr)}; 

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox(); 
        auto fabtyp = (flags) ? (*flags)[mfi].getType(box) :FabType::regular;
        AMREX_D_TERM(Box const& fbx = mfi.nodaltilebox(0);,
                     Box const& fby = mfi.nodaltilebox(1);,
                     Box const& fbz = mfi.nodaltilebox(2););
        AMREX_D_TERM(const auto& gx = grad[0]->array(mfi);,
                     const auto& gy = grad[1]->array(mfi);,
                     const auto& gz = grad[2]->array(mfi););
        const auto& s = sol.const_array(mfi);
        if (fabtyp == FabType::covered) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbx, ncomp, i, j, k, n,
            {
                gx(i,j,k,n) = 0.0;
            });
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fby, ncomp, i, j, k, n,
            {
                gy(i,j,k,n) = 0.0;
            });
#if (AMREX_SPACEDIM == 3)
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbz, ncomp, i, j, k, n,
            {
                gz(i,j,k,n) = 0.0;
            });
#endif
        } else if(fabtyp == FabType::regular) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( fbx, ncomp, i, j, k, n,
            {
                gx(i,j,k,n) = dxi*(s(i,j,k,n) - s(i-1,j,k,n));
            });
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( fby, ncomp, i, j, k, n,
            {
                gy(i,j,k,n) = dyi*(s(i,j,k,n) - s(i,j-1,k,n));
            });
#if (AMREX_SPACEDIM == 3)
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( fbz, ncomp, i, j, k, n,
            {
                gz(i,j,k,n) = dzi*(s(i,j,k,n) - s(i,j,k-1,n));
            });
#endif
        } else if (at_centroid) {
            AMREX_D_TERM(Array4<Real const> const& apx = area[0]->const_array(mfi);,
                         Array4<Real const> const& apy = area[1]->const_array(mfi);,
                         Array4<Real const> const& apz = area[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););
            Array4<int const> const& msk = ccmask.const_array(mfi);
            Array4<EBCellFlag const> flg = flags->const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                fbx, txbx,
                {
                    mlebabeclap_grad_x(txbx, gx, s, apx, fcx, msk, flg, dxi, ncomp);
                }
                , fby, tybx,
                {
                    mlebabeclap_grad_y(tybx, gy, s, apy, fcy, msk, flg, dyi, ncomp);
                }
#if (AMREX_SPACEDIM == 3)
                , fbz, tzbx,
                {
                    mlebabeclap_grad_z(tzbx, gz, s, apz, fcz, msk, flg, dzi, ncomp);
                }
#endif
            );
        } else {
            AMREX_D_TERM(Array4<Real const> const& ax = area[0]->const_array(mfi);,
                         Array4<Real const> const& ay = area[1]->const_array(mfi);,
                         Array4<Real const> const& az = area[2]->const_array(mfi););
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                fbx, txbx,
                {
                    mlebabeclap_grad_x_0(txbx, gx, s, ax, dxi, ncomp);
                }
                , fby, tybx,
                {
                    mlebabeclap_grad_y_0(tybx, gy, s, ay, dyi, ncomp);
                }
#if (AMREX_SPACEDIM == 3)
                , fbz, tzbx,
                {
                    mlebabeclap_grad_z_0(tzbx, gz, s, az, dzi, ncomp);
                }
#endif
            );
        }
    }
}

void
MLEBABecLap::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
    const iMultiFab& ccmask = m_cc_mask[amrlev][mglev];

    const auto dxinvarray = m_geom[amrlev][mglev].InvCellSizeArray();
    AMREX_D_TERM(Real dhx = m_b_scalar*dxinvarray[0]*dxinvarray[0];,
                 Real dhy = m_b_scalar*dxinvarray[1]*dxinvarray[1];,
                 Real dhz = m_b_scalar*dxinvarray[2]*dxinvarray[2];);

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* barea = (factory) ? &(factory->getBndryArea()) : nullptr;
    const MultiCutFab* bcent = (factory) ? &(factory->getBndryCent()) : nullptr;

    bool is_eb_dirichlet =  isEBDirichlet();

    Array4<Real const> foo;

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;
    const int ncomp = getNComp();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& fab = mf.array(mfi);
        Array4<Real const> const& afab = acoef.const_array(mfi);
        AMREX_D_TERM(Array4<Real const> const& bxfab = bxcoef.const_array(mfi);,
                     Array4<Real const> const& byfab = bycoef.const_array(mfi);,
                     Array4<Real const> const& bzfab = bzcoef.const_array(mfi););

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;

        if (fabtyp == FabType::regular)
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
            {
                mlabeclap_normalize(bx, fab, afab, AMREX_D_DECL(bxfab, byfab, bzfab),
                                    dxinvarray, ascalar, bscalar, ncomp);
            });
        }
        else if (fabtyp == FabType::singlevalued)
        {
            Array4<Real const> const& bebfab
                = (is_eb_dirichlet) ? m_eb_b_coeffs[amrlev][mglev]->const_array(mfi) : foo;
            Array4<int const> const& ccmfab = ccmask.const_array(mfi);
            Array4<EBCellFlag const> const& flagfab = flags->const_array(mfi);
            Array4<Real const> const& vfracfab = vfrac->const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const& apxfab = area[0]->const_array(mfi);,
                         Array4<Real const> const& apyfab = area[1]->const_array(mfi);,
                         Array4<Real const> const& apzfab = area[2]->const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const& fcxfab = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcyfab = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fczfab = fcent[2]->const_array(mfi););
            Array4<Real const> const& bafab = barea->const_array(mfi);
            Array4<Real const> const& bcfab = bcent->const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlebabeclap_normalize(bx, fab, ascalar, afab,
                                      AMREX_D_DECL(dhx, dhy, dhz),
                                      AMREX_D_DECL(bxfab, byfab, bzfab),
                                      ccmfab, flagfab, vfracfab,
                                      AMREX_D_DECL(apxfab,apyfab,apzfab),
                                      AMREX_D_DECL(fcxfab,fcyfab,fczfab),
                                      bafab, bcfab, bebfab, is_eb_dirichlet, ncomp);
            });
        }
    }
}

void
MLEBABecLap::restriction (int, int, MultiFab& crse, MultiFab& fine) const
{
    const int ncomp = getNComp();
    amrex::EB_average_down(fine, crse, 0, ncomp, 2);
}

void
MLEBABecLap::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLEBABecLap::interpolation()");

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][fmglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;

    const int ncomp = getNComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;

        Array4<Real const> const& cfab = crse.const_array(mfi);
        Array4<Real> const& ffab = fine.array(mfi);

        if (fabtyp == FabType::regular)
        {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                int ic = amrex::coarsen(i,2);
                int jc = amrex::coarsen(j,2);
                int kc = amrex::coarsen(k,2);
                ffab(i,j,k,n) += cfab(ic,jc,kc,n);
            });
        }
        else if (fabtyp == FabType::singlevalued)
        {
            Array4<EBCellFlag const> const& flg = flags->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                if (!flg(i,j,k).isCovered()) {
                    int ic = amrex::coarsen(i,2);
                    int jc = amrex::coarsen(j,2);
                    int kc = amrex::coarsen(k,2);
                    ffab(i,j,k,n) += cfab(ic,jc,kc,n);
                }
            });
        }
    }
}

void
MLEBABecLap::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                     const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto amrrr = AMRRefRatio(camrlev);
    const int ncomp = getNComp();
    amrex::EB_average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
    amrex::EB_average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

void
MLEBABecLap::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                      const MLMGBndry* bndry, bool skip_fillboundary) const
{
    BL_PROFILE("MLEBABecLap::applyBC()");

    // No coarsened boundary values, cannot apply inhomog at mglev>0.
    BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
    BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

    const int ncomp = getNComp();
    if (!skip_fillboundary) {
        const int cross = false;
        in.FillBoundary(0, ncomp, m_geom[amrlev][mglev].periodicity(),cross);
    }

    int m_is_inhomog = bc_mode == BCMode::Inhomogeneous;
    int flagbc = m_is_inhomog;
    m_is_eb_inhomog = s_mode == StateMode::Solution;
    const int imaxorder = maxorder;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(imaxorder <= 4, "MLEBABecLap::applyBC: maxorder too high");

    const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);
    const Real dyi = (AMREX_SPACEDIM >= 2) ? m_geom[amrlev][mglev].InvCellSize(1) : 1.0;
    const Real dzi = (AMREX_SPACEDIM == 3) ? m_geom[amrlev][mglev].InvCellSize(2) : 1.0;

    const auto& maskvals = m_maskvals[amrlev][mglev];
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    const auto& ccmask = m_cc_mask[amrlev][mglev];
    
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    
    FArrayBox foofab(Box::TheUnitBox(),ncomp);
    const auto& foo = foofab.array();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(in, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx   = mfi.validbox();
        const auto& iofab = in.array(mfi);

        auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;
        if (fabtyp != FabType::covered)
        {
            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const Orientation olo(idim,Orientation::low);
                const Orientation ohi(idim,Orientation::high);
                Box blo = amrex::adjCellLo(vbx, idim);
                Box bhi = amrex::adjCellHi(vbx, idim);
                if (fabtyp != FabType::regular) {
                    blo.grow(IntVect(1)-IntVect::TheDimensionVector(idim));
                    bhi.grow(IntVect(1)-IntVect::TheDimensionVector(idim));
                }
                const int blen = vbx.length(idim);
                const auto& mlo = maskvals[olo].array(mfi);
                const auto& mhi = maskvals[ohi].array(mfi);
                const auto& bvlo = (bndry != nullptr) ? bndry->bndryValues(olo).array(mfi) : foo;
                const auto& bvhi = (bndry != nullptr) ? bndry->bndryValues(ohi).array(mfi) : foo;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    const BoundCond bctlo = bdcv[icomp][olo];
                    const BoundCond bcthi = bdcv[icomp][ohi];
                    const Real bcllo = bdlv[icomp][olo];
                    const Real bclhi = bdlv[icomp][ohi];
                    if (fabtyp == FabType::regular)
                    {
                        if (idim == 0) {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mllinop_apply_bc_x(0, tboxlo, blen, iofab, mlo,
                                               bctlo, bcllo, bvlo,
                                               imaxorder, dxi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mllinop_apply_bc_x(1, tboxhi, blen, iofab, mhi,
                                               bcthi, bclhi, bvhi,
                                               imaxorder, dxi, flagbc, icomp);
                            });
                        } else if (idim == 1) {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mllinop_apply_bc_y(0, tboxlo, blen, iofab, mlo,
                                               bctlo, bcllo, bvlo,
                                               imaxorder, dyi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mllinop_apply_bc_y(1, tboxhi, blen, iofab, mhi,
                                               bcthi, bclhi, bvhi,
                                               imaxorder, dyi, flagbc, icomp);
                            });
                        } else {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mllinop_apply_bc_z(0, tboxlo, blen, iofab, mlo,
                                               bctlo, bcllo, bvlo,
                                               imaxorder, dzi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mllinop_apply_bc_z(1, tboxhi, blen, iofab, mhi,
                                               bcthi, bclhi, bvhi,
                                               imaxorder, dzi, flagbc, icomp);
                            });
                        }
                    }
                    else // irregular
                    {
                        const auto& ap = area[idim]->const_array(mfi);
                        const auto& mask = ccmask.const_array(mfi);
                        if (idim == 0) {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mlebabeclap_apply_bc_x(0, tboxlo, blen, iofab, mask, ap,
                                                   bctlo, bcllo, bvlo,
                                                   imaxorder, dxi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mlebabeclap_apply_bc_x(1, tboxhi, blen, iofab, mask, ap,
                                                   bcthi, bclhi, bvhi,
                                                   imaxorder, dxi, flagbc, icomp);
                            });
                        } else if (idim == 1) {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mlebabeclap_apply_bc_y(0, tboxlo, blen, iofab, mask, ap,
                                                   bctlo, bcllo, bvlo,
                                                   imaxorder, dyi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mlebabeclap_apply_bc_y(1, tboxhi, blen, iofab, mask, ap,
                                                   bcthi, bclhi, bvhi,
                                                   imaxorder, dyi, flagbc, icomp);
                            });
                        } else {
                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                            blo, tboxlo, {
                            mlebabeclap_apply_bc_z(0, tboxlo, blen, iofab, mask, ap,
                                                   bctlo, bcllo, bvlo,
                                                   imaxorder, dzi, flagbc, icomp);
                            },
                            bhi, tboxhi, {
                            mlebabeclap_apply_bc_z(1, tboxhi, blen, iofab, mask, ap,
                                                   bcthi, bclhi, bvhi,
                                                   imaxorder, dzi, flagbc, icomp);
                            });
                        }
                    }
                }
            }
        }
    }
}

void
MLEBABecLap::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    StateMode s_mode, const MLMGBndry* bndry) const
{
    BL_PROFILE("MLEBABecLap::apply()");
    applyBC(amrlev, mglev, in, bc_mode, s_mode, bndry);
    Fapply(amrlev, mglev, out, in);
}

void
MLEBABecLap::update ()
{
    if (MLCellABecLap::needsUpdate()) MLCellABecLap::update();

    averageDownCoeffs();

    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, false);
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == m_lobc[0].end() && ithi == m_hibc[0].end() && !isEBDirichlet())
    {  // No Dirichlet
        for (int alev = 0; alev < m_num_amr_levels; ++alev)
        {
            if (m_domain_covered[alev])
            {
                if (m_a_scalar == 0.0)
                {
                    m_is_singular[alev] = true;
                }
                else
                {
                    Real asum = m_a_coeffs[alev].back().sum();
                    Real amax = m_a_coeffs[alev].back().norm0();
                    m_is_singular[alev] = (asum <= amax * 1.e-12);
                }
            }
        }
    }

    m_needs_update = false;
}

#ifdef AMREX_USE_HYPRE
std::unique_ptr<Hypre>
MLEBABecLap::makeHypre (Hypre::Interface hypre_interface) const
{
    auto hypre_solver = MLCellABecLap::makeHypre(hypre_interface);
    auto ijmatrix_solver = dynamic_cast<HypreABecLap3*>(hypre_solver.get());
    ijmatrix_solver->setEBDirichlet(m_eb_b_coeffs[0].back().get());
    return hypre_solver;
}
#endif

#ifdef AMREX_USE_PETSC
std::unique_ptr<PETScABecLap>
MLEBABecLap::makePETSc () const
{
    auto petsc_solver = MLCellABecLap::makePETSc();
    petsc_solver->setEBDirichlet(m_eb_b_coeffs[0].back().get());
    return petsc_solver;
}
#endif


}
