#include <AMReX_MLNodeCrossLaplacian.H>

#if (AMREX_SPACEDIM == 2)

#include <AMReX_MLNodeCrossLap_K.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_Scan.H>

namespace amrex {

namespace {
    // Stencil coefficients of D sigma G: corners, x neighbors (y neighbors
    // have the opposite sign) and center.
    struct CrossCoef { Real cc, cx, c0; };
    CrossCoef cross_coef (Geometry const& geom, Real sigma)
    {
        auto const& dxinv = geom.InvCellSizeArray();
        Real const ax = sigma*dxinv[0]*dxinv[0];
        Real const ay = sigma*dxinv[1]*dxinv[1];
        return {.cc = Real(0.25)*(ax+ay), .cx = Real(0.5)*(ax-ay), .c0 = -(ax+ay)};
    }
}

MLNodeCrossLaplacian::MLNodeCrossLaplacian (const Vector<Geometry>& a_geom,
                                            const Vector<BoxArray>& a_grids,
                                            const Vector<DistributionMapping>& a_dmap,
                                            const LPInfo& a_info,
                                            Real a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_const_sigma);
}

void
MLNodeCrossLaplacian::define (const Vector<Geometry>& a_geom,
                              const Vector<BoxArray>& a_grids,
                              const Vector<DistributionMapping>& a_dmap,
                              const LPInfo& a_info,
                              Real a_const_sigma)
{
    BL_PROFILE("MLNodeCrossLaplacian::define()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_geom.size() == 1,
        "MLNodeCrossLaplacian: only a single AMR level is supported");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!a_geom[0].IsRZ(),
        "MLNodeCrossLaplacian: RZ is not supported");

    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    m_const_sigma = a_const_sigma;

    // Geometric coarsening does not work for this operator.
    LPInfo lpinfo = a_info;
    lpinfo.setMaxCoarseningLevel(0);

    m_coarsening_strategy = CoarseningStrategy::Sigma; // fills nodes outside Neumann BC
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, lpinfo);
}

void
MLNodeCrossLaplacian::restriction (int, int, MultiFab&, MultiFab&) const
{
    amrex::Abort("MLNodeCrossLaplacian::restriction: single MG level only");
}

void
MLNodeCrossLaplacian::interpolation (int, int, MultiFab&, const MultiFab&) const
{
    amrex::Abort("MLNodeCrossLaplacian::interpolation: single MG level only");
}

void
MLNodeCrossLaplacian::averageDownSolutionRHS (int, MultiFab&, MultiFab&,
                                              const MultiFab&, const MultiFab&)
{
    amrex::Abort("MLNodeCrossLaplacian::averageDownSolutionRHS: single AMR level only");
}

void
MLNodeCrossLaplacian::reflux (int, MultiFab&, const MultiFab&, const MultiFab&,
                              MultiFab&, MultiFab&, const MultiFab&) const
{
    amrex::Abort("MLNodeCrossLaplacian::reflux: single AMR level only");
}

void
MLNodeCrossLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeCrossLaplacian::prepareForSolve()");
    MLNodeLinOp::prepareForSolve();
    buildMasks();
}

void
MLNodeCrossLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel)
{
    BL_PROFILE("MLNodeCrossLaplacian::compRHS()");

    AMREX_ALWAYS_ASSERT(rhs.size() == 1 && vel.size() == 1);

    if (!m_masks_built) { buildMasks(); }

    const Geometry& geom = m_geom[0][0];
    const auto lobc = LoBC();
    const auto hibc = HiBC();
    const Box& nddom = amrex::surroundingNodes(geom.Domain());
    const auto dxinvarr = geom.InvCellSizeArray();

    AMREX_ASSERT(vel[0]->nComp() >= AMREX_SPACEDIM);
    AMREX_ASSERT(vel[0]->nGrowVect().allGE(1));

    bool has_inflow = false;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        has_inflow = has_inflow || (lobc[idim] == LinOpBCType::inflow ||
                                    hibc[idim] == LinOpBCType::inflow);
    }

    if (has_inflow) { // Zero out transverse velocity so that it's not seen.
        Box domain = geom.Domain();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (lobc[idim] != LinOpBCType::inflow) {
                domain.growLo(idim,1);
            }
            if (hibc[idim] != LinOpBCType::inflow) {
                domain.growHi(idim,1);
            }
        }
        const auto dlo = domain.smallEnd();
        const auto dhi = domain.bigEnd();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(1);
            Array4<Real> const& vfab = vel[0]->array(mfi);
            if ( ! domain.contains(bx) ) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    IntVect cell(AMREX_D_DECL(i,j,k));
                    for (int in = 0; in < AMREX_SPACEDIM; ++in) {
                        for (int it = 0; it < AMREX_SPACEDIM; ++it) {
                            if (it != in) {
                                if (cell[in] < dlo[in] || cell[in] > dhi[in]) {
                                    vfab(i,j,k,it) = Real(0.0);
                                }
                            }
                        }
                    }
                });
            }
        }
    }

    vel[0]->FillBoundary(0, AMREX_SPACEDIM, IntVect(1), geom.periodicity());

    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& rhsarr = rhs[0]->array(mfi);
        Array4<Real const> const& velarr = vel[0]->const_array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
        {
            mlndlap_divu(i,j,k,rhsarr,velarr,dmskarr,dxinvarr,nddom,lobc,hibc,false);
        });
        mlndlap_impose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);
    }
}

void
MLNodeCrossLaplacian::updateVelocity (const Vector<MultiFab*>& vel,
                                      const Vector<MultiFab const*>& sol) const
{
    AMREX_ALWAYS_ASSERT(vel.size() == 1 && sol.size() == 1);
    const auto dxinv = m_geom[0][0].InvCellSizeArray();
    const Real sigma = m_const_sigma;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& varr = vel[0]->array(mfi);
        Array4<Real const> const& solarr = sol[0]->const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            mlndlap_mknewu_c(i,j,k,varr,solarr,sigma,dxinv,false);
        });
    }
}

void
MLNodeCrossLaplacian::compGradPhi (const Vector<MultiFab*>& grad,
                                   const Vector<MultiFab const*>& sol) const
{
    AMREX_ALWAYS_ASSERT(grad.size() == 1 && sol.size() == 1);
    const auto dxinv = m_geom[0][0].InvCellSizeArray();
    const Real facx = Real(0.5)*dxinv[0];
    const Real facy = Real(0.5)*dxinv[1];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*grad[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& g = grad[0]->array(mfi);
        Array4<Real const> const& p = sol[0]->const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            g(i,j,k,0) = facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
            g(i,j,k,1) = facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
        });
    }
}

void
MLNodeCrossLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLNodeCrossLaplacian::Fapply()");

    auto const coef = cross_coef(m_geom[amrlev][mglev], m_const_sigma);
    Real const cc = coef.cc, cx = coef.cx, c0 = coef.c0;

    auto const& in_a = in.const_arrays();
    auto const& out_a = out.arrays();
    auto const& dmsk_a = m_dirichlet_mask[amrlev][mglev]->const_arrays();

    amrex::ParallelFor(out,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        if (dmsk_a[box_no](i,j,k)) {
            out_a[box_no](i,j,k) = Real(0.0);
        } else {
            out_a[box_no](i,j,k) = mlndcrosslap_adotx(i, j, k, in_a[box_no], cc, cx, c0);
        }
    });
    if (!Gpu::inNoSyncRegion()) {
        Gpu::streamSynchronize();
    }
}

void
MLNodeCrossLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeCrossLaplacian::Fsmooth()");

    auto const coef = cross_coef(m_geom[amrlev][mglev], m_const_sigma);
    Real const cc = coef.cc, cx = coef.cx, c0 = coef.c0;

    auto const& dmsk = *m_dirichlet_mask[amrlev][mglev];

    // Nodes with the same index parities are not coupled by the stencil.
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& sol_a = sol.arrays();
        auto const& rhs_a = rhs.const_arrays();
        auto const& dmsk_a = dmsk.const_arrays();
        for (int color = 0; color < 4; ++color) {
            ParallelForStrided(sol, IntVect(2), multicolor_offset(color),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                mlndcrosslap_gauss_seidel(i, j, k, sol_a[box_no], rhs_a[box_no], dmsk_a[box_no],
                                          cc, cx, c0);
            });
        }
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
        for (int color = 0; color < 4; ++color) {
            IntVect const offset = multicolor_offset(color);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(sol,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& sol_a = sol.array(mfi);
                auto const& rhs_a = rhs.const_array(mfi);
                auto const& dmsk_a = dmsk.const_array(mfi);
                auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                lo.x += (offset[0] - lo.x) & 1;
                lo.y += (offset[1] - lo.y) & 1;
                for (int j = lo.y; j <= hi.y; j += 2) {
                for (int i = lo.x; i <= hi.x; i += 2) {
                    mlndcrosslap_gauss_seidel(i, j, 0, sol_a, rhs_a, dmsk_a, cc, cx, c0);
                }}
            }
        }
    }
}

void
MLNodeCrossLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    amrex::ignore_unused(amrlev,mglev,mf);
}

bool
MLNodeCrossLaplacian::checkerboardIsNull () const
{
    if (!isSingular(0)) { return false; }
    const Geometry& geom = m_geom[0][0];
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim) && (geom.Domain().length(idim) % 2 != 0)) { return false; }
    }
    return true;
}

Vector<Real>
MLNodeCrossLaplacian::getSolvabilityOffset (int amrlev, int mglev, MultiFab const& rhs) const
{
    if (!checkerboardIsNull()) {
        return MLNodeLinOp::getSolvabilityOffset(amrlev, mglev, rhs);
    }

    // Weighted means of rhs and of checkerboard * rhs. The two null vectors
    // are orthogonal under these weights.
    const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
    const auto& mask_ma = mask.const_arrays();
    const auto& rhs_ma = rhs.const_arrays();
    auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{},
                       TypeList<Real,Real,Real>{},
                       rhs, IntVect(0),
                       [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                           -> GpuTuple<Real,Real,Real>
                       {
                           Real const w = mask_ma[box_no](i,j,k);
                           Real const chi = ((i+j) & 1) ? Real(-1.0) : Real(1.0);
                           return { w * rhs_ma[box_no](i,j,k),
                                    w * chi * rhs_ma[box_no](i,j,k), w };
                       });
    Real s1 = amrex::get<0>(r);
    Real s2 = amrex::get<1>(r);
    Real s3 = amrex::get<2>(r);
    ParallelAllReduce::Sum<Real>({s1,s2,s3}, ParallelContext::CommunicatorSub());
    return {s1/s3, s2/s3};
}

void
MLNodeCrossLaplacian::fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                              Vector<Real> const& offset) const
{
    if (offset.size() < 2) {
        MLNodeLinOp::fixSolvabilityByOffset(amrlev, mglev, rhs, offset);
        return;
    }
    Real const c0 = offset[0];
    Real const c1 = offset[1];
    auto const& rhs_ma = rhs.arrays();
    ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        rhs_ma[box_no](i,j,k) -= ((i+j) & 1) ? c0 - c1 : c0 + c1;
    });
    Gpu::streamSynchronize();
}

void
MLNodeCrossLaplacian::postSolve (Vector<MultiFab*> const& sol) const
{
    // Also in preconditioner mode: GMRESMLMG builds phi from these results.
    if (!checkerboardIsNull()) { return; }

    BL_PROFILE("MLNodeCrossLaplacian::postSolve()");

    MultiFab& phi = *sol[0];
    auto const& phi_a = phi.const_arrays();
    auto const& omsk_a = m_owner_mask_top->const_arrays();
    auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                       phi, IntVect(0),
                       [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                           -> GpuTuple<Real,Real>
                       {
                           Real const m = omsk_a[box_no](i,j,0) ? Real(1.0) : Real(0.0);
                           Real const chi = ((i+j) & 1) ? Real(-1.0) : Real(1.0);
                           return { m * chi * phi_a[box_no](i,j,0), m };
                       });
    Real s1 = amrex::get<0>(r);
    Real s2 = amrex::get<1>(r);
    ParallelAllReduce::Sum<Real>({s1,s2}, ParallelContext::CommunicatorSub());
    Real const c = s1/s2;

    auto const& phi_ma = phi.arrays();
    ParallelFor(phi, phi.nGrowVect(), [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
    {
        phi_ma[box_no](i,j,0) -= ((i+j) & 1) ? -c : c;
    });
    Gpu::streamSynchronize();
}

void
MLNodeCrossLaplacian::fillAlgMatrix (int mglev, MFIter const& mfi,
                                     Array4<Long const> const& gid,
                                     Array4<int const> const& lid,
                                     Long* ncols, Long* cols, Real* mat) const
{
    const int amrlev = 0;
    // The AlgMG matrix keeps only the corner couplings. With dx == dy it is
    // the operator itself; otherwise it is a preconditioner for it.
    Real const cc = cross_coef(m_geom[amrlev][mglev], m_const_sigma).cc;
    const Box& ccdom = m_geom[amrlev][mglev].growPeriodicDomain(1);
    const Box& ndbx = mfi.validbox();
    constexpr int nslots = 5;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE
            (static_cast<Long>(ndbx.numPts())*nslots <
             static_cast<Long>(std::numeric_limits<int>::max()),
             "The Box is too big.  We could use Long here, but it would much slower.");
        const int nmax = static_cast<int>(ndbx.numPts()) * nslots;
        const auto ndlo = amrex::lbound(ndbx);
        const int lenx = ndbx.length(0);
        amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 int const inode = offset / nslots;
                 int const i = ndlo.x + inode % lenx;
                 int const j = ndlo.y + inode / lenx;
                 Long col = 0; Real val = 0;
                 return (lid(i,j,0) >= 0 &&
                         mlndcrosslap_matrix_entry(i, j, 0, offset - inode*nslots, gid, ccdom,
                                                   cc, col, val));
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 int const inode = offset / nslots;
                 int const i = ndlo.x + inode % lenx;
                 int const j = ndlo.y + inode / lenx;
                 int const slot = offset - inode*nslots;
                 if (lid(i,j,0) < 0) { return; }
                 Long col = 0; Real val = 0;
                 if (mlndcrosslap_matrix_entry(i, j, 0, slot, gid, ccdom, cc, col, val)) {
                     cols[ps] = col;
                     mat[ps] = val;
                 }
                 if (slot == 0) {
                     Long n = 0;
                     for (int s = 0; s < nslots; ++s) {
                         if (mlndcrosslap_matrix_entry(i, j, 0, s, gid, ccdom, cc, col, val)) { ++n; }
                     }
                     ncols[lid(i,j,0)] = n;
                 }
             },
             amrex::Scan::Type::exclusive);
    } else
#endif
    {
        Long nelems = 0;
        amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
        {
            if (lid(i,j,k) >= 0) {
                Long const nelems_old = nelems;
                for (int s = 0; s < nslots; ++s) {
                    Long col = 0; Real val = 0;
                    if (mlndcrosslap_matrix_entry(i, j, k, s, gid, ccdom, cc, col, val)) {
                        cols[nelems] = col;
                        mat[nelems] = val;
                        ++nelems;
                    }
                }
                ncols[lid(i,j,k)] = nelems - nelems_old;
            }
        });
    }
}

void
MLNodeCrossLaplacian::fillRHS (int mglev, MFIter const& mfi, Array4<int const> const& lid,
                               Real* rhs, Array4<Real const> const& bfab) const
{
    // The matrix rows at Neumann boundaries are half stencils, so undo the
    // doubling of the right-hand side there.
    const Box& nddom = amrex::surroundingNodes(Geom(0,mglev).Domain());
    const Box& bx = mfi.validbox();
    const auto lobc = LoBC();
    const auto hibc = HiBC();
    GpuArray<int,AMREX_SPACEDIM> neumann_lo{std::numeric_limits<int>::lowest(),
                                            std::numeric_limits<int>::lowest()};
    GpuArray<int,AMREX_SPACEDIM> neumann_hi{std::numeric_limits<int>::lowest(),
                                            std::numeric_limits<int>::lowest()};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (lobc[idim] == LinOpBCType::Neumann || lobc[idim] == LinOpBCType::inflow) {
            if (bx.smallEnd(idim) == nddom.smallEnd(idim)) {
                neumann_lo[idim] = bx.smallEnd(idim);
            }
        }
        if (hibc[idim] == LinOpBCType::Neumann || hibc[idim] == LinOpBCType::inflow) {
            if (bx.bigEnd(idim) == nddom.bigEnd(idim)) {
                neumann_hi[idim] = bx.bigEnd(idim);
            }
        }
    }

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
    {
        if (lid(i,j,k) >= 0) {
            Real fac = Real(1.0);
            if ((neumann_lo[0] == i) || neumann_hi[0] == i) { fac *= Real(0.5); }
            if ((neumann_lo[1] == j) || neumann_hi[1] == j) { fac *= Real(0.5); }
            rhs[lid(i,j,k)] = fac * bfab(i,j,k);
        }
    });
}

}

#endif
