#include <AMReX_MLCurlCurl.H>

namespace amrex {

MLCurlCurl::MLCurlCurl (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

void MLCurlCurl::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info)
{
    MLLinOpT<MF>::define(a_geom, a_grids, a_dmap, a_info, {});

    m_dotmask.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_dotmask[amrlev].resize(this->m_num_mg_levels[amrlev]);
    }

    m_lusolver.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_lusolver[amrlev].resize(this->m_num_mg_levels[amrlev]);
    }
}

void MLCurlCurl::setScalars (RT a_alpha, RT a_beta) noexcept
{
    m_alpha = a_alpha;
    m_beta = a_beta;
    AMREX_ASSERT(m_beta > RT(0));
}

void MLCurlCurl::prepareRHS (Vector<MF*> const& rhs) const
{
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (auto& mf : *rhs[amrlev]) {
            mf.OverrideSync(m_geom[amrlev][0].periodicity());
        }
    }
}

void MLCurlCurl::setDirichletNodesToZero (int amrlev, int mglev, MF& a_mf) const
{
    MFItInfo mfi_info{};
#ifdef AMREX_USE_GPU
    Vector<Array4BoxTag<RT>> tags;
    mfi_info.DisableDeviceSync();
#endif

    for (auto& mf : a_mf)
    {
        auto const idxtype = mf.ixType();
        Box const domain = amrex::convert(m_geom[amrlev][mglev].Domain(), idxtype);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,mfi_info); mfi.isValid(); ++mfi) {
            auto const& vbx = mfi.validbox();
            auto const& a = mf.array(mfi);
            for (OrientationIter oit; oit; ++oit) {
                Orientation const face = oit();
                int const idim = face.coordDir();
                bool is_dirichlet = face.isLow()
                    ? m_lobc[0][idim] == LinOpBCType::Dirichlet
                    : m_hibc[0][idim] == LinOpBCType::Dirichlet;
                if (is_dirichlet && domain[face] == vbx[face] &&
                    idxtype.nodeCentered(idim))
                {
                    Box b = vbx;
                    b.setRange(idim, vbx[face], 1);
#ifdef AMREX_USE_GPU
                    tags.emplace_back(Array4BoxTag<RT>{a,b});
#else
                    amrex::LoopOnCpu(b, [&] (int i, int j, int k)
                    {
                        a(i,j,k) = RT(0.0);
                    });
#endif
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    ParallelFor(tags,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, Array4BoxTag<RT> const& tag) noexcept
    {
        tag.dfab(i,j,k) = RT(0.0);
    });
#endif
}

void MLCurlCurl::setLevelBC (int amrlev, const MF* levelbcdata, // TODO
                             const MF* robinbc_a, const MF* robinbc_b,
                             const MF* robinbc_f)
{
    amrex::ignore_unused(amrlev, levelbcdata, robinbc_a, robinbc_b, robinbc_f);
}

void MLCurlCurl::restriction (int amrlev, int cmglev, MF& crse, MF& fine) const
{
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[cmglev-1];
    AMREX_ALWAYS_ASSERT(ratio == 2);

    applyBC(amrlev, cmglev-1, fine, CurlCurlStateType::r);

    auto dinfo = getDirichletInfo(amrlev,cmglev-1);

    for (int idim = 0; idim < 3; ++idim) {
        bool need_parallel_copy = !amrex::isMFIterSafe(crse[idim], fine[idim]);
        MultiFab cfine;
        if (need_parallel_copy) {
            BoxArray const& ba = amrex::coarsen(fine[idim].boxArray(), 2);
            cfine.define(ba, fine[idim].DistributionMap(), 1, 0);
        }

        MultiFab* pcrse = (need_parallel_copy) ? &cfine : &(crse[idim]);

        auto const& crsema = pcrse->arrays();
        auto const& finema = fine[idim].const_arrays();
        ParallelFor(*pcrse, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
        {
            mlcurlcurl_restriction(idim,i,j,k,crsema[bno],finema[bno],dinfo);
        });
        Gpu::streamSynchronize();

        if (need_parallel_copy) {
            crse[idim].ParallelCopy(cfine);
        }
    }
}

void MLCurlCurl::interpolation (int amrlev, int fmglev, MF& fine,
                                const MF& crse) const
{
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];
    AMREX_ALWAYS_ASSERT(ratio == 2);

    auto dinfo = getDirichletInfo(amrlev,fmglev);

    for (int idim = 0; idim < 3; ++idim) {
        bool need_parallel_copy = !amrex::isMFIterSafe(crse[idim], fine[idim]);
        MultiFab cfine;
        MultiFab const* cmf = &(crse[idim]);
        if (need_parallel_copy) {
            BoxArray const& ba = amrex::coarsen(fine[idim].boxArray(), 2);
            cfine.define(ba, fine[idim].DistributionMap(), 1, 0);
            cfine.ParallelCopy(crse[idim]);
            cmf = &cfine;
        }
        auto const& finema = fine[idim].arrays();
        auto const& crsema = cmf->const_arrays();
        ParallelFor(fine[idim], [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
        {
            if (!dinfo.is_dirichlet_edge(idim,i,j,k)) {
                mlcurlcurl_interpadd(idim,i,j,k,finema[bno],crsema[bno]);
            }
        });
    }
    Gpu::streamSynchronize();
}

void
MLCurlCurl::apply (int amrlev, int mglev, MF& out, MF& in, BCMode /*bc_mode*/,
                   StateMode /*s_mode*/, const MLMGBndryT<MF>* /*bndry*/) const
{
    applyBC(amrlev, mglev, in, CurlCurlStateType::x);

    auto adxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        adxinv[idim] *= std::sqrt(m_alpha);
    }
    auto const b = m_beta;

    auto dinfo = getDirichletInfo(amrlev,mglev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& xbx = mfi.tilebox(out[0].ixType().toIntVect());
        Box const& ybx = mfi.tilebox(out[1].ixType().toIntVect());
        Box const& zbx = mfi.tilebox(out[2].ixType().toIntVect());
        auto const& xout = out[0].array(mfi);
        auto const& yout = out[1].array(mfi);
        auto const& zout = out[2].array(mfi);
        auto const& xin = in[0].array(mfi);
        auto const& yin = in[1].array(mfi);
        auto const& zin = in[2].array(mfi);
        amrex::ParallelFor(xbx, ybx, zbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_x_edge(i,j,k)) {
                xout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_x(i,j,k,xout,xin,yin,zin,b,adxinv);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_y_edge(i,j,k)) {
                yout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_y(i,j,k,yout,xin,yin,zin,b,adxinv);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_z_edge(i,j,k)) {
                zout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_z(i,j,k,zout,xin,yin,zin,b,adxinv);
            }
        });
    }
}

void MLCurlCurl::smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                         bool skip_fillboundary) const
{
    AMREX_ASSERT(rhs[0].nGrowVect().allGE(IntVect(1)));

    applyBC(amrlev, mglev, const_cast<MF&>(rhs), CurlCurlStateType::b);

    for (int color = 0; color < 4; ++color) {
        if (!skip_fillboundary) {
            applyBC(amrlev, mglev, sol, CurlCurlStateType::x);
        }
        skip_fillboundary = false;
        smooth4(amrlev, mglev, sol, rhs, color);
    }
}

void MLCurlCurl::smooth4 (int amrlev, int mglev, MF& sol, MF const& rhs,
                          int color) const
{
    auto const& ex = sol[0].arrays();
    auto const& ey = sol[1].arrays();
    auto const& ez = sol[2].arrays();
    auto const& rhsx = rhs[0].const_arrays();
    auto const& rhsy = rhs[1].const_arrays();
    auto const& rhsz = rhs[2].const_arrays();

#if (AMREX_SPACEDIM == 2)
    auto b = m_beta;
#endif

    auto adxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        adxinv[idim] *= std::sqrt(m_alpha);
    }

    auto* plusolver = m_lusolver[amrlev][mglev]->dataPtr();

    auto dinfo = getDirichletInfo(amrlev,mglev);
    auto sinfo = getSymmetryInfo(amrlev,mglev);

    MultiFab nmf(amrex::convert(rhs[0].boxArray(),IntVect(1)),
                 rhs[0].DistributionMap(), 1, 0, MFInfo().SetAlloc(false));
    ParallelFor(nmf, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
    {
        mlcurlcurl_gs4(i,j,k,ex[bno],ey[bno],ez[bno],rhsx[bno],rhsy[bno],rhsz[bno],
#if (AMREX_SPACEDIM == 2)
                       b,
#endif
                       adxinv,color,*plusolver,dinfo,sinfo);
    });
    Gpu::streamSynchronize();
}

void MLCurlCurl::solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                   const MF* /*crse_bcdata*/)
{
    BL_PROFILE("MLCurlCurl::solutionResidual()");
    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);
    compresid(amrlev, mglev, resid, b);
}

void MLCurlCurl::correctionResidual (int amrlev, int mglev, MF& resid, MF& x,
                                     const MF& b, BCMode bc_mode,
                                     const MF* crse_bcdata)
{
    AMREX_ALWAYS_ASSERT(bc_mode != BCMode::Inhomogeneous && crse_bcdata == nullptr);
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
    compresid(amrlev, mglev, resid, b);
}

void MLCurlCurl::compresid (int amrlev, int mglev, MF& resid, MF const& b) const
{
    auto dinfo = getDirichletInfo(amrlev,mglev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resid[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& xbx = mfi.tilebox(resid[0].ixType().toIntVect());
        Box const& ybx = mfi.tilebox(resid[1].ixType().toIntVect());
        Box const& zbx = mfi.tilebox(resid[2].ixType().toIntVect());
        auto const& resx = resid[0].array(mfi);
        auto const& resy = resid[1].array(mfi);
        auto const& resz = resid[2].array(mfi);
        auto const& bx = b[0].array(mfi);
        auto const& by = b[1].array(mfi);
        auto const& bz = b[2].array(mfi);
        amrex::ParallelFor(xbx, ybx, zbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_x_edge(i,j,k)) {
                resx(i,j,k) = Real(0.0);
            } else {
                resx(i,j,k) = bx(i,j,k) - resx(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_y_edge(i,j,k)) {
                resy(i,j,k) = Real(0.0);
            } else {
                resy(i,j,k) = by(i,j,k) - resy(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (dinfo.is_dirichlet_z_edge(i,j,k)) {
                resz(i,j,k) = Real(0.0);
            } else {
                resz(i,j,k) = bz(i,j,k) - resz(i,j,k);
            }
        });
    }
}

void MLCurlCurl::prepareForSolve ()
{
    for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            auto const& dxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
            Real dxx = dxinv[0]*dxinv[0];
            Real dyy = dxinv[1]*dxinv[1];
            Real dxy = dxinv[0]*dxinv[1];
#if (AMREX_SPACEDIM == 2)
            Array2D<Real,0,3,0,3,Order::C> A
                {m_alpha*dyy*Real(2.0) + m_beta,
                 Real(0.0),
                -m_alpha*dxy,
                 m_alpha*dxy,
                 //
                 Real(0.0),
                 m_alpha*dyy*Real(2.0) + m_beta,
                 m_alpha*dxy,
                -m_alpha*dxy,
                 //
                -m_alpha*dxy,
                 m_alpha*dxy,
                 m_alpha*dxx*Real(2.0) + m_beta,
                 Real(0.0),
                 //
                 m_alpha*dxy,
                -m_alpha*dxy,
                 Real(0.0),
                 m_alpha*dxx*Real(2.0) + m_beta};
#else
            Real dzz = dxinv[2]*dxinv[2];
            Real dxz = dxinv[0]*dxinv[2];
            Real dyz = dxinv[1]*dxinv[2];

            Array2D<Real,0,5,0,5,Order::C> A
                {m_alpha*(dyy+dzz)*Real(2.0) + m_beta,
                 Real(0.0),
                -m_alpha*dxy,
                 m_alpha*dxy,
                -m_alpha*dxz,
                 m_alpha*dxz,
                 //
                 Real(0.0),
                 m_alpha*(dyy+dzz)*Real(2.0) + m_beta,
                 m_alpha*dxy,
                -m_alpha*dxy,
                 m_alpha*dxz,
                -m_alpha*dxz,
                 //
                -m_alpha*dxy,
                 m_alpha*dxy,
                 m_alpha*(dxx+dzz)*Real(2.0) + m_beta,
                 Real(0.0),
                -m_alpha*dyz,
                 m_alpha*dyz,
                 //
                 m_alpha*dxy,
                -m_alpha*dxy,
                 Real(0.0),
                 m_alpha*(dxx+dzz)*Real(2.0) + m_beta,
                 m_alpha*dyz,
                -m_alpha*dyz,
                 //
                -m_alpha*dxz,
                 m_alpha*dxz,
                -m_alpha*dyz,
                 m_alpha*dyz,
                 m_alpha*(dxx+dyy)*Real(2.0) + m_beta,
                 Real(0.0),
                 //
                 m_alpha*dxz,
                -m_alpha*dxz,
                 m_alpha*dyz,
                -m_alpha*dyz,
                 Real(0.0),
                 m_alpha*(dxx+dyy)*Real(2.0) + m_beta};
#endif

            m_lusolver[amrlev][mglev]
                = std::make_unique<Gpu::DeviceScalar
                                   <LUSolver<AMREX_SPACEDIM*2,RT>>>(A);
        }
    }
}

Real MLCurlCurl::xdoty (int amrlev, int mglev, const MF& x, const MF& y,
                        bool local) const
{
    auto result = Real(0.0);
    for (int idim = 0; idim < 3; ++idim) {
        auto rtmp = MultiFab::Dot(getDotMask(amrlev,mglev,idim),
                                  x[idim], 0, y[idim], 0, 1, 0, true);
        result += rtmp;
    }
    if (!local) {
        ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    }
    return result;
}

Real MLCurlCurl::normInf (int /*amrlev*/, MF const& mf, bool local) const
{
    return amrex::norminf(mf, 0, m_ncomp, IntVect(0), local);
}

void MLCurlCurl::averageDownAndSync (Vector<MF>& sol) const
{
    BL_PROFILE("MLCurlCurl::averageDownAndSync()");
    AMREX_ALWAYS_ASSERT(sol.size() == 1);
    const int amrlev = 0;
    const int mglev = 0;
    for (int idim = 0; idim < 3; ++idim) {
        amrex::OverrideSync(sol[amrlev][idim], getDotMask(amrlev,mglev,idim),
                            this->m_geom[amrlev][mglev].periodicity());
    }
}

void MLCurlCurl::make (Vector<Vector<MF> >& mf, IntVect const& ng) const
{
    MLLinOpT<MF>::make(mf, ng);
}

Array<MultiFab,3>
MLCurlCurl::make (int amrlev, int mglev, IntVect const& ng) const
{
    MF r;
    for (int idim = 0; idim < 3; ++idim) {
        r[idim].define(amrex::convert(this->m_grids[amrlev][mglev], m_etype[idim]),
                       this->m_dmap[amrlev][mglev], m_ncomp, ng, MFInfo(),
                       *(this->m_factory)[amrlev][mglev]);
    }
    return r;
}

Array<MultiFab,3>
MLCurlCurl::makeAlias (MF const& mf) const
{
    MF r;
    for (int idim = 0; idim < 3; ++idim) {
        r[idim] = MultiFab(mf[idim], amrex::make_alias, 0, mf[idim].nComp());
    }
    return r;
}

Array<MultiFab,3>
MLCurlCurl::makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    BoxArray cba = this->m_grids[amrlev][mglev];
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[mglev];
    cba.coarsen(ratio);

    MF r;
    for (int idim = 0; idim < 3; ++idim) {
        r[idim].define(amrex::convert(cba, m_etype[idim]),
                       this->m_dmap[amrlev][mglev], m_ncomp, ng);
    }
    return r;
}

Array<MultiFab,3>
MLCurlCurl::makeCoarseAmr (int famrlev, IntVect const& ng) const
{
    BoxArray cba = this->m_grids[famrlev][0];
    IntVect ratio(this->AMRRefRatio(famrlev-1));
    cba.coarsen(ratio);

    MF r;
    for (int idim = 0; idim < 3; ++idim) {
        r[idim].define(amrex::convert(cba, m_etype[idim]),
                       this->m_dmap[famrlev][0], m_ncomp, ng);
    }
    return r;
}

void MLCurlCurl::applyBC (int amrlev, int mglev, MF& in, CurlCurlStateType type) const
{
    int nmfs = 3;
#if (AMREX_SPACEDIM == 2)
    if (CurlCurlStateType::b == type) {
        nmfs = 2; // no need to applyBC on Ez
    }
#endif
    Vector<MultiFab*> mfs(nmfs);
    for (int imf = 0; imf < nmfs; ++imf) {
        mfs[imf] = in.data() + imf;
    }
    FillBoundary(mfs, this->m_geom[amrlev][mglev].periodicity());
    for (auto* mf : mfs) {
        applyPhysBC(amrlev, mglev, *mf, type);
    }
}

#ifdef AMREX_USE_GPU
struct MLCurlCurlBCTag {
    Array4<Real> fab;
    Box bx;
    Orientation face;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};

struct MLCurlCurlEdgeBCTag {
    Array4<Real> fab;
    Box bx;
    Dim3 offset;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};
#endif

void MLCurlCurl::applyPhysBC (int amrlev, int mglev, MultiFab& mf, CurlCurlStateType type) const
{
    if (CurlCurlStateType::b == type) { return; }

    auto const idxtype = mf.ixType();
    Box const domain = amrex::convert(this->m_geom[amrlev][mglev].Domain(), idxtype);
    Box const gdomain = amrex::convert
        (this->m_geom[amrlev][mglev].growPeriodicDomain(1), idxtype);

    MFItInfo mfi_info{};

#ifdef AMREX_USE_GPU
    Vector<MLCurlCurlBCTag> tags;
    mfi_info.DisableDeviceSync();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,mfi_info); mfi.isValid(); ++mfi) {
        auto const& vbx = mfi.validbox();
        auto const& a = mf.array(mfi);
        for (OrientationIter oit; oit; ++oit) {
            Orientation const face = oit();
            int const idim = face.coordDir();
            bool is_symmetric = face.isLow()
                ? m_lobc[0][idim] == LinOpBCType::symmetry
                : m_hibc[0][idim] == LinOpBCType::symmetry;
            if (domain[face] == vbx[face] && is_symmetric &&
                ((type == CurlCurlStateType::x) ||
                 (type == CurlCurlStateType::r && idxtype.nodeCentered(idim)))) // transverse direction only
            {
                Box b = vbx;
                for (int jdim = 0; jdim < AMREX_SPACEDIM; ++jdim) {
                    if (jdim == idim) {
                        int shift = face.isLow() ? -1 : 1;
                        b.setRange(jdim, domain[face] + shift, 1);
                    } else {
                        if (b.smallEnd(jdim) > gdomain.smallEnd(jdim)) {
                            b.growLo(jdim);
                        }
                        if (b.bigEnd(jdim) < gdomain.bigEnd(jdim)) {
                            b.growHi(jdim);
                        }
                    }
                }
#ifdef AMREX_USE_GPU
                tags.emplace_back(MLCurlCurlBCTag{a,b,face});
#else
                amrex::LoopOnCpu(b, [&] (int i, int j, int k)
                {
                    mlcurlcurl_bc_symmetry(i, j, k, face, idxtype, a);
                });
#endif
            }
        }
    }

#ifdef AMREX_USE_GPU
    ParallelFor(tags,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, MLCurlCurlBCTag const& tag) noexcept
    {
        mlcurlcurl_bc_symmetry(i, j, k, tag.face, idxtype, tag.fab);
    });
#endif

    if (CurlCurlStateType::r == type) { // fix domain edges
        auto sinfo = getSymmetryInfo(amrlev,mglev);

#ifdef AMREX_USE_GPU
        Vector<MLCurlCurlEdgeBCTag> tags2;
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,mfi_info); mfi.isValid(); ++mfi) {
            auto const& vbx = mfi.validbox();
            auto const& a = mf.array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM-1; ++idim) {
                for (int jdim = idim+1; jdim < AMREX_SPACEDIM; ++jdim) {
                    if (idxtype.nodeCentered(idim) &&
                        idxtype.nodeCentered(jdim))
                    {
                        for (int iside = 0; iside < 2; ++iside) {
                            int ii = (iside == 0) ? vbx.smallEnd(idim) : vbx.bigEnd(idim);
                            for (int jside = 0; jside < 2; ++jside) {
                                int jj = (jside == 0) ? vbx.smallEnd(jdim) : vbx.bigEnd(jdim);
                                if (sinfo.is_symmetric(idim,iside,ii) &&
                                    sinfo.is_symmetric(jdim,jside,jj))
                                {
                                    IntVect oiv(0);
                                    oiv[idim] = (iside == 0) ? 2 : -2;
                                    oiv[jdim] = (jside == 0) ? 2 : -2;
                                    Dim3 offset = oiv.dim3();

                                    Box b = vbx;
                                    if (iside == 0) {
                                        b.setRange(idim,vbx.smallEnd(idim)-1);
                                    } else {
                                        b.setRange(idim,vbx.bigEnd(idim)+1);
                                    }
                                    if (jside == 0) {
                                        b.setRange(jdim,vbx.smallEnd(jdim)-1);
                                    } else {
                                        b.setRange(jdim,vbx.bigEnd(jdim)+1);
                                    }
#ifdef AMREX_USE_GPU
                                    tags2.emplace_back(MLCurlCurlEdgeBCTag{a,b,offset});
#else
                                    amrex::LoopOnCpu(b, [&] (int i, int j, int k)
                                    {
                                        a(i,j,k) = a(i+offset.x,j+offset.y,k+offset.z);
                                    });
#endif
                                }
                            }
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_GPU
        ParallelFor(tags2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, MLCurlCurlEdgeBCTag const& tag)
        {
            tag.fab(i,j,k) = tag.fab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z);
        });
#endif
    }
}

iMultiFab const& MLCurlCurl::getDotMask (int amrlev, int mglev, int idim) const
{
    if (m_dotmask[amrlev][mglev][idim] == nullptr) {
        MultiFab tmp(amrex::convert(this->m_grids[amrlev][mglev], m_etype[idim]),
                     this->m_dmap[amrlev][mglev], 1, 0, MFInfo().SetAlloc(false));
        m_dotmask[amrlev][mglev][idim] =
            tmp.OwnerMask(this->m_geom[amrlev][mglev].periodicity());
    }
    return *m_dotmask[amrlev][mglev][idim];
}

CurlCurlDirichletInfo MLCurlCurl::getDirichletInfo (int amrlev, int mglev) const
{

    auto helper = [&] (int idim, int face) -> int
    {
#if (AMREX_SPACEDIM == 2)
        if (idim == 2) {
            return std::numeric_limits<int>::lowest();
        }
#endif

        if (face == 0) {
            if (m_lobc[0][idim] == LinOpBCType::Dirichlet) {
                return m_geom[amrlev][mglev].Domain().smallEnd(idim);
            } else {
                return std::numeric_limits<int>::lowest();
            }
        } else {
            if (m_hibc[0][idim] == LinOpBCType::Dirichlet) {
                return m_geom[amrlev][mglev].Domain().bigEnd(idim) + 1;
            } else {
                return std::numeric_limits<int>::max();
            }
        }
    };

    return CurlCurlDirichletInfo{IntVect(AMREX_D_DECL(helper(0,0),
                                                      helper(1,0),
                                                      helper(2,0))),
                                 IntVect(AMREX_D_DECL(helper(0,1),
                                                      helper(1,1),
                                                      helper(2,1)))};
}

CurlCurlSymmetryInfo MLCurlCurl::getSymmetryInfo (int amrlev, int mglev) const
{

    auto helper = [&] (int idim, int face) -> int
    {
#if (AMREX_SPACEDIM == 2)
        if (idim == 2) {
            return std::numeric_limits<int>::lowest();
        }
#endif

        if (face == 0) {
            if (m_lobc[0][idim] == LinOpBCType::symmetry) {
                return m_geom[amrlev][mglev].Domain().smallEnd(idim);
            } else {
                return std::numeric_limits<int>::lowest();
            }
        } else {
            if (m_hibc[0][idim] == LinOpBCType::symmetry) {
                return m_geom[amrlev][mglev].Domain().bigEnd(idim) + 1;
            } else {
                return std::numeric_limits<int>::max();
            }
        }
    };

    return CurlCurlSymmetryInfo{IntVect(AMREX_D_DECL(helper(0,0),
                                                     helper(1,0),
                                                     helper(2,0))),
                                IntVect(AMREX_D_DECL(helper(0,1),
                                                     helper(1,1),
                                                     helper(2,1)))};
}

}
