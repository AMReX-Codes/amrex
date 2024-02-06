
#include <AMReX_MLCurlCurl.H>
#include <AMReX_MLCurlCurl_K.H>

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
}

void MLCurlCurl::setScalars (RT a_alpha, RT a_beta) noexcept
{
    m_alpha = a_alpha;
    m_beta = a_beta;
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
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        amrex::average_down_edges(fine[idim], crse[idim], ratio);
    }
#if (AMREX_SPACEDIM == 2)
    amrex::average_down_nodal(fine[2], crse[2], ratio);
#endif
}

void MLCurlCurl::interpolation (int amrlev, int fmglev, MF& fine,
                                const MF& crse) const
{
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];
    AMREX_ALWAYS_ASSERT(ratio == 2);

    for (int idim = 0; idim < 3; ++idim) {
        auto const& finema = fine[idim].arrays();
        auto const& crsema = crse[idim].const_arrays();
        ParallelFor(fine[idim], [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
        {
            mlcurlcurl_interpadd(idim,i,j,k,finema[bno],crsema[bno]);
        });
    }
    Gpu::streamSynchronize();
}

void
MLCurlCurl::apply (int amrlev, int mglev, MF& out, MF& in, BCMode /*bc_mode*/,
                   StateMode /*s_mode*/, const MLMGBndryT<MF>* /*bndry*/) const
{
    applyBC(amrlev, mglev, in);

    auto const& dxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
    auto const a = m_alpha;
    auto const b = m_beta;

    int const dirichlet_xlo = getDirichlet(amrlev, mglev, 0, 0);
    int const dirichlet_xhi = getDirichlet(amrlev, mglev, 0, 1);
    int const dirichlet_ylo = getDirichlet(amrlev, mglev, 1, 0);
    int const dirichlet_yhi = getDirichlet(amrlev, mglev, 1, 1);
    int const dirichlet_zlo = getDirichlet(amrlev, mglev, 2, 0);
    int const dirichlet_zhi = getDirichlet(amrlev, mglev, 2, 1);

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
            if (j == dirichlet_ylo || j == dirichlet_yhi ||
                k == dirichlet_zlo || k == dirichlet_zhi) {
                xout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_x(i,j,k,xout,xin,yin,zin,a,b,dxinv);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i == dirichlet_xlo || i == dirichlet_xhi ||
                k == dirichlet_zlo || k == dirichlet_zhi) {
                yout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_y(i,j,k,yout,xin,yin,zin,a,b,dxinv);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i == dirichlet_xlo || i == dirichlet_xhi ||
                j == dirichlet_ylo || j == dirichlet_yhi) {
                zout(i,j,k) = Real(0.0);
            } else {
                mlcurlcurl_adotx_z(i,j,k,zout,xin,yin,zin,a,b,dxinv);
            }
        });
    }
}

void MLCurlCurl::smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                         bool skip_fillboundary) const
{
    if (!skip_fillboundary) {
        applyBC(amrlev, mglev, sol);
    }

    smooth(amrlev, mglev, sol, rhs[0], 0); // Ex red
    applyBC(amrlev, mglev, sol[0]);

    smooth(amrlev, mglev, sol, rhs[1], 0); // Ey red
    applyBC(amrlev, mglev, sol[1]);

    smooth(amrlev, mglev, sol, rhs[2], 0); // Ez red
    applyBC(amrlev, mglev, sol[2]);

    smooth(amrlev, mglev, sol, rhs[0], 1); // Ex black
    applyBC(amrlev, mglev, sol[0]);

    smooth(amrlev, mglev, sol, rhs[1], 1); // Ey black
#if (AMREX_SPACEDIM == 3)
    applyBC(amrlev, mglev, sol[1]);
#endif

    smooth(amrlev, mglev, sol, rhs[2], 1); // Ez black

    for (int idim = 0; idim < 3; ++idim) {
        amrex::OverrideSync(sol[idim], getDotMask(amrlev,mglev,idim),
                            this->m_geom[amrlev][mglev].periodicity());
    }
}

void MLCurlCurl::smooth (int amrlev, int mglev, MF& sol, MultiFab const& rhs,
                         int redblack) const
{
    auto const& dxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
    auto const a = m_alpha;
    auto const b = m_beta;

    int const dirichlet_xlo = getDirichlet(amrlev, mglev, 0, 0);
    int const dirichlet_xhi = getDirichlet(amrlev, mglev, 0, 1);
    int const dirichlet_ylo = getDirichlet(amrlev, mglev, 1, 0);
    int const dirichlet_yhi = getDirichlet(amrlev, mglev, 1, 1);
    int const dirichlet_zlo = getDirichlet(amrlev, mglev, 2, 0);
    int const dirichlet_zhi = getDirichlet(amrlev, mglev, 2, 1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto const& rh = rhs.const_array(mfi);
        if (rhs.ixType() == sol[0].ixType()) {
            auto const& ex = sol[0].array(mfi);
            auto const& ey = sol[1].const_array(mfi);
            auto const& ez = sol[2].const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (j != dirichlet_ylo && j != dirichlet_yhi &&
                    k != dirichlet_zlo && k != dirichlet_zhi) {
                    mlcurlcurl_gsrb_x(i,j,k,ex,ey,ez,rh,a,b,dxinv,redblack);
                }
            });
        } else if (rhs.ixType() == sol[1].ixType()) {
            auto const& ex = sol[0].const_array(mfi);
            auto const& ey = sol[1].array(mfi);
            auto const& ez = sol[2].const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (i != dirichlet_xlo && i != dirichlet_xhi &&
                    k != dirichlet_zlo && k != dirichlet_zhi) {
                    mlcurlcurl_gsrb_y(i,j,k,ex,ey,ez,rh,a,b,dxinv,redblack);
                }
            });
        } else {
            auto const& ex = sol[0].const_array(mfi);
            auto const& ey = sol[1].const_array(mfi);
            auto const& ez = sol[2].array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (i != dirichlet_xlo && i != dirichlet_xhi &&
                    j != dirichlet_ylo && j != dirichlet_yhi) {
                    mlcurlcurl_gsrb_z(i,j,k,ex,ey,ez,rh,a,b,dxinv,redblack);
                }
            });
        }
    }
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
    int const dirichlet_xlo = getDirichlet(amrlev, mglev, 0, 0);
    int const dirichlet_xhi = getDirichlet(amrlev, mglev, 0, 1);
    int const dirichlet_ylo = getDirichlet(amrlev, mglev, 1, 0);
    int const dirichlet_yhi = getDirichlet(amrlev, mglev, 1, 1);
    int const dirichlet_zlo = getDirichlet(amrlev, mglev, 2, 0);
    int const dirichlet_zhi = getDirichlet(amrlev, mglev, 2, 1);

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
            if (j == dirichlet_ylo || j == dirichlet_yhi ||
                k == dirichlet_zlo || k == dirichlet_zhi) {
                resx(i,j,k) = Real(0.0);
            } else {
                resx(i,j,k) = bx(i,j,k) - resx(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i == dirichlet_xlo || i == dirichlet_xhi ||
                k == dirichlet_zlo || k == dirichlet_zhi) {
                resy(i,j,k) = Real(0.0);
            } else {
                resy(i,j,k) = by(i,j,k) - resy(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i == dirichlet_xlo || i == dirichlet_xhi ||
                j == dirichlet_ylo || j == dirichlet_yhi) {
                resz(i,j,k) = Real(0.0);
            } else {
                resz(i,j,k) = bz(i,j,k) - resz(i,j,k);
            }
        });
    }
}

void MLCurlCurl::prepareForSolve ()
{
}

Real MLCurlCurl::xdoty (int amrlev, int mglev, const MF& x, const MF& y,
                        bool local) const
{
    auto result = Real(0.0);
    for (int idim = 0; idim < 3; ++idim) {
        auto rtmp = MultiFab::Dot(getDotMask(amrlev,mglev,idim),
                                  x[idim], 0, y[idim], 0, 1, 0, false);
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

void MLCurlCurl::applyBC (int amrlev, int mglev, MF& in) const
{
    Vector<MultiFab*> mfs{in.data(),&(in[1]),&(in[2])};
    FillBoundary(mfs, this->m_geom[amrlev][mglev].periodicity());
    for (auto& mf : in) {
        applyPhysBC(amrlev, mglev, mf);
    }
}

void MLCurlCurl::applyBC (int amrlev, int mglev, MultiFab& mf) const
{
    mf.FillBoundary(this->m_geom[amrlev][mglev].periodicity());
    applyPhysBC(amrlev, mglev, mf);
}

#ifdef AMREX_USE_GPU
struct MLCurlCurlBCTag {
    Array4<Real> fab;
    Box bx;
    Orientation face;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};
#endif

void MLCurlCurl::applyPhysBC (int amrlev, int mglev, MultiFab& mf) const
{
    auto const idxtype = mf.ixType();
    Box const domain = amrex::convert(this->m_geom[amrlev][mglev].Domain(), idxtype);

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
            if (domain[face] == vbx[face] && is_symmetric) {
                Box b = vbx;
                int shift = face.isLow() ? -1 : 1;
                b.setRange(idim, domain[face] + shift, 1);
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

int MLCurlCurl::getDirichlet (int amrlev, int mglev, int idim, int face) const
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
}

}
