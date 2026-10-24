#include <AMReX_HypreMLNodeLap.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLinOp_K.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Scan.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <array>
#include <numeric>

namespace amrex {

namespace {

constexpr int nnodes_cell = AMREX_D_TERM(2,*2,*2); // corners of a cell

MFInfo host_info ()
{
#ifdef AMREX_USE_GPU
    return MFInfo().SetArena(The_Pinned_Arena());
#else
    return MFInfo();
#endif
}

// Copy a (device) iMultiFab into a host-accessible one on the same layout
std::unique_ptr<iMultiFab> host_copy (iMultiFab const& src)
{
    auto dst = std::make_unique<iMultiFab>(src.boxArray(), src.DistributionMap(),
                                           src.nComp(), src.nGrowVect(), host_info());
    iMultiFab::Copy(*dst, src, 0, 0, src.nComp(), src.nGrowVect());
    Gpu::streamSynchronize();
    return dst;
}

std::unique_ptr<MultiFab> host_copy (MultiFab const& src, int ngrow)
{
    auto dst = std::make_unique<MultiFab>(src.boxArray(), src.DistributionMap(),
                                          src.nComp(), ngrow, host_info());
    MultiFab::Copy(*dst, src, 0, 0, src.nComp(), ngrow);
    Gpu::streamSynchronize();
    return dst;
}

// Multilinear interpolation of fine node n from the coarse nodes of ratio rr:
// fills the coarse node indices and weights, returns their number.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int interp_stencil (IntVect const& n, int rr, IntVect* cnodes, Real* wts)
{
    int m = 1;
    cnodes[0] = IntVect(0);
    wts[0] = Real(1.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int lo = amrex::coarsen(n[d], rr);
        const int rem = n[d] - lo*rr;
        if (rem == 0) {
            for (int k = 0; k < m; ++k) { cnodes[k][d] = lo; }
        } else {
            const Real f = Real(rem)/Real(rr);
            for (int k = 0; k < m; ++k) {
                cnodes[k+m] = cnodes[k];
                wts[k+m] = wts[k]*f;
                cnodes[k][d] = lo;
                cnodes[k+m][d] = lo+1;
                wts[k] *= (Real(1.0)-f);
            }
            m *= 2;
        }
    }
    return m;
}

// Q1 element stiffness of a cell with per-direction sigma:
// K = sum_d sigma_d * (x)_{d' != d} M1(h_d') (x) S1(h_d),
// M1 = h [1/3 1/6; 1/6 1/3], S1 = (1/h) [1 -1; -1 1].
// Corner c has bit d set if it is at the high side in direction d.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void element_matrix (GpuArray<Real,AMREX_SPACEDIM> const& h,
                     GpuArray<Real,AMREX_SPACEDIM> const& sig,
                     Real (&K)[nnodes_cell][nnodes_cell])
{
    for (int a = 0; a < nnodes_cell; ++a) {
        for (int b = 0; b < nnodes_cell; ++b) {
            Real kab = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real t = sig[d];
                for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                    const bool same = ((a >> e) & 1) == ((b >> e) & 1);
                    if (e == d) {
                        t *= (same ? Real(1.0) : Real(-1.0)) / h[e];
                    } else {
                        t *= h[e] * (same ? Real(1.0/3.0) : Real(1.0/6.0));
                    }
                }
                kab += t;
            }
            K[a][b] = kab;
        }
    }
}

#ifdef AMREX_USE_EB
// Moment of x^px y^py z^pz (0 <= p <= 2, at most two nonzero powers) over the
// fluid part of a cell, in cell coordinates in [-1/2,1/2] (MLNodeLaplacian's
// integrals; the zeroth moment is the volume fraction).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real cell_moment (Array4<Real const> const& intg, Real vf, IntVect const& c,
                  IntVect const& p)
{
#if (AMREX_SPACEDIM == 2)
    switch (p[0] + 3*p[1]) {
    case 0: return vf;
    case 1: return intg(c,0);   // x
    case 3: return intg(c,1);   // y
    case 2: return intg(c,2);   // x2
    case 6: return intg(c,3);   // y2
    case 4: return intg(c,4);   // xy
    default: return Real(0.0);
    }
#else
    switch (p[0] + 3*p[1] + 9*p[2]) {
    case  0: return vf;
    case  1: return intg(c, 0); // x
    case  3: return intg(c, 1); // y
    case  9: return intg(c, 2); // z
    case  2: return intg(c, 3); // x2
    case  6: return intg(c, 4); // y2
    case 18: return intg(c, 5); // z2
    case  4: return intg(c, 6); // x y
    case 10: return intg(c, 7); // x z
    case 12: return intg(c, 8); // y z
    case  5: return intg(c, 9); // x2 y
    case 11: return intg(c,10); // x2 z
    case  7: return intg(c,11); // x y2
    case 15: return intg(c,12); // y2 z
    case 19: return intg(c,13); // x z2
    case 21: return intg(c,14); // y z2
    case  8: return intg(c,15); // x2 y2
    case 20: return intg(c,16); // x2 z2
    case 24: return intg(c,17); // y2 z2
    default: return Real(0.0);
    }
#endif
}

// Q1 element stiffness of a cut cell with scalar sigma, integrated over the
// fluid part of the cell: with cell coordinates x_e in [-1/2,1/2] and
// s_ae = -1/+1 for corner a at the low/high side in direction e,
// K_ab = sigma sum_d s_ad s_bd |cell|/h_d^2 int prod_{e != d} (1/2 + s_ae x_e)(1/2 + s_be x_e).
// This is the per-cell contribution to MLNodeLaplacian's EB stencil.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void element_matrix_eb (GpuArray<Real,AMREX_SPACEDIM> const& h, Real sig,
                        Array4<Real const> const& intg, Real vf, IntVect const& c,
                        Real (&K)[nnodes_cell][nnodes_cell])
{
    const Real vol = AMREX_D_TERM(h[0],*h[1],*h[2]);
    for (int a = 0; a < nnodes_cell; ++a) {
        for (int b = 0; b < nnodes_cell; ++b) {
            Real kab = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // 1D factors c0 + c1 x + c2 x^2 in the other directions
                int dirs[2] = {-1, -1};
                Real cf[2][3] = {{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
                int nd = 0;
                for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                    if (e == d) { continue; }
                    const int sa = ((a >> e) & 1) ? 1 : -1;
                    const int sb = ((b >> e) & 1) ? 1 : -1;
                    dirs[nd] = e;
                    cf[nd][0] = Real(0.25);
                    cf[nd][1] = (sa == sb) ? Real(sa) : Real(0.0);
                    cf[nd][2] = (sa == sb) ? Real(1.0) : Real(-1.0);
                    ++nd;
                }
                Real integral = 0.0;
                for (int p = 0; p < 3; ++p) {
                    for (int q = 0; q < ((nd > 1) ? 3 : 1); ++q) {
                        const Real coef = cf[0][p] * cf[1][q];
                        if (coef == Real(0.0)) { continue; }
                        IntVect pw(0);
                        pw[dirs[0]] = p;
                        if (nd > 1) { pw[dirs[1]] = q; }
                        integral += coef * cell_moment(intg, vf, c, pw);
                    }
                }
                const int sad = ((a >> d) & 1) ? 1 : -1;
                const int sbd = ((b >> d) & 1) ? 1 : -1;
                kab += Real(sad*sbd) * vol / (h[d]*h[d]) * integral;
            }
            K[a][b] = sig * kab;
        }
    }
}
#endif

#ifdef AMREX_USE_EB
// Host copy of the fluid volume fraction of the nodes' control volumes
// (cells outside dom do not count)
std::unique_ptr<MultiFab> nodal_fluid_fraction (BoxArray const& nd_ba,
                                                DistributionMapping const& dm,
                                                EBFArrayBoxFactory const& fact,
                                                Box const& dom)
{
    MultiFab ndv(nd_ba, dm, 1, 0);
    auto const& ndv_ma = ndv.arrays();
    auto const& vf_ma = fact.getVolFrac().const_arrays();
    constexpr Real fac = Real(1.0)/Real(nnodes_cell);
    ParallelFor(ndv, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
    {
        Real v = 0.0;
        for (int a = 0; a < nnodes_cell; ++a) {
            IntVect c(AMREX_D_DECL(i,j,k));
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { c[d] -= 1 - ((a >> d) & 1); }
            if (dom.contains(c)) { v += fac * vf_ma[bno](c); }
        }
        ndv_ma[bno](i,j,k) = v;
    });
    Gpu::streamSynchronize();
    return host_copy(ndv, 0);
}

// EB factory of level lev, or nullptr without EB or if all regular
EBFArrayBoxFactory const* eb_factory (FabFactory<FArrayBox> const* factory)
{
    auto const* f = dynamic_cast<EBFArrayBoxFactory const*>(factory);
    return (f && !f->isAllRegular()) ? f : nullptr;
}
#endif

} // namespace

HypreMLNodeLap::HypreMLNodeLap (MLNodeLaplacian const& a_linop, std::string a_prefix,
                                int a_verbose)
    : m_linop(&a_linop), m_prefix(std::move(a_prefix)), m_verbose(a_verbose)
{
    BL_PROFILE("HypreMLNodeLap::HypreMLNodeLap()");

    static_assert(std::is_same_v<Real, HYPRE_Real>, "amrex::Real != HYPRE_Real");

    // RAP (EB builds): the same element sum; it halves the Neumann
    // right-hand side instead of doubling those rows, so the right-hand side
    // map is as for Sigma.
    const bool rap = (m_linop->m_coarsening_strategy == MLNodeLaplacian::CoarseningStrategy::RAP);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!rap || !m_linop->m_use_mapped,
        "HypreMLNodeLap: mesh mapping requires the Sigma coarsening strategy");

    m_comm = m_linop->Communicator();
    m_nprocs = ParallelDescriptor::NProcs(m_comm);
    m_myproc = ParallelDescriptor::MyProc(m_comm);

    m_nlevels = m_linop->NAMRLevels();
    for (int lev = 0; lev < m_nlevels; ++lev) {
        m_geom.push_back(m_linop->m_geom[lev][0]);
        m_nd_grids.push_back(amrex::convert(m_linop->m_grids[lev][0], IntVect::TheNodeVector()));
        m_dmap.push_back(m_linop->m_dmap[lev][0]);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_geom[lev].IsRZ(),
                                         "HypreMLNodeLap: RZ is not supported");
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!m_geom[0].isPeriodic(d) && (m_linop->m_lobc[0][d] == LinOpBCType::Dirichlet ||
                                         m_linop->m_hibc[0][d] == LinOpBCType::Dirichlet)) {
            m_has_dirichlet = true;
        }
    }
    for (int lev = 0; lev < m_nlevels-1; ++lev) {
        const int rr = m_linop->AMRRefRatio(lev);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rr == 2 || rr == 4,
                                         "HypreMLNodeLap: refinement ratio must be 2 or 4");
        m_ref_ratio.push_back(rr);
    }

    {
        ParmParse pp(m_prefix);
        pp.query("eb_scaling", m_eb_scaling);
        std::string backend = "hypre";
        pp.query("backend", backend);
        if (backend == "algmg") {
            m_backend = Backend::algmg;
        } else if (backend != "hypre") {
            amrex::Abort("HypreMLNodeLap: unknown " + m_prefix + ".backend = " + backend +
                         " (hypre or algmg)");
        }
        pp.query("dump_matrix", m_dump_max);
        m_dump_prefix = m_prefix;
        pp.query("dump_prefix", m_dump_prefix);
        pp.query("dump_gather_max", m_dump_gather_max);
    }

    numberUnknowns();

    // Row ids are HYPRE_BigInt on the device, and HypreIJIface takes the row
    // bounds as HYPRE_Int; both are 32-bit by default.
    using IJInt = HypreIJIface::HypreIntType;
    Long nrows_max = Long(std::numeric_limits<BigInt>::max());
    if (m_backend == Backend::hypre) {
        nrows_max = std::min(nrows_max, Long(std::numeric_limits<IJInt>::max()));
    }
    // m_nglobal has already been narrowed to BigInt: check the Long total
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_rank_offset.back() <= nrows_max,
        "HypreMLNodeLap: too many unknowns for this hypre build's integer type; "
        "build hypre with --enable-bigint");

    if (m_backend == Backend::hypre) {
        m_ij = std::make_unique<HypreIJIface>(m_comm, IJInt(m_ilower), IJInt(m_iupper), m_verbose);
        m_ij->parse_inputs(m_prefix);
        m_A = m_ij->A();
        m_b = m_ij->b();
        m_x = m_ij->x();
    }

    buildMatrix();
}

HypreMLNodeLap::~HypreMLNodeLap () = default;

FabFactory<FArrayBox> const* HypreMLNodeLap::factory (int lev) const
{
    return m_linop->Factory(lev);
}

bool HypreMLNodeLap::onDirichletFace (int lev, IntVect const& iv) const
{
    const Box& domain = m_geom[lev].Domain();
    const auto& lobc = m_linop->m_lobc[0];
    const auto& hibc = m_linop->m_hibc[0];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom[lev].isPeriodic(d)) { continue; }
        if (iv[d] == domain.smallEnd(d)   && lobc[d] == LinOpBCType::Dirichlet) { return true; }
        if (iv[d] == domain.bigEnd(d)+1   && hibc[d] == LinOpBCType::Dirichlet) { return true; }
    }
    return false;
}

int HypreMLNodeLap::neumannScale (int lev, IntVect const& iv) const
{
    const Box& domain = m_geom[lev].Domain();
    const auto& lobc = m_linop->m_lobc[0];
    const auto& hibc = m_linop->m_hibc[0];
    int ns = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom[lev].isPeriodic(d)) { continue; }
        if (iv[d] == domain.smallEnd(d) &&
            (lobc[d] == LinOpBCType::Neumann || lobc[d] == LinOpBCType::inflow)) { ns *= 2; }
        if (iv[d] == domain.bigEnd(d)+1 &&
            (hibc[d] == LinOpBCType::Neumann || hibc[d] == LinOpBCType::inflow)) { ns *= 2; }
    }
    return ns;
}

void HypreMLNodeLap::numberUnknowns ()
{
    BL_PROFILE("HypreMLNodeLap::numberUnknowns()");

    const int finest = m_nlevels-1;

    m_gid.resize(m_nlevels);
    m_owner.resize(m_nlevels);
    m_dmsk.resize(m_nlevels);
    m_ndfine.resize(m_nlevels);
    m_cgid_on_f.resize(m_nlevels);
    m_cndfine_on_f.resize(m_nlevels);

    // Local numbering per level: count the unknowns this rank owns
    Long nlocal = 0;
    m_ndvol.resize(m_nlevels);
    m_scale.resize(m_nlevels);
    m_cscale_on_f.resize(m_nlevels);
    for (int lev = 0; lev < m_nlevels; ++lev) {
        const auto& period = m_geom[lev].periodicity();
        m_gid[lev] = std::make_unique<IdFabArray>(m_nd_grids[lev], m_dmap[lev], 1, 0, host_info());
        {
            auto owner = amrex::OwnerMask(*m_gid[lev], period);
            m_owner[lev] = host_copy(*owner);
        }
        m_dmsk[lev] = host_copy(*m_linop->m_dirichlet_mask[lev][0]);
        if (lev < finest) {
            m_ndfine[lev] = host_copy(*m_linop->m_nd_fine_mask[lev]);
        }
#ifdef AMREX_USE_EB
        // Fluid volume fraction of each node's control volume (the weights
        // of MLNodeLaplacian's solvability offset); nodes without fluid
        // have no equation and are not unknowns.
        if (auto const* fact = eb_factory(factory(lev))) {
            Box dom = m_geom[lev].Domain();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                auto const blo = m_linop->m_lobc[0][d];
                auto const bhi = m_linop->m_hibc[0][d];
                if (blo != LinOpBCType::Neumann && blo != LinOpBCType::inflow) { dom.growLo(d, 10); }
                if (bhi != LinOpBCType::Neumann && bhi != LinOpBCType::inflow) { dom.growHi(d, 10); }
            }
            m_ndvol[lev] = nodal_fluid_fraction(m_nd_grids[lev], m_dmap[lev], *fact, dom);
            if (m_eb_scaling) {
                // Unknowns scaled by 1/sqrt(fluid fraction), roughly the
                // inverse square root of the diagonal, so that the rows of
                // nodes next to small cut cells count in the residual norm
                m_scale[lev] = std::make_unique<MultiFab>(m_nd_grids[lev], m_dmap[lev], 1, 0,
                                                          host_info());
                for (MFIter mfi(*m_scale[lev]); mfi.isValid(); ++mfi) {
                    auto const& sc = m_scale[lev]->array(mfi);
                    auto const& nv = m_ndvol[lev]->const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
                    {
                        sc(i,j,k) = (nv(i,j,k) > Real(0.0)) ? Real(1.0)/std::sqrt(nv(i,j,k))
                                                            : Real(1.0);
                    });
                }
            }
        }
#endif

        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& gid = m_gid[lev]->array(mfi);
            auto const& own = m_owner[lev]->const_array(mfi);
            auto const& dm = m_dmsk[lev]->const_array(mfi);
            Array4<int const> ndf;
            if (lev < finest) { ndf = m_ndfine[lev]->const_array(mfi); }
            Array4<Real const> ndv;
            if (m_ndvol[lev]) { ndv = m_ndvol[lev]->const_array(mfi); }
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
            {
                const bool unknown = own(i,j,k) && !dm(i,j,k) &&
                    (lev == finest || ndf(i,j,k) != nodelap_detail::fine_node) &&
                    (!ndv || ndv(i,j,k) > Real(0.0));
                gid(i,j,k) = unknown ? nlocal++ : Long(-1);
            });
        }
    }

    // Rank offsets
    Vector<Long> counts(m_nprocs, nlocal);
#ifdef AMREX_USE_MPI
    MPI_Allgather(&nlocal, 1, ParallelDescriptor::Mpi_typemap<Long>::type(),
                  counts.data(), 1, ParallelDescriptor::Mpi_typemap<Long>::type(), m_comm);
#endif
    m_rank_offset.assign(m_nprocs+1, 0);
    for (int p = 0; p < m_nprocs; ++p) { m_rank_offset[p+1] = m_rank_offset[p] + counts[p]; }
    const Long os = m_rank_offset[m_myproc];
    m_ilower = static_cast<BigInt>(os);
    m_iupper = static_cast<BigInt>(os + nlocal - 1);
    m_nglobal = static_cast<BigInt>(m_rank_offset[m_nprocs]);

    // Global ids, shared with the non-owned copies of each node
    for (int lev = 0; lev < m_nlevels; ++lev) {
        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                if (gid(i,j,k) >= 0) { gid(i,j,k) += os; }
            });
        }
        // Non-owners: take the owner's id (all copies of a node agree afterwards)
        auto tmp = std::make_unique<IdFabArray>(m_nd_grids[lev], m_dmap[lev], 1, 0, host_info());
        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->array(mfi);
            auto const& t = tmp->array(mfi);
            auto const& own = m_owner[lev]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                // store gid+1 at owners so that "not an unknown" becomes 0
                t(i,j,k) = own(i,j,k) ? gid(i,j,k)+1 : 0;
            });
        }
        m_gid[lev]->setVal(0);
        m_gid[lev]->ParallelAdd(*tmp, 0, 0, 1, m_geom[lev].periodicity());
        Gpu::streamSynchronize();
        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                gid(i,j,k) -= 1;
            });
        }
    }

    // Coarse ids and nodal fine mask on the coarsened fine layout
    for (int lev = 1; lev < m_nlevels; ++lev) {
        const int rr = m_ref_ratio[lev-1];
        const BoxArray cba = amrex::coarsen(m_nd_grids[lev], rr);
        m_cgid_on_f[lev] = std::make_unique<IdFabArray>(cba, m_dmap[lev], 1, 0, host_info());
        m_cgid_on_f[lev]->setVal(-1);
        m_cgid_on_f[lev]->ParallelCopy(*m_gid[lev-1], 0, 0, 1, IntVect(0), IntVect(0),
                                       m_geom[lev-1].periodicity());
        m_cndfine_on_f[lev] = std::make_unique<iMultiFab>(cba, m_dmap[lev], 1, 0, host_info());
        m_cndfine_on_f[lev]->setVal(0);
        m_cndfine_on_f[lev]->ParallelCopy(*m_ndfine[lev-1], 0, 0, 1, IntVect(0), IntVect(0),
                                          m_geom[lev-1].periodicity());
        if (m_scale[lev-1] || m_scale[lev]) {
            m_cscale_on_f[lev] = std::make_unique<MultiFab>(cba, m_dmap[lev], 1, 0, host_info());
            m_cscale_on_f[lev]->setVal(1.0);
            if (m_scale[lev-1]) {
                m_cscale_on_f[lev]->ParallelCopy(*m_scale[lev-1], 0, 0, 1, IntVect(0), IntVect(0),
                                                 m_geom[lev-1].periodicity());
            }
        }
        Gpu::streamSynchronize();
    }
    m_ndfine.clear();   // only needed for the numbering

    if (m_verbose > 0) {
        amrex::Print() << "HypreMLNodeLap: " << m_nglobal << " composite unknowns on "
                       << m_nlevels << " level(s)\n";
    }
}

void HypreMLNodeLap::buildMatrix ()
{
    BL_PROFILE("HypreMLNodeLap::buildMatrix()");

    // Every cell emits the (row, column, value) triplets of its element
    // matrix mapped through the hanging-node constraints, on the device.
    // They go to hypre (which sums duplicates and routes off-rank rows), to
    // the AlgMG matrix (built on the device) and, for dumps, to an owned-row
    // CSR on the host, one AMR level at a time.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_dbc.empty(), "HypreMLNodeLap::buildMatrix called twice");
    const int finest = m_nlevels-1;
    m_pin_gid = (m_linop->isBottomSingular() && m_nglobal > 0) ? 0 : -1;
    Long const pin = m_pin_gid;
    const bool need_csr = composite::dump_pending(m_dump_prefix, m_dump_max);
    const bool algmg_dev = (m_backend == Backend::algmg) && !need_csr;

    struct Triplets {
        Gpu::DeviceVector<BigInt> rows, cols;
        Gpu::DeviceVector<Real> vals;
    };
    Vector<Triplets> chunks;
    Vector<composite::Triplet> host_trips;
    Gpu::DeviceVector<Long> alg_rows, alg_cols;
    Gpu::DeviceVector<Real> alg_vals;
    auto flush = [&] ()
    {
        Long ntot = 0;
        for (auto const& t : chunks) { ntot += Long(t.rows.size()); }
        Gpu::DeviceVector<BigInt> rows(ntot), cols(ntot);
        Gpu::DeviceVector<Real> vals(ntot);
        Long off = 0;
        for (auto const& t : chunks) {
            Gpu::copyAsync(Gpu::deviceToDevice, t.rows.begin(), t.rows.end(), rows.begin()+off);
            Gpu::copyAsync(Gpu::deviceToDevice, t.cols.begin(), t.cols.end(), cols.begin()+off);
            Gpu::copyAsync(Gpu::deviceToDevice, t.vals.begin(), t.vals.end(), vals.begin()+off);
            off += Long(t.rows.size());
        }
        Gpu::streamSynchronize();
        chunks.clear();
        if (m_backend == Backend::hypre && ntot > 0) {
            Gpu::DeviceVector<HYPRE_Int> ones(ntot, 1);
            // HYPRE_IJMatrixAddToValues takes the count as HYPRE_Int
            constexpr Long max_chunk = std::numeric_limits<HYPRE_Int>::max();
            for (Long n0 = 0; n0 < ntot; n0 += max_chunk) {
                auto const n = HYPRE_Int(std::min(max_chunk, ntot-n0));
                HYPRE_IJMatrixAddToValues(m_A, n, ones.data()+n0, rows.data()+n0,
                                          cols.data()+n0, vals.data()+n0);
                Gpu::hypreSynchronize();
            }
        }
        if (algmg_dev && ntot > 0) {
            const auto n0 = Long(alg_rows.size());
            alg_rows.resize(n0+ntot);
            alg_cols.resize(n0+ntot);
            alg_vals.resize(n0+ntot);
            auto const* r = rows.data();
            auto const* c = cols.data();
            auto const* v = vals.data();
            auto* ar = alg_rows.data() + n0;
            auto* ac = alg_cols.data() + n0;
            auto* av = alg_vals.data() + n0;
            amrex::ParallelFor(ntot, [=] AMREX_GPU_DEVICE (Long n) noexcept
            {
                ar[n] = r[n];
                ac[n] = c[n];
                av[n] = v[n];
            });
            Gpu::streamSynchronize();
        }
        if (need_csr && ntot > 0) {
            Vector<BigInt> hr(ntot), hc(ntot);
            Vector<Real> hv(ntot);
            Gpu::copyAsync(Gpu::deviceToHost, rows.begin(), rows.end(), hr.begin());
            Gpu::copyAsync(Gpu::deviceToHost, cols.begin(), cols.end(), hc.begin());
            Gpu::copyAsync(Gpu::deviceToHost, vals.begin(), vals.end(), hv.begin());
            Gpu::streamSynchronize();
            const Long n0 = host_trips.size();
            host_trips.resize(n0 + ntot);
            for (Long n = 0; n < ntot; ++n) {
                auto& t = host_trips[n0+n];
                t.row = hr[n];
                t.col = hc[n];
                t.val = hv[n];
            }
        }
    };

    for (int lev = 0; lev < m_nlevels; ++lev) {
        // Device copies of the ids
        IdFabArray gid_d(m_nd_grids[lev], m_dmap[lev], 1, 0);
        amrex::Copy(gid_d, *m_gid[lev], 0, 0, 1, 0);
        IdFabArray cgid_d;
        if (lev > 0) {
            cgid_d.define(m_cgid_on_f[lev]->boxArray(), m_dmap[lev], 1, 0);
            amrex::Copy(cgid_d, *m_cgid_on_f[lev], 0, 0, 1, 0);
        }
        // Device copies of the unknown scaling (S' = s S s)
        MultiFab sc_d, csc_d;
        const bool has_sc = static_cast<bool>(m_scale[lev]);
        const bool has_csc = (lev > 0) && static_cast<bool>(m_cscale_on_f[lev]);
        if (has_sc) {
            sc_d.define(m_nd_grids[lev], m_dmap[lev], 1, 0);
            MultiFab::Copy(sc_d, *m_scale[lev], 0, 0, 1, 0);
        }
        if (has_csc) {
            csc_d.define(m_cscale_on_f[lev]->boxArray(), m_dmap[lev], 1, 0);
            MultiFab::Copy(csc_d, *m_cscale_on_f[lev], 0, 0, 1, 0);
        }

        const auto dx = m_geom[lev].CellSizeArray();
        const Real const_sigma = m_linop->m_const_sigma;
        const bool per_dir = m_linop->m_use_mapped ||
            (const_sigma == Real(0.0) && m_linop->m_sigma[lev][0][1] != nullptr);
        const int rr = (lev > 0) ? m_ref_ratio[lev-1] : 1;
        const bool has_fine = (lev < finest);

        // Dirichlet domain faces, for the device
        Box const domain = m_geom[lev].Domain();
        GpuArray<int,AMREX_SPACEDIM> dlo{}, dhi{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bool const per = m_geom[lev].isPeriodic(d);
            dlo[d] = (!per && m_linop->m_lobc[0][d] == LinOpBCType::Dirichlet) ? 1 : 0;
            dhi[d] = (!per && m_linop->m_hibc[0][d] == LinOpBCType::Dirichlet) ? 1 : 0;
        }

        for (MFIter mfi(gid_d); mfi.isValid(); ++mfi) {
            const Box ccbx = amrex::enclosedCells(mfi.validbox());
            auto const& gid = gid_d.const_array(mfi);
            auto const& dm = m_linop->m_dirichlet_mask[lev][0]->const_array(mfi);
            Array4<int const> ccf;
            if (has_fine) { ccf = m_linop->m_cc_fine_mask[lev]->const_array(mfi); }
            Array4<Long const> cgid;
            if (lev > 0) { cgid = cgid_d.const_array(mfi); }
            Array4<Real const> sc, csc;
            if (has_sc) { sc = sc_d.const_array(mfi); }
            if (has_csc) { csc = csc_d.const_array(mfi); }
            GpuArray<Array4<Real const>,AMREX_SPACEDIM> sig{};
            if (const_sigma == Real(0.0)) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (d == 0 || per_dir) { sig[d] = m_linop->m_sigma[lev][0][d]->const_array(mfi); }
                }
            }
#ifdef AMREX_USE_EB
            // Cut cells: element matrices from the fluid-part integrals
            Array4<EBCellFlag const> flag;
            Array4<Real const> vfrac, intg;
            bool eb = false;
            if (auto const* fact = eb_factory(factory(lev))) {
                if (fact->getMultiEBCellFlagFab()[mfi].getType(ccbx) != FabType::regular) {
                    eb = true;
                    flag = fact->getMultiEBCellFlagFab().const_array(mfi);
                    vfrac = fact->getVolFrac().const_array(mfi);
                    intg = m_linop->m_integral[lev]->const_array(mfi);
                }
            }
            auto skip_cell = [=] AMREX_GPU_DEVICE (IntVect const& c) -> bool
            {
                return eb && flag(c).isCovered();
            };
            auto cell_matrix = [=] AMREX_GPU_DEVICE (IntVect const& c,
                                                     GpuArray<Real,AMREX_SPACEDIM> const& s,
                                                     Real (&K)[nnodes_cell][nnodes_cell])
            {
                constexpr Real almostone = Real(1.) - Real(100.)*std::numeric_limits<Real>::epsilon();
                if (eb && !flag(c).isRegular() && vfrac(c) < almostone) {
                    element_matrix_eb(dx, s[0], intg, vfrac(c), c, K);
                } else {
                    element_matrix(dx, s, K);
                }
            };
#else
            auto skip_cell = [=] AMREX_GPU_DEVICE (IntVect const&) -> bool { return false; };
            auto cell_matrix = [=] AMREX_GPU_DEVICE (IntVect const&,
                                                     GpuArray<Real,AMREX_SPACEDIM> const& s,
                                                     Real (&K)[nnodes_cell][nnodes_cell])
            {
                element_matrix(dx, s, K);
            };
#endif

            auto on_dirichlet = [=] AMREX_GPU_DEVICE (IntVect const& n) -> bool
            {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if ((dlo[d] && n[d] == domain.smallEnd(d)) ||
                        (dhi[d] && n[d] == domain.bigEnd(d)+1)) { return true; }
                }
                return false;
            };

            // Dirichlet domain nodes (fine index, weight) of a corner node: the
            // node itself, or the coarse Dirichlet nodes that a hanging node
            // interpolates from (qmap keeps only the unknowns), as fine nodes
            // at the same place, inside this box since fine boxes are
            // coarsenable
            auto dmap_nodes = [=] AMREX_GPU_DEVICE (IntVect const& n, IntVect* dn, Real* dw) -> int
            {
                if (on_dirichlet(n)) {
                    dn[0] = n;
                    dw[0] = Real(1.0);
                    return 1;
                }
                if (rr == 1 || gid(n) >= 0 || !dm(n)) { return 0; }
                IntVect cn[nnodes_cell];
                Real cw[nnodes_cell];
                const int m = interp_stencil(n, rr, cn, cw);
                int nd = 0;
                for (int k = 0; k < m; ++k) {
                    const IntVect fn = cn[k]*rr;
                    if (cgid(cn[k]) < 0 && on_dirichlet(fn)) {
                        dn[nd] = fn;
                        dw[nd] = cw[k];
                        ++nd;
                    }
                }
                return nd;
            };

            // Composite unknowns (id, weight) of a corner node
            auto qmap = [=] AMREX_GPU_DEVICE (IntVect const& n, Long* g, Real* w,
                                              bool scaled) -> int
            {
                if (gid(n) >= 0) {
                    g[0] = gid(n);
                    w[0] = (sc && scaled) ? sc(n) : Real(1.0);
                    return 1;
                }
                if (rr == 1 || !dm(n) || on_dirichlet(n)) { return 0; }
                IntVect cn[nnodes_cell];
                Real cw[nnodes_cell];
                const int m = interp_stencil(n, rr, cn, cw);
                int nq = 0;
                for (int k = 0; k < m; ++k) {
                    const Long cg = cgid(cn[k]);
                    if (cg >= 0) {
                        g[nq] = cg;
                        w[nq] = (csc && scaled) ? cw[k]*csc(cn[k]) : cw[k];
                        ++nq;
                    }
                }
                return nq;
            };

            // Pass 1: number of triplets per cell
            const Long ncells = ccbx.numPts();
            Gpu::DeviceVector<Long> count(ncells), offset(ncells);
            Long* pcount = count.data();
            const auto lo = amrex::lbound(ccbx);
            const auto len = amrex::length(ccbx);
            amrex::ParallelFor(ccbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Long idx = (Long(k-lo.z)*len.y + (j-lo.y))*len.x + (i-lo.x);
                const IntVect cell(AMREX_D_DECL(i,j,k));
                if ((has_fine && ccf(i,j,k)) || skip_cell(cell)) { pcount[idx] = 0; return; }
                int ntot = 0;
                for (int a = 0; a < nnodes_cell; ++a) {
                    IntVect n = cell;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { n[d] += (a >> d) & 1; }
                    Long g[nnodes_cell];
                    Real w[nnodes_cell];
                    ntot += qmap(n, g, w, true);
                }
                pcount[idx] = Long(ntot)*ntot;
            });
            const Long ntrip = Scan::ExclusiveSum(ncells, pcount, offset.data(), Scan::retSum);
            if (ntrip == 0) { Gpu::streamSynchronize(); continue; }

            // Pass 2: the triplets
            Triplets& t = chunks.emplace_back();
            t.rows.resize(ntrip); t.cols.resize(ntrip); t.vals.resize(ntrip);
            BigInt* prow = t.rows.data();
            BigInt* pcol = t.cols.data();
            Real* pval = t.vals.data();
            Long const* poff = offset.data();
            amrex::ParallelFor(ccbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const IntVect cell(AMREX_D_DECL(i,j,k));
                if ((has_fine && ccf(i,j,k)) || skip_cell(cell)) { return; }
                const Long idx = (Long(k-lo.z)*len.y + (j-lo.y))*len.x + (i-lo.x);
                GpuArray<Real,AMREX_SPACEDIM> s;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    s[d] = (const_sigma != Real(0.0)) ? const_sigma
                         : (per_dir ? sig[d](cell) : sig[0](cell));
                }
                Real K[nnodes_cell][nnodes_cell];
                cell_matrix(cell, s, K);
                Long g[nnodes_cell][nnodes_cell];
                Real w[nnodes_cell][nnodes_cell];
                int nq[nnodes_cell];
                for (int a = 0; a < nnodes_cell; ++a) {
                    IntVect n = cell;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { n[d] += (a >> d) & 1; }
                    nq[a] = qmap(n, g[a], w[a], true);
                }
                Long e = poff[idx];
                for (int a = 0; a < nnodes_cell; ++a) {
                    for (int qa = 0; qa < nq[a]; ++qa) {
                        for (int b = 0; b < nnodes_cell; ++b) {
                            for (int qb = 0; qb < nq[b]; ++qb) {
                                prow[e] = BigInt(g[a][qa]);
                                pcol[e] = BigInt(g[b][qb]);
                                // pinned unknown: identity row (added below),
                                // zero column (x_pin = 0 keeps S symmetric)
                                pval[e] = (g[a][qa] == pin || g[b][qb] == pin) ? Real(0.0)
                                        : w[a][qa]*w[b][qb]*K[a][b];
                                ++e;
                            }
                        }
                    }
                }
            });
            // sync: the 'continue's below free count/offset
            Gpu::streamSynchronize();

            // Couplings of the unknowns to Dirichlet domain nodes, directly or
            // through the interpolation of hanging nodes, whose values (taken
            // from the initial guess at each solve) move to the right-hand
            // side.  Unscaled; loadRHS scales the right-hand side.
            bool touches_dirichlet = false;
            {
                Box const ndbx = mfi.validbox();
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if ((dlo[d] && ndbx.smallEnd(d) == domain.smallEnd(d)) ||
                        (dhi[d] && ndbx.bigEnd(d) == domain.bigEnd(d)+1)) { touches_dirichlet = true; }
                }
            }
            if (!touches_dirichlet) { continue; }
            amrex::ParallelFor(ccbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Long idx = (Long(k-lo.z)*len.y + (j-lo.y))*len.x + (i-lo.x);
                const IntVect cell(AMREX_D_DECL(i,j,k));
                pcount[idx] = 0;
                if ((has_fine && ccf(i,j,k)) || skip_cell(cell)) { return; }
                int nd = 0, nu = 0;
                for (int a = 0; a < nnodes_cell; ++a) {
                    IntVect n = cell;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { n[d] += (a >> d) & 1; }
                    IntVect dn[nnodes_cell];
                    Real dw[nnodes_cell];
                    nd += dmap_nodes(n, dn, dw);
                    Long g[nnodes_cell];
                    Real w[nnodes_cell];
                    nu += qmap(n, g, w, false);
                }
                pcount[idx] = Long(nd)*nu;
            });
            const Long nbc = Scan::ExclusiveSum(ncells, pcount, offset.data(), Scan::retSum);
            if (nbc == 0) { Gpu::streamSynchronize(); continue; }
            Gpu::DeviceVector<Long> bc_row(nbc);
            Gpu::DeviceVector<Real> bc_coef(nbc);
            Gpu::DeviceVector<IntVect> bc_node(nbc);
            Long* pbr = bc_row.data();
            Real* pbc = bc_coef.data();
            IntVect* pbn = bc_node.data();
            amrex::ParallelFor(ccbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const IntVect cell(AMREX_D_DECL(i,j,k));
                if ((has_fine && ccf(i,j,k)) || skip_cell(cell)) { return; }
                const Long idx = (Long(k-lo.z)*len.y + (j-lo.y))*len.x + (i-lo.x);
                if (pcount[idx] == 0) { return; }
                GpuArray<Real,AMREX_SPACEDIM> s;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    s[d] = (const_sigma != Real(0.0)) ? const_sigma
                         : (per_dir ? sig[d](cell) : sig[0](cell));
                }
                Real K[nnodes_cell][nnodes_cell];
                cell_matrix(cell, s, K);
                Long e = poff[idx];
                for (int b = 0; b < nnodes_cell; ++b) {
                    IntVect nb = cell;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { nb[d] += (b >> d) & 1; }
                    IntVect dn[nnodes_cell];
                    Real dw[nnodes_cell];
                    const int nd = dmap_nodes(nb, dn, dw);
                    if (nd == 0) { continue; }
                    for (int a = 0; a < nnodes_cell; ++a) {
                        IntVect na = cell;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) { na[d] += (a >> d) & 1; }
                        Long g[nnodes_cell];
                        Real w[nnodes_cell];
                        const int nq = qmap(na, g, w, false);
                        for (int q = 0; q < nq; ++q) {
                            for (int p = 0; p < nd; ++p) {
                                pbr[e] = g[q];
                                pbc[e] = (g[q] == pin) ? Real(0.0) : w[q]*K[a][b]*dw[p];
                                pbn[e] = dn[p];
                                ++e;
                            }
                        }
                    }
                }
            });
            Vector<Long> h_row(nbc);
            Vector<Real> h_coef(nbc);
            Vector<IntVect> h_node(nbc);
            Gpu::copyAsync(Gpu::deviceToHost, bc_row.begin(), bc_row.end(), h_row.begin());
            Gpu::copyAsync(Gpu::deviceToHost, bc_coef.begin(), bc_coef.end(), h_coef.begin());
            Gpu::copyAsync(Gpu::deviceToHost, bc_node.begin(), bc_node.end(), h_node.begin());
            Gpu::streamSynchronize();
            for (Long n = 0; n < nbc; ++n) {
                AMREX_ALWAYS_ASSERT(mfi.validbox().contains(h_node[n]));
                if (h_coef[n] != Real(0.0)) {
                    m_dbc.push_back(DirichletCoupling{.row = h_row[n], .coef = h_coef[n], .lev = lev,
                                                      .li = mfi.LocalIndex(), .node = h_node[n]});
                }
            }
        }
        // gid_d/cgid_d are freed at the end of this iteration
        Gpu::streamSynchronize();
        flush();
    }

    const bool pin_here = (m_pin_gid >= m_ilower && m_pin_gid <= m_iupper);
    if (m_backend == Backend::hypre) {
        if (pin_here) {
            Gpu::DeviceVector<BigInt> r(1, BigInt(m_pin_gid));
            Gpu::DeviceVector<Real> v(1, Real(1.0));
            Gpu::DeviceVector<HYPRE_Int> o(1, 1);
            HYPRE_IJMatrixAddToValues(m_A, 1, o.data(), r.data(), r.data(), v.data());
            Gpu::hypreSynchronize();
        }
        HYPRE_IJMatrixAssemble(m_A);
    }

    if (need_csr) {
        BL_PROFILE("HypreMLNodeLap::buildCSR()");
        if (pin_here) {
            auto& t = host_trips.emplace_back();
            t.row = m_pin_gid;
            t.col = m_pin_gid;
            t.val = Real(1.0);
        }
        host_trips = composite::exchange_triplets(std::move(host_trips), m_rank_offset, m_comm);
        m_csr.define(m_ilower, m_iupper - m_ilower + 1, std::move(host_trips));
    }

    if (m_backend == Backend::algmg) {
        m_spA = std::make_unique<SpMatrix<Real>>();
        if (algmg_dev) {
            if (pin_here) {
                const auto n0 = Long(alg_rows.size());
                alg_rows.resize(n0+1);
                alg_cols.resize(n0+1);
                alg_vals.resize(n0+1);
                auto* ar = alg_rows.data();
                auto* ac = alg_cols.data();
                auto* av = alg_vals.data();
                amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
                {
                    ar[n0] = pin;
                    ac[n0] = pin;
                    av[n0] = Real(1.0);
                });
            }
            composite::triplets_to_spmatrix(alg_rows, alg_cols, alg_vals, m_ilower,
                                            m_iupper - m_ilower + 1, m_rank_offset, m_comm, *m_spA);
        } else {
            m_csr.toSpMatrix(*m_spA, m_rank_offset);
        }
        m_algmg = std::make_unique<AlgMG<Real>>(*m_spA);
        m_algmg->setVerbose(m_verbose > 1 ? m_verbose-1 : 0);
        // V-cycles alone can stall on composite AMR systems; S is symmetric
        composite::algmg_parse_options(*m_algmg, m_prefix, AlgMGKrylovSolver::pcg);
    }
    if (!composite::dump_pending(m_dump_prefix, m_dump_max)) { m_csr = composite::OwnedCSR{}; }
}

void HypreMLNodeLap::loadRHS (Vector<MultiFab const*> const& a_rhs,
                              Vector<MultiFab*> const& a_sol)
{
    BL_PROFILE("HypreMLNodeLap::loadRHS()");

    // Owned unknowns are exactly this rank's rows; contributions to rows of
    // other ranks (coarse/fine and Dirichlet couplings) are sent to them.
    const Long nrows = m_iupper - m_ilower + 1;
    Vector<Real> bvec(nrows, Real(0.0));
    Vector<Real> mvec(nrows, Real(0.0));   // solvability weights
    Vector<Real> svec(nrows, Real(1.0));   // unknown scaling
    Vector<composite::Triplet> remote;     // (row, unused, value)
    auto add = [&] (Long row, Real v)
    {
        if (row >= m_ilower && row <= m_iupper) {
            bvec[row - m_ilower] += v;
        } else {
            auto& t = remote.emplace_back();
            t.row = row;
            t.col = 0;
            t.val = v;
        }
    };

    for (int lev = 0; lev < m_nlevels; ++lev) {
        const auto dx = m_geom[lev].CellSizeArray();
        const Real hd = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
        auto rhs = host_copy(*a_rhs[lev], 0);
        const int rr = (lev > 0) ? m_ref_ratio[lev-1] : 1;
        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->const_array(mfi);
            auto const& own = m_owner[lev]->const_array(mfi);
            auto const& b = rhs->const_array(mfi);
            Array4<Long const> cgid;
            Array4<int const> cndf;
            if (lev > 0) {
                cgid = m_cgid_on_f[lev]->const_array(mfi);
                cndf = m_cndfine_on_f[lev]->const_array(mfi);
            }
            Array4<Real const> ndv;
            if (m_ndvol[lev]) { ndv = m_ndvol[lev]->const_array(mfi); }
            Array4<Real const> sc;
            if (m_scale[lev]) { sc = m_scale[lev]->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, [[maybe_unused]] int k)
            {
                const IntVect n(AMREX_D_DECL(i,j,k));
                if (!own(n) || gid(n) < 0) { return; }
                const Long r = gid(n) - m_ilower;
                const Real M = hd / Real(neumannScale(lev, n));
                // MLMG row p times -M_p is the finite-element row
                bvec[r] -= M*b(n);
                mvec[r] = ndv ? hd*ndv(n) : M;   // fluid volume with EB
                if (sc) { svec[r] = sc(n); }
                if (lev > 0) {
                    // Undo the fine equations that MLMG adds into the coarse/fine
                    // rows: add P_nI M_n b_n to each coarse/fine unknown I with
                    // |n - rr I|_inf <= rr-1.
                    IntVect clo, chi;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        clo[d] = amrex::coarsen(n[d], rr);
                        chi[d] = (n[d] == clo[d]*rr) ? clo[d] : clo[d]+1;
                    }
                    amrex::LoopOnCpu(Box(clo,chi), [&] (int ci, int cj, [[maybe_unused]] int ck)
                    {
                        const IntVect I(AMREX_D_DECL(ci,cj,ck));
                        if (cndf(I) != nodelap_detail::crse_fine_node || cgid(I) < 0) { return; }
                        Real P = 1.0;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            P *= Real(1.0) - Real(std::abs(n[d] - I[d]*rr))/Real(rr);
                        }
                        if (P > Real(0.0)) { add(cgid(I), P*M*b(n)); }
                    });
                }
            });
        }
    }

    // Dirichlet domain values (from the initial guess): K_UD x_D to the
    // right-hand side, -M A x = M rhs  =>  K_UU x_U = -M rhs - K_UD x_D
    if (!m_dbc.empty()) {
        Vector<std::unique_ptr<MultiFab>> hsol(m_nlevels);
        for (auto const& c : m_dbc) {
            if (!hsol[c.lev]) { hsol[c.lev] = host_copy(*a_sol[c.lev], 0); }
            auto const& x = hsol[c.lev]->atLocalIdx(c.li).const_array();
            add(c.row, -c.coef * x(c.node));
        }
    }

    for (auto const& t : composite::exchange_triplets(std::move(remote), m_rank_offset, m_comm)) {
        bvec[t.row - m_ilower] += t.val;
    }

    // Singular problems: make the right-hand side exactly consistent
    if (m_linop->isBottomSingular()) {
        Real bsum = 0.0, msum = 0.0;
        for (Long r = 0; r < nrows; ++r) { bsum += bvec[r]; msum += mvec[r]; }
        ParallelAllReduce::Sum(bsum, m_comm);
        ParallelAllReduce::Sum(msum, m_comm);
        const Real shift = bsum / msum;
        for (Long r = 0; r < nrows; ++r) { bvec[r] -= shift*mvec[r]; }
        if (m_pin_gid >= m_ilower && m_pin_gid <= m_iupper) {
            bvec[m_pin_gid - m_ilower] = 0.0;
        }
    }

    // Scaled system: b' = s b
    Real bnorm2 = 0.0;
    for (Long r = 0; r < nrows; ++r) {
        bvec[r] *= svec[r];
        bnorm2 += bvec[r]*bvec[r];
    }
    ParallelAllReduce::Sum(bnorm2, m_comm);
    m_bnorm = std::sqrt(bnorm2);

    m_bvec = std::move(bvec);
    m_xvec.assign(nrows, Real(0.0));
}

void HypreMLNodeLap::getSolution (Vector<MultiFab*> const& a_sol)
{
    BL_PROFILE("HypreMLNodeLap::getSolution()");

    Vector<Real> const& xvec = m_xvec;

    // Per level: owned unknowns, then shared copies (ParallelAdd of a sparse field)
    Vector<std::unique_ptr<MultiFab>> hsol(m_nlevels);
    for (int lev = 0; lev < m_nlevels; ++lev) {
        MultiFab tmp(m_nd_grids[lev], m_dmap[lev], 1, 0, host_info());
        tmp.setVal(0.0);
        Gpu::streamSynchronize();
        for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->const_array(mfi);
            auto const& own = m_owner[lev]->const_array(mfi);
            auto const& t = tmp.array(mfi);
            Array4<Real const> sc;
            if (m_scale[lev]) { sc = m_scale[lev]->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                if (own(i,j,k) && gid(i,j,k) >= 0) {
                    // x = s x'
                    t(i,j,k) = xvec[gid(i,j,k) - m_ilower] * (sc ? sc(i,j,k) : Real(1.0));
                }
            });
        }
        hsol[lev] = std::make_unique<MultiFab>(m_nd_grids[lev], m_dmap[lev], 1, 0, host_info());
        hsol[lev]->setVal(0.0);
        hsol[lev]->ParallelAdd(tmp, 0, 0, 1, m_geom[lev].periodicity());
        Gpu::streamSynchronize();

        // Dirichlet domain nodes keep their given values
        if (m_has_dirichlet) {
            auto orig = host_copy(*a_sol[lev], 0);
            for (MFIter mfi(*hsol[lev]); mfi.isValid(); ++mfi) {
                auto const& s = hsol[lev]->array(mfi);
                auto const& o = orig->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, [[maybe_unused]] int k)
                {
                    if (onDirichletFace(lev, IntVect(AMREX_D_DECL(i,j,k)))) { s(i,j,k) = o(i,j,k); }
                });
            }
        }

        // Coarse/fine nodes of this level: interpolate the coarser solution
        if (lev > 0) {
            const int rr = m_ref_ratio[lev-1];
            const BoxArray cba = amrex::coarsen(m_nd_grids[lev], rr);
            MultiFab csol(cba, m_dmap[lev], 1, 0, host_info());
            csol.setVal(0.0);
            csol.ParallelCopy(*hsol[lev-1], 0, 0, 1, IntVect(0), IntVect(0),
                              m_geom[lev-1].periodicity());
            Gpu::streamSynchronize();
            for (MFIter mfi(*hsol[lev]); mfi.isValid(); ++mfi) {
                auto const& gid = m_gid[lev]->const_array(mfi);
                auto const& dm = m_dmsk[lev]->const_array(mfi);
                auto const& s = hsol[lev]->array(mfi);
                auto const& cs = csol.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, [[maybe_unused]] int k)
                {
                    const IntVect n(AMREX_D_DECL(i,j,k));
                    if (gid(n) >= 0 || !dm(n) || onDirichletFace(lev, n)) { return; }
                    IntVect cn[nnodes_cell];
                    Real cw[nnodes_cell];
                    const int m = interp_stencil(n, rr, cn, cw);
                    Real v = 0.0;
                    for (int q = 0; q < m; ++q) { v += cw[q]*cs(cn[q]); }
                    s(n) = v;
                });
            }
        }
    }

    for (int lev = 0; lev < m_nlevels; ++lev) {
        a_sol[lev]->setVal(0.0);
        MultiFab::Copy(*a_sol[lev], *hsol[lev], 0, 0, 1, 0);
    }
    Gpu::streamSynchronize();   // hsol is freed on return

    // Covered coarse nodes: inject the finer solution (finest to coarsest)
    for (int lev = m_nlevels-1; lev > 0; --lev) {
        amrex::average_down_nodal(*a_sol[lev], *a_sol[lev-1], IntVect(m_ref_ratio[lev-1]));
    }
    for (int lev = 0; lev < m_nlevels; ++lev) {
        a_sol[lev]->FillBoundary(m_geom[lev].periodicity());
    }
}

void HypreMLNodeLap::fail (std::string const& msg) const
{
    if (m_throw_exception) { throw std::runtime_error(msg); }
    amrex::Abort(msg);
}

void HypreMLNodeLap::solve (Vector<MultiFab*> const& a_sol, Vector<MultiFab const*> const& a_rhs,
                            Real a_reltol, Real a_abstol, int a_maxiter)
{
    BL_PROFILE("HypreMLNodeLap::solve()");
    AMREX_ALWAYS_ASSERT(static_cast<int>(a_sol.size()) == m_nlevels &&
                        static_cast<int>(a_rhs.size()) == m_nlevels);
    // The unknowns of a box are found by its index
    for (int lev = 0; lev < m_nlevels; ++lev) {
        for (auto const* mf : {a_rhs[lev], static_cast<MultiFab const*>(a_sol[lev])}) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf->boxArray() == m_nd_grids[lev] &&
                                             mf->DistributionMap() == m_dmap[lev],
                "HypreMLNodeLap: rhs and solution must have the grids of the operator");
        }
    }

    loadRHS(a_rhs, a_sol);
    // BoomerAMG can report convergence with non-finite data.  Only the
    // assembled right-hand side counts (MLMG allows NaNs at Dirichlet nodes);
    // it includes the Dirichlet values, and m_bnorm is reduced over m_comm.
    if (!std::isfinite(m_bnorm)) {
        fail("HypreMLNodeLap: right-hand side or Dirichlet values are not finite");
    }
    const Long nrows = m_iupper - m_ilower + 1;
    const int dump_k = (m_dump_max > 0) ? composite::next_dump_index(m_dump_prefix, m_dump_max) : -1;
    if (dump_k >= 0) { dumpSystem(dump_k, false); }

    // MLMG rows = FE rows / (Neumann-reduced) nodal volume >= h^D/2^D, and
    // the scaling s >= 1 only enlarges rows: this bounds every MLMG row.
    const Real fine_vol = AMREX_D_TERM(m_geom.back().CellSize(0),
                                       *m_geom.back().CellSize(1),
                                       *m_geom.back().CellSize(2));
    const Real abstol = (a_abstol > Real(0.0))
        ? a_abstol * fine_vol / Real(AMREX_D_TERM(2,*2,*2)) : Real(0.0);

    if (m_backend == Backend::hypre) {
        HYPRE_IJVectorInitialize(m_b);
        HYPRE_IJVectorInitialize(m_x);
        if (nrows > 0) {
            Vector<BigInt> rowids(nrows);
            std::iota(rowids.begin(), rowids.end(), static_cast<BigInt>(m_ilower));
            Gpu::DeviceVector<BigInt> d_rows(nrows);
            Gpu::DeviceVector<Real> d_b(nrows), d_x(nrows);
            Gpu::copyAsync(Gpu::hostToDevice, rowids.begin(), rowids.end(), d_rows.begin());
            Gpu::copyAsync(Gpu::hostToDevice, m_bvec.begin(), m_bvec.end(), d_b.begin());
            Gpu::copyAsync(Gpu::hostToDevice, m_xvec.begin(), m_xvec.end(), d_x.begin());
            Gpu::streamSynchronize();
            HYPRE_IJVectorSetValues(m_b, static_cast<HYPRE_Int>(nrows), d_rows.data(), d_b.data());
            HYPRE_IJVectorSetValues(m_x, static_cast<HYPRE_Int>(nrows), d_rows.data(), d_x.data());
            Gpu::hypreSynchronize();
        }
        HYPRE_IJVectorAssemble(m_b);
        HYPRE_IJVectorAssemble(m_x);

        m_ij->setThrowException(m_throw_exception);
        m_ij->solve(a_reltol, abstol, a_maxiter);
        m_num_iters = static_cast<int>(m_ij->getNumIters());
        m_final_res = m_ij->getFinalResidualNorm();

        if (nrows > 0) {
            Vector<BigInt> rowids(nrows);
            std::iota(rowids.begin(), rowids.end(), static_cast<BigInt>(m_ilower));
            Gpu::DeviceVector<BigInt> d_rows(nrows);
            Gpu::DeviceVector<Real> d_x(nrows);
            Gpu::copyAsync(Gpu::hostToDevice, rowids.begin(), rowids.end(), d_rows.begin());
            Gpu::streamSynchronize();
            HYPRE_IJVectorGetValues(m_x, static_cast<HYPRE_Int>(nrows), d_rows.data(), d_x.data());
            Gpu::hypreSynchronize();
            Gpu::copyAsync(Gpu::deviceToHost, d_x.begin(), d_x.end(), m_xvec.begin());
            Gpu::streamSynchronize();
        }
    } else {
        BL_PROFILE("HypreMLNodeLap::solve_algmg");
        AlgPartition part(m_rank_offset);
        AlgVector<Real> b(part), x(part);
        Gpu::copyAsync(Gpu::hostToDevice, m_bvec.begin(), m_bvec.end(), b.data());
        Gpu::copyAsync(Gpu::hostToDevice, m_xvec.begin(), m_xvec.end(), x.data());
        Gpu::streamSynchronize();
        m_algmg->setThrowException(m_throw_exception);
        m_algmg->setRelTol(a_reltol);
        m_algmg->setAbsTol(abstol);
        m_algmg->setMaxIter(a_maxiter);
        m_algmg->solve(x, b);
        m_num_iters = m_algmg->getNumIters();
        m_final_res = (m_bnorm > Real(0.0)) ? m_algmg->getResidualNorm() / m_bnorm : Real(0.0);
        Gpu::copyAsync(Gpu::deviceToHost, x.data(), x.data()+nrows, m_xvec.begin());
        Gpu::streamSynchronize();
    }

    if (m_verbose > 0) {
        amrex::Print() << "HypreMLNodeLap(" << ((m_backend == Backend::hypre) ? "hypre" : "algmg")
                       << "): " << m_num_iters
                       << " iterations, relative residual " << m_final_res << '\n';
    }

    const bool converged = (m_final_res <= a_reltol) ||
        (abstol > Real(0.0) && m_final_res*m_bnorm <= abstol);
    if (!converged) {
        fail("HypreMLNodeLap: failed to converge after " + std::to_string(m_num_iters) +
             " iterations; final relative residual " + composite::to_sci(m_final_res));
    }

    if (dump_k >= 0) { dumpSystem(dump_k, true); }

    getSolution(a_sol);
}

void HypreMLNodeLap::dumpSystem (int k, bool solution)
{
    BL_PROFILE("HypreMLNodeLap::dumpSystem()");
    const std::string f = m_dump_prefix + "_" + std::to_string(k);
    if (solution) {
        composite::write_vector_market(f + "_x.mtx", m_xvec, m_nglobal, m_ilower, m_comm,
                                       m_dump_gather_max);
        return;
    }
    composite::write_matrix_market(f + "_A.mtx", m_csr, m_nglobal, m_comm, m_dump_gather_max);
    composite::write_vector_market(f + "_b.mtx", m_bvec, m_nglobal, m_ilower, m_comm,
                                   m_dump_gather_max);

    // Row map: level and node index of each row, and its scaling
    const Long nrows = m_iupper - m_ilower + 1;
    Vector<std::array<Long,4>> info(nrows, std::array<Long,4>{-1,-1,-1,-1});
    Vector<Real> scl(nrows, Real(1.0));
    for (int lev = 0; lev < m_nlevels; ++lev) {
        for (MFIter mfi(*m_gid[lev]); mfi.isValid(); ++mfi) {
            auto const& gid = m_gid[lev]->const_array(mfi);
            auto const& own = m_owner[lev]->const_array(mfi);
            Array4<Real const> sc;
            if (m_scale[lev]) { sc = m_scale[lev]->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int kk)
            {
                if (!own(i,j,kk) || gid(i,j,kk) < 0) { return; }
                const Long r = gid(i,j,kk) - m_ilower;
                info[r] = std::array<Long,4>{lev, i, j, kk};
                if (sc) { scl[r] = sc(i,j,kk); }
            });
        }
    }
    composite::write_row_map(f + "_rows.txt",
                             "HypreMLNodeLap composite system: " + std::to_string(m_nglobal) +
                             " rows on " + std::to_string(m_nlevels) + " level(s), " +
                             std::to_string(AMREX_SPACEDIM) + "D nodal",
                             info, scl, m_ilower, m_comm, m_dump_gather_max);
}

}
