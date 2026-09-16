//
// State redistribution regression test.
//
// State redistribution rests on one invariant: every neighbor that MakeITracker
// records for a cell is a cell whose volume MakeStateRedistUtils folds into
// nbhd_vol/alpha/cent_hat *and* whose state StateRedistribute folds into Qhat.
// If a neighborhood reaches outside a non-periodic domain boundary, or reaches a
// periodic image that Qhat leaves out, the volume is counted but the state is not,
// and conservation is silently lost.
//
// This test checks the invariant two ways:
//
//   1. Directly, by calling MakeITracker on hand-built area and volume fractions
//      as well as on real cut-cell geometry, and verifying that no neighborhood
//      reaches outside a non-periodic domain.
//
//   2. End to end, by redistributing a *constant* state with a zero update.  The
//      answer must come back unchanged to roundoff; a cell whose volume is counted
//      in nbhd_vol but whose state is left out of Qhat shows up as a nonzero update.
//
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <limits>

using namespace amrex;

namespace {

#if (AMREX_SPACEDIM == 2)
constexpr int itracker_ncomp = 4;
constexpr int nbor_id_max = 9;
#else
constexpr int itracker_ncomp = 8;
constexpr int nbor_id_max = 27;
#endif

constexpr Real target_volfrac = Real(0.5);

struct NborMaps
{
    GpuArray<int,nbor_id_max> imap;
    GpuArray<int,nbor_id_max> jmap;
    GpuArray<int,nbor_id_max> kmap;
};

NborMaps make_maps ()
{
#if (AMREX_SPACEDIM == 2)
    return NborMaps{.imap = {0,-1, 0, 1,-1, 1,-1, 0, 1},
                    .jmap = {0,-1,-1,-1, 0, 0, 1, 1, 1},
                    .kmap = {0, 0, 0, 0, 0, 0, 0, 0, 0}};
#else
    return NborMaps{.imap = {0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1},
                    .jmap = {0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1},
                    .kmap = {0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
#endif
}

// Number of neighbors of (i,j,k) that fall outside a non-periodic domain boundary.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int count_outside (Array4<int const> const& itracker, int i, int j, int k,
                   Box const& domain, GpuArray<int,AMREX_SPACEDIM> const& is_per,
                   NborMaps const& m)
{
    int nbad = 0;
    for (int n = 1; n <= itracker(i,j,k,0); ++n) {
        const int idx = itracker(i,j,k,n);
        const IntVect nbor(AMREX_D_DECL(i + m.imap[idx],
                                        j + m.jmap[idx],
                                        k + m.kmap[idx]));
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!is_per[d] && (nbor[d] < domain.smallEnd(d) || nbor[d] > domain.bigEnd(d))) {
                ++nbad;
                break;
            }
        }
    }
    return nbad;
}

//
// Test 1: hand-built fractions.
//
// Make one cell the only small cell in an otherwise regular domain and give it a
// prescribed EB normal.  This pins down the domain-edge cases -- in particular a
// cell on a domain face or corner whose normal points out of the domain -- which
// are awkward to hit with a real implicit function.
//
int synthetic_case (int n_cell, IntVect const& cell, Real vfrac_of_cell,
                    RealArray const& normal, std::string const& label)
{
    const Box domain(IntVect(0), IntVect(n_cell-1));
    const RealBox rb({AMREX_D_DECL(Real(0),Real(0),Real(0))},
                     {AMREX_D_DECL(Real(1),Real(1),Real(1))});
    const Array<int,AMREX_SPACEDIM> is_per{AMREX_D_DECL(0,0,0)};
    const Geometry geom(domain, rb, 0, is_per);

    // MakeITracker looks 4 cells out from the box it is given and reads area
    // fractions one further out than that.
    const Box gbx = amrex::grow(domain,6);

    FArrayBox vfrac_fab(gbx, 1, The_Async_Arena());
    vfrac_fab.setVal<RunOn::Device>(Real(1));

    AMREX_D_TERM(FArrayBox apx_fab(amrex::surroundingNodes(gbx,0), 1, The_Async_Arena());,
                 FArrayBox apy_fab(amrex::surroundingNodes(gbx,1), 1, The_Async_Arena());,
                 FArrayBox apz_fab(amrex::surroundingNodes(gbx,2), 1, The_Async_Arena()););
    AMREX_D_TERM(apx_fab.setVal<RunOn::Device>(Real(1));,
                 apy_fab.setVal<RunOn::Device>(Real(1));,
                 apz_fab.setVal<RunOn::Device>(Real(1)););

    // MakeITracker builds the normal from the jump in area fraction across the cell,
    // so setting the two opposing faces of one cell is enough to prescribe it.
    const Real base = Real(0.2);
    const GpuArray<Real,AMREX_SPACEDIM> nrm{AMREX_D_DECL(normal[0],normal[1],normal[2])};
    auto const& vf = vfrac_fab.array();
    AMREX_D_TERM(auto const& ax = apx_fab.array();,
                 auto const& ay = apy_fab.array();,
                 auto const& az = apz_fab.array(););
    amrex::ParallelFor(Box(cell,cell),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vf(i,j,k) = vfrac_of_cell;
        AMREX_D_TERM(ax(i,j,k) = base; ax(i+1,j,k) = base + nrm[0];,
                     ay(i,j,k) = base; ay(i,j+1,k) = base + nrm[1];,
                     az(i,j,k) = base; az(i,j,k+1) = base + nrm[2];);
    });

    IArrayBox itracker_fab(gbx, itracker_ncomp, The_Async_Arena());
    itracker_fab.setVal<RunOn::Device>(0);

    MakeITracker(domain, AMREX_D_DECL(apx_fab.const_array(),apy_fab.const_array(),apz_fab.const_array()),
                 vfrac_fab.const_array(), itracker_fab.array(), geom, target_volfrac);

    // Bring the one cell we care about back to the host.
    IArrayBox host_fab(Box(cell,cell), itracker_ncomp, The_Pinned_Arena());
    host_fab.copy<RunOn::Device>(itracker_fab, Box(cell,cell));
    Gpu::streamSynchronize();

    const auto maps = make_maps();
    const GpuArray<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,0,0)};
    const int ci = cell[0];
    const int cj = cell[1];
#if (AMREX_SPACEDIM == 3)
    const int ck = cell[2];
#else
    const int ck = 0;
#endif
    auto const& host_itracker = host_fab.const_array();
    const int nbad = count_outside(host_itracker, ci, cj, ck, domain, isp, maps);
    const int nnbor = host_itracker(ci,cj,ck,0);

    amrex::Print() << "  " << label << " : " << nnbor << " neighbors, "
                   << nbad << " outside the domain" << (nbad > 0 ? "   ** FAILED **" : "") << '\n';
    return (nbad > 0) ? 1 : 0;
}

//
// Test 2: cut-cell geometry.
//
// Cut the domain with a plane, check the same invariant everywhere, then
// redistribute a constant state and check that it comes back unchanged.
//
int geometry_case (int n_cell, Array<int,AMREX_SPACEDIM> const& is_per,
                   RealArray const& point, RealArray const& normal,
                   int max_grid_size, std::string const& label)
{
    const Box domain(IntVect(0), IntVect(n_cell-1));
    const RealBox rb({AMREX_D_DECL(Real(0),Real(0),Real(0))},
                     {AMREX_D_DECL(Real(1),Real(1),Real(1))});
    const Geometry geom(domain, rb, 0, is_per);

    EB2::Build(EB2::makeShop(EB2::PlaneIF(point,normal)), geom, 0, 0);

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    const DistributionMapping dm(ba);

    EBFArrayBoxFactory factory(EB2::IndexSpace::top().getLevel(geom), geom, ba, dm,
                               {5,5,5}, EBSupport::full);

    const int ncomp = 1;
    MultiFab U_in    (ba, dm, ncomp, 5, MFInfo(), factory);
    MultiFab dUdt_in (ba, dm, ncomp, 4, MFInfo(), factory);
    MultiFab dUdt_out(ba, dm, ncomp, 0, MFInfo(), factory);
    MultiFab scratch (ba, dm, ncomp, 4, MFInfo(), factory);

    U_in.setVal(Real(1));
    dUdt_in.setVal(Real(0));
    dUdt_out.setVal(Real(0));
    scratch.setVal(Real(0));
    U_in.FillBoundary(geom.periodicity());

    Vector<BCRec> bcrec(ncomp);
    for (int n = 0; n < ncomp; ++n) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bcrec[n].setLo(d, is_per[d] ? BCType::int_dir : BCType::foextrap);
            bcrec[n].setHi(d, is_per[d] ? BCType::int_dir : BCType::foextrap);
        }
    }
    Gpu::DeviceVector<BCRec> d_bcrec(ncomp);
    Gpu::copy(Gpu::hostToDevice, bcrec.begin(), bcrec.end(), d_bcrec.begin());

    const MultiFab& vfrac = factory.getVolFrac();
    const auto& areafrac  = factory.getAreaFrac();
    const auto& facecent  = factory.getFaceCent();
    const MultiCutFab& ccent = factory.getCentroid();

    const auto maps = make_maps();
    const GpuArray<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(is_per[0],is_per[1],is_per[2])};

    int nbad = 0;
    Real max_update = Real(0);

    for (MFIter mfi(U_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const auto& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        const auto ftype = flagfab.getType(amrex::grow(bx,4));
        if (ftype == FabType::covered || ftype == FabType::regular) { continue; }

        IArrayBox itracker_fab(amrex::grow(bx,5), itracker_ncomp, The_Async_Arena());
        itracker_fab.setVal<RunOn::Device>(0);

        MakeITracker(bx, AMREX_D_DECL(areafrac[0]->const_array(mfi),
                                      areafrac[1]->const_array(mfi),
                                      areafrac[2]->const_array(mfi)),
                     vfrac.const_array(mfi), itracker_fab.array(), geom, target_volfrac);

        auto const& itracker = itracker_fab.const_array();
        {
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(amrex::grow(bx,4) & domain, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
            {
                return { count_outside(itracker, i, j, k, domain, isp, maps) };
            });
            nbad += amrex::get<0>(reduce_data.value());
        }

        ApplyRedistribution(bx, ncomp, dUdt_out.array(mfi), dUdt_in.array(mfi),
                            U_in.const_array(mfi), scratch.array(mfi), flagfab.const_array(),
                            AMREX_D_DECL(areafrac[0]->const_array(mfi),
                                         areafrac[1]->const_array(mfi),
                                         areafrac[2]->const_array(mfi)),
                            vfrac.const_array(mfi),
                            AMREX_D_DECL(facecent[0]->const_array(mfi),
                                         facecent[1]->const_array(mfi),
                                         facecent[2]->const_array(mfi)),
                            ccent.const_array(mfi), d_bcrec.data(), geom,
                            Real(1), "StateRedist");

        // Covered cells are filled with a large sentinel by design, so only the
        // uncovered cells carry a meaningful update.
        auto const& dudt = dUdt_out.const_array(mfi);
        auto const& flag = flagfab.const_array();
        {
            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
            {
                return { flag(i,j,k).isCovered() ? Real(0) : std::abs(dudt(i,j,k,0)) };
            });
            max_update = amrex::max(max_update, amrex::get<0>(reduce_data.value()));
        }
    }

    ParallelDescriptor::ReduceIntSum(nbad);
    ParallelDescriptor::ReduceRealMax(max_update);

    // A constant state must come back to roundoff, in either precision.
    const Real tol = Real(100.)*std::numeric_limits<Real>::epsilon();
    const bool failed = (nbad > 0) || (max_update > tol);
    amrex::Print() << "  " << label << " : " << nbad << " neighbors outside the domain, "
                   << "max |update| of a constant state = " << max_update
                   << (failed ? "   ** FAILED **" : "") << '\n';
    return failed ? 1 : 0;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // In single precision a plane can cut a cell more than twice; cover those
        // cells instead of aborting.
        ParmParse pp("eb2");
        pp.add("cover_multiple_cuts", 1);

        int nfail = 0;

        amrex::Print() << "Neighborhoods built on prescribed normals at domain edges:\n";
#if (AMREX_SPACEDIM == 2)
        // Corner cell whose fluid opens out of the domain corner.  Each coordinate
        // override has to leave the replacement inside the domain, otherwise the
        // x-override and the y-override undo each other.
        nfail += synthetic_case(8, IntVect(7,7), 0.3_rt, {0.7_rt,0.5_rt},   "hi-hi corner, n = ( .7, .5)");
        nfail += synthetic_case(8, IntVect(0,0), 0.3_rt, {-0.7_rt,-0.5_rt}, "lo-lo corner, n = (-.7,-.5)");
        nfail += synthetic_case(8, IntVect(7,0), 0.3_rt, {0.7_rt,-0.5_rt},  "hi-lo corner, n = ( .7,-.5)");
        nfail += synthetic_case(8, IntVect(0,7), 0.3_rt, {-0.7_rt,0.5_rt},  "lo-hi corner, n = (-.7, .5)");
#else
        nfail += synthetic_case(8, IntVect(7,7,4), 0.3_rt, {0.7_rt,0.5_rt,0.2_rt},    "hi-hi edge,   n = ( .7, .5, .2)");
        nfail += synthetic_case(8, IntVect(0,0,4), 0.3_rt, {-0.7_rt,-0.5_rt,-0.2_rt}, "lo-lo edge,   n = (-.7,-.5,-.2)");
        nfail += synthetic_case(8, IntVect(7,7,7), 0.3_rt, {0.7_rt,0.5_rt,0.2_rt},    "hi corner,    n = ( .7, .5, .2)");
        nfail += synthetic_case(8, IntVect(0,0,0), 0.3_rt, {-0.7_rt,-0.5_rt,-0.2_rt}, "lo corner,    n = (-.7,-.5,-.2)");
        // Equal normal components break symmetry and force the second merge, which
        // has to respect the domain boundary just like the first one.
        nfail += synthetic_case(8, IntVect(4,7,4), 0.3_rt, {0.5_rt,0.5_rt,0.2_rt}, "hi-y face,    n = ( .5, .5, .2)");
        nfail += synthetic_case(8, IntVect(4,7,4), 0.3_rt, {0.2_rt,0.5_rt,0.5_rt}, "hi-y face,    n = ( .2, .5, .5)");
        nfail += synthetic_case(8, IntVect(7,4,4), 0.3_rt, {0.5_rt,0.2_rt,0.5_rt}, "hi-x face,    n = ( .5, .2, .5)");
        nfail += synthetic_case(8, IntVect(4,4,7), 0.3_rt, {0.2_rt,0.5_rt,0.5_rt}, "hi-z face,    n = ( .2, .5, .5)");
#endif

        amrex::Print() << "Redistributing a constant state across a cut domain:\n";
#if (AMREX_SPACEDIM == 2)
        nfail += geometry_case(17, {0,0}, {0.30_rt,0.30_rt}, {0.9_rt,0.8_rt}, 17, "non-periodic,   1 box ");
        nfail += geometry_case(17, {0,0}, {0.42_rt,0.42_rt}, {0.9_rt,0.4_rt},  9, "non-periodic,   4 boxes");
        nfail += geometry_case(17, {1,0}, {0.30_rt,0.30_rt}, {0.9_rt,0.8_rt},  9, "x-periodic,     4 boxes");
        nfail += geometry_case(16, {1,1}, {0.30_rt,0.30_rt}, {0.9_rt,0.8_rt},  4, "fully periodic, 16 boxes");
#else
        nfail += geometry_case(17, {0,0,0}, {0.30_rt,0.30_rt,0.30_rt}, {0.9_rt,0.8_rt,0.7_rt}, 17, "non-periodic,   1 box ");
        nfail += geometry_case(17, {0,0,0}, {0.42_rt,0.42_rt,0.42_rt}, {0.9_rt,0.8_rt,0.7_rt}, 17, "non-periodic,   1 box ");
        nfail += geometry_case(17, {0,0,0}, {0.30_rt,0.30_rt,0.30_rt}, {0.9_rt,-0.8_rt,0.7_rt}, 9, "non-periodic,   8 boxes");
        nfail += geometry_case(17, {0,0,0}, {0.30_rt,0.30_rt,0.30_rt}, {0.9_rt,0.4_rt,0.2_rt}, 17, "non-periodic,   1 box ");
        nfail += geometry_case(17, {1,0,0}, {0.30_rt,0.30_rt,0.30_rt}, {0.9_rt,0.8_rt,0.7_rt},  9, "x-periodic,     8 boxes");
        nfail += geometry_case(16, {1,1,1}, {0.42_rt,0.42_rt,0.42_rt}, {0.9_rt,0.8_rt,0.7_rt},  8, "fully periodic, 8 boxes");
        nfail += geometry_case(16, {1,1,1}, {0.55_rt,0.55_rt,0.55_rt}, {0.6_rt,0.5_rt,0.4_rt},  4, "fully periodic, 64 boxes");
#endif

        if (nfail > 0) {
            amrex::Abort("StateRedist test FAILED: " + std::to_string(nfail) + " failing case(s)");
        }
        amrex::Print() << "StateRedist test PASSED\n";
    }
    amrex::Finalize();
    return 0;
}
