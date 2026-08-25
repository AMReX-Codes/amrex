/**
 * \file main.cpp
 * \brief Test for face-centered EB factory capability (EB2::BuildFC)
 *
 * Validates: (1) EB2::BuildFC() workflow, (2) FC factory creation for x/y/z directions,
 * (3) EB data bounds, (4) sphere volume vs analytical (5% tolerance)
 */

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <cmath>

using namespace amrex;

void main_main()
{
    int nerrors = 0;

    // Parse input parameters
    int n_cell = 64;
    int max_grid_size = 32;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
    }

    // Domain: [0,1]³
    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};

    Geometry geom(Box(IntVect(AMREX_D_DECL(0,0,0)),
                     IntVect(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1))),
                  rb, 0, is_periodic);

    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    // Get sphere parameters
    Real sphere_radius = 0.25;
    Vector<Real> sphere_center(AMREX_SPACEDIM, 0.5);
    int sphere_has_fluid_inside = 1;
    {
        ParmParse pp;
        pp.query("sphere_radius", sphere_radius);
        pp.queryarr("sphere_center", sphere_center);
        pp.query("sphere_has_fluid_inside", sphere_has_fluid_inside);
    }

    // Build CC EB geometry
    EB2::SphereIF sf(sphere_radius,
                     {AMREX_D_DECL(sphere_center[0], sphere_center[1], sphere_center[2])},
                     sphere_has_fluid_inside);
    auto gshop = EB2::makeShop(sf);
    EB2::Build(gshop, geom, 0, 0);

    // Build FC EB data
    EB2::BuildFC();

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);

    // Create CC factory for validation
    auto cc_factory = makeEBFabFactory(&eb_level, ba, dm,
                                       Vector<int>{1,1,1}, EBSupport::full);
    if (cc_factory->faceDir() != -1) {
        ++nerrors;
        amrex::Print() << "ERROR: CC factory faceDir() = " << cc_factory->faceDir() << "\n";
    }

    // Create FC factories and validate
    std::array<std::unique_ptr<EBFArrayBoxFactory>, AMREX_SPACEDIM> fc_factories;

    for (int face_dir = 0; face_dir < AMREX_SPACEDIM; ++face_dir) {
        fc_factories[face_dir] = std::make_unique<EBFArrayBoxFactory>(
            eb_level, geom, ba, dm, Vector<int>{1,1,1}, EBSupport::full, face_dir);

        if (fc_factories[face_dir]->faceDir() != face_dir) {
            ++nerrors;
            amrex::Print() << "ERROR: FC factory dir=" << face_dir << " faceDir()="
                           << fc_factories[face_dir]->faceDir() << "\n";
        }

        const MultiFab& volfrac = fc_factories[face_dir]->getVolFrac();

        // Check volume fraction bounds
        Real vf_min = volfrac.min(0), vf_max = volfrac.max(0);
        if (vf_min < 0.0 || vf_min > 1.0 || vf_max < 0.0 || vf_max > 1.0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << face_dir << " volfrac range ["
                           << vf_min << ", " << vf_max << "]\n";
        }

        // Count cut cells
        const auto& ebflag = fc_factories[face_dir]->getMultiEBCellFlagFab();
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename ReduceData<int>::Type;

        for (MFIter mfi(ebflag); mfi.isValid(); ++mfi) {
            const auto& flag_fab = ebflag[mfi];
            const Box& bx = mfi.validbox();
            auto const& flag_arr = flag_fab.const_array();

            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                    return {flag_arr(i,j,k).isSingleValued() ? 1 : 0};
                });
        }

        int num_cut = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceIntSum(num_cut);

        if (num_cut == 0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << face_dir << " num_cut=0\n";
        }
    }

    // Numerical accuracy validation - check volume for each FC factory
    Real cell_volume = AMREX_D_TERM(geom.CellSize(0),*geom.CellSize(1),*geom.CellSize(2));
    const Real pi = std::numbers::pi_v<Real>;
#if (AMREX_SPACEDIM == 3)
    Real vol_analytical = Real(4.0/3.0) * pi * std::pow(sphere_radius, 3);
#else
    Real vol_analytical = pi * (sphere_radius*sphere_radius);
#endif
    Real vol_tol = Real(1.0e-3);
    {
        ParmParse pp;
        pp.query("vol_tol", vol_tol);
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const MultiFab& volfrac_fc = fc_factories[dir]->getVolFrac();
        Real vol_fluid = volfrac_fc.sum_unique(0) * cell_volume;
        Real vol_error = std::abs(vol_fluid - vol_analytical) / vol_analytical;
        if (vol_error > vol_tol) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << dir << " vol_error=" << vol_error
                           << " exceeds " << vol_tol << "\n";
        }
    }

    // Covered regions. With eb2.max_grid_size small enough that whole boxes fall inside the
    // body, the EB level has covered boxes, and the face-centered data must report them as
    // covered rather than as regular fluid. A face-centered cell is covered exactly when both
    // of the cell-centered cells it straddles are, so classify by the cell-centered flags
    // rather than by the face-centered volume fraction - the failure being guarded against
    // makes that volume fraction 1, which would hide the very cells at issue.
    const auto& ccflag = cc_factory->getMultiEBCellFlagFab();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const MultiFab& volfrac_fc = fc_factories[dir]->getVolFrac();
        const auto& ebflag = fc_factories[dir]->getMultiEBCellFlagFab();
        auto const& afrac = fc_factories[dir]->getAreaFrac();
        const IntVect fdir = IntVect::TheDimensionVector(dir);

        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
        ReduceData<int,int,int,int> reduce_data(reduce_op);
        using ReduceTuple = typename ReduceData<int,int,int,int>::Type;

        for (MFIter mfi(volfrac_fc); mfi.isValid(); ++mfi) {
            auto const& vf    = volfrac_fc.const_array(mfi);
            auto const& flag  = ebflag[mfi].const_array();
            auto const& cflag = ccflag[mfi].const_array();
            // the area fractions carry the cell-centered type, so enclosedCells of the
            // face-centered valid box is exactly their valid range, and those indices are
            // valid in the face-centered arrays too
            const Box abx = amrex::enclosedCells(mfi.validbox());
            const bool have_ap = afrac[0]->ok(mfi);
            AMREX_D_TERM(auto const& apx = have_ap ? afrac[0]->const_array(mfi) : Array4<Real const>{};,
                         auto const& apy = have_ap ? afrac[1]->const_array(mfi) : Array4<Real const>{};,
                         auto const& apz = have_ap ? afrac[2]->const_array(mfi) : Array4<Real const>{};);

            reduce_op.eval(abx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                    const IntVect iv{AMREX_D_DECL(i,j,k)};
                    if (!(cflag(iv-fdir).isCovered() && cflag(iv).isCovered())) {
                        return {0,0,0,0};
                    }
                    int nbad_vf   = (vf(i,j,k) != Real(0.0))     ? 1 : 0;
                    int nbad_flag = (!flag(i,j,k).isCovered())   ? 1 : 0;
                    int nbad_ap   = 0;
                    if (have_ap) {
                        if (AMREX_D_TERM(apx(i,j,k) != Real(0.0),
                                      || apy(i,j,k) != Real(0.0),
                                      || apz(i,j,k) != Real(0.0))) { nbad_ap = 1; }
                    }
                    return {1, nbad_vf, nbad_flag, nbad_ap};
                });
        }

        auto rv = reduce_data.value();
        int ncovered  = amrex::get<0>(rv);
        int nbad_vf   = amrex::get<1>(rv);
        int nbad_flag = amrex::get<2>(rv);
        int nbad_ap   = amrex::get<3>(rv);
        ParallelDescriptor::ReduceIntSum(ncovered);
        ParallelDescriptor::ReduceIntSum(nbad_vf);
        ParallelDescriptor::ReduceIntSum(nbad_flag);
        ParallelDescriptor::ReduceIntSum(nbad_ap);

        if (ncovered == 0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << dir << " no covered cells; set eb2.max_grid_size "
                           << "small enough that whole boxes fall inside the body\n";
        }
        if (nbad_vf > 0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << dir << " " << nbad_vf << " of " << ncovered
                           << " covered cells have a nonzero volume fraction\n";
        }
        if (nbad_flag > 0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << dir << " " << nbad_flag << " of " << ncovered
                           << " covered cells are not flagged covered\n";
        }
        if (nbad_ap > 0) {
            ++nerrors;
            amrex::Print() << "ERROR: dir=" << dir << " " << nbad_ap << " of " << ncovered
                           << " covered cells have a nonzero area fraction\n";
        }
    }

    // Report
    if (nerrors > 0) {
        amrex::Print() << "FCFactory test FAILED with " << nerrors << " errors\n";
        amrex::Abort();
    } else {
        amrex::Print() << "FCFactory test PASSED\n";
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");
        main_main();
    }
    amrex::Finalize();
    return 0;
}
