/*--------------------------------------------------------------------
  associated include
  --------------------------------------------------------------------*/
#include "helpers.H"

/*--------------------------------------------------------------------
  standard includes
  --------------------------------------------------------------------*/
#include <AMReX_BCUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

/*--------------------------------------------------------------------
  defines and static variables
  --------------------------------------------------------------------*/
// Helpers for indexing velocity components
static constexpr int U = 0;
static constexpr int V = 1;
static constexpr int W = 2;

/*--------------------------------------------------------------------
  private free function declarations
  --------------------------------------------------------------------*/
static void SetTemperatureProfile (  //
    amrex::MultiFab& temp,
    const amrex::Geometry& geom );

static void UpdateTemperatureGhosts (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec );

static void FillCoarseFineGhosts (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec,
    int level );

static double Interpolate (  //
    const std::map<double, double>& data,
    double x );

/*--------------------------------------------------------------------
  free function definitions
  --------------------------------------------------------------------*/
amrex::Vector<amrex::Geometry> DefineGeometry ( int nx, int ny, int nz, double dx )
{
    const amrex::RealBox real_box( { 0.0, 0.0, 0.0 },  //
                                   { nx * dx, ny * dx, nz * dx } );
    constexpr amrex::CoordSys::CoordType coord = amrex::CoordSys::CoordType::cartesian;
    const amrex::IntArray is_periodic{ 0, 0, 0 };

    amrex::Vector<amrex::Geometry> geoms;
    for ( int level = 0; level < 3; ++level ) {
        const int multiplier = std::pow( 2, level );
        const amrex::Box domain( amrex::IntVect( 0, 0, 0 ),  //
                                 amrex::IntVect( multiplier * nx - 1,
                                                 multiplier * ny - 1,
                                                 multiplier * nz - 1 ) );
        geoms.emplace_back( domain, real_box, coord, is_periodic );
    }
    return geoms;
}


amrex::Vector<amrex::BoxArray> DefineBoxArrays ( int nx, int ny, int nz )
{
    // Level zero is the full domain
    // Next level is the left half of the full domain, refined 2x
    // Final level is the left quarter of the full domain, refined another 2x

    amrex::Vector<amrex::BoxArray> bas;
    for ( int level = 0; level < 3; ++level ) {
        const int x_multiplier = 1;
        const int yz_multiplier = std::pow( 2, level );
        bas.emplace_back(                               //
            amrex::Box(                                 //
                amrex::IntVect( 0, 0, 0 ),              //
                amrex::IntVect( x_multiplier * nx - 1,  //
                                yz_multiplier * ny - 1,
                                yz_multiplier * nz - 1 ) ) );
    }
    return bas;
}


amrex::Vector<amrex::DistributionMapping> DefineDMs (
    const amrex::Vector<amrex::BoxArray>& bas )
{
    amrex::Vector<amrex::DistributionMapping> dms;
    for ( const auto& ba : bas ) {
        dms.emplace_back( ba );
    }
    return dms;
}


void DefineFABs (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::BoxArray>& bas,
    const amrex::Vector<amrex::DistributionMapping>& dms )
{
    temps.clear();
    vels.clear();
    for ( int level = 0; level < bas.size(); ++level ) {
        constexpr int SingleComp = 1;
        constexpr int TempGhosts = 1;
        constexpr int VelGhosts = 0;
        const auto& ba = bas.at( level );
        const auto Uba = amrex::convert( ba, amrex::IntVect::TheDimensionVector( U ) );
        const auto Vba = amrex::convert( ba, amrex::IntVect::TheDimensionVector( V ) );
        const auto Wba = amrex::convert( ba, amrex::IntVect::TheDimensionVector( W ) );
        const auto& dm = dms.at( level );
        temps.emplace_back( ba, dm, SingleComp, TempGhosts );
        vels.emplace_back( std::array{ amrex::MultiFab( Uba, dm, SingleComp, VelGhosts ),
                                       amrex::MultiFab( Vba, dm, SingleComp, VelGhosts ),
                                       amrex::MultiFab( Wba, dm, SingleComp, VelGhosts ) } );
    }
}


void CalculateVelocities (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels )
{
    // For each level
    for ( auto& vel : vels ) {
        // Initialize all velocities to zero m/s.
        vel.at( U ).setVal( 0.0 );
        vel.at( V ).setVal( 0.0 );
        vel.at( W ).setVal( 0.0 );
    }
}


void CalculateTemperatures (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec )
{
    for ( int level = 0; level < temps.size(); ++level ) {
        auto& temp = temps.at( level );
        const auto& geom = geoms.at( level );

        // Set all non-ghost (valid cell) temperatures based on an analytic expression.
        SetTemperatureProfile( temp, geom );
    }

    // Fill ghosts
    UpdateTemperatureGhosts( temps, geoms, bcrec );

} /* end: CalculateTemperatures() */


static void SetTemperatureProfile (  //
    amrex::MultiFab& temp,
    const amrex::Geometry& geom )
{
    for ( amrex::MFIter mfi( temp ); mfi.isValid(); ++mfi ) {
        const amrex::Box& box = mfi.validbox();
        const auto& temp_arr = temp.array( mfi );

        const auto lo = amrex::lbound( box );
        const auto hi = amrex::ubound( box );
        for ( int i = lo.x; i <= hi.x; ++i ) {
            for ( int j = lo.y; j <= hi.y; ++j ) {
                for ( int k = lo.z; k <= hi.z; ++k ) {
                    // Set "valid" cells with a linear gradient.
                    const amrex::Real x = geom.CellCenter( i, U );
                    temp_arr( i, j, k ) = 300.0 + 100.0 * x;
                }
            }
        }
    }
} /* end: SetTemperatureProfile() */


static void UpdateTemperatureGhosts (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec )
{
    for ( int level = 0; level < temps.size(); ++level ) {
        auto& temp = temps.at( level );
        const auto& geom = geoms.at( level );

        if ( level > 0 ) {
            // Initialize the coarse/fine ghosts from coarse data.
            // At level 0, we don't have a coarser solution to look at, but then we
            // also don't have any coarse/fine ghosts (only internal and external
            // ghosts).
            FillCoarseFineGhosts( temps, geoms, bcrec, level );
        }

        // Fill "internal" ghosts
        // Note: This is not necessary when we only have one FAB per level.
        // TODO: Is this already handled by one of the other ghost fills?
        temp.FillBoundary();

        // Fill "external" ghosts.
        // TODO: Does FillCoarseFineGhosts() already handle this?
        const amrex::Vector<amrex::BCRec> bcrecs( temp.nComp(), bcrec );
        amrex::FillDomainBoundary( temp, geom, bcrecs );
    }
} /* end: UpdateTemperatureGhosts() */


static void FillCoarseFineGhosts (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec,
    int level )
{
    // "Interpolate from coarse" to initialize coarse/fine ghosts on the fine MultiFab

    // MultiFabs
    auto& fine_mf = temps.at( level );
    const auto& coarse_mf = temps.at( level - 1 );

    // Geometry
    const auto& fine_geom = geoms.at( level );
    const auto& coarse_geom = geoms.at( level - 1 );

    // Unused time (passed to PhysBCFunc and otherwise not used?)
    constexpr double time_unused = 0.0;

    // Source, destination, and BC components
    constexpr int scomp = 0;
    constexpr int dcomp = 0;
    constexpr int ncomp = 1;
    constexpr int cbccomp = 0;
    constexpr int fbccomp = 0;
    constexpr int bcscomp = 0;

    // Refinement ratio, per dimension
    const amrex::IntVect ref_ratio( 2 );

    // Domain boundary conditions, per-component
    const amrex::Vector<amrex::BCRec> bcrecs{ bcrec };

    // Select the interpolator:
    amrex::Interpolater* interp = nullptr;

    // Piece-wise constant
    // interp = &amrex::pc_interp; // Does not generate the desired values

    // "Bilinear interpolation on cell centered data."
    interp = &amrex::cell_bilinear_interp;  // Works

    // "Linear conservative interpolation on cell centered data", limited
    // interp = &amrex::lincc_interp;  // Works

    // "Linear conservative interpolation on cell centered data", unlimited
    // interp = &amrex::cell_cons_interp;  // Works

    // Set up BC functors to handle external and coarse/fine ghosts.
    assert( !amrex::Gpu::inLaunchRegion() );  // Use GpuBndryFuncFab for GPU support?
    amrex::CpuBndryFuncFab null_bndry_func = nullptr;
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc( coarse_geom, bcrecs, null_bndry_func );
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc( fine_geom, bcrecs, null_bndry_func );

    if constexpr ( true ) {
        // TODO: I don't believe this function will modify the coarse MultiFab, but
        //       the function signature seems to require a Vector of non-const
        //       pointers.
        auto& non_const_coarse_mf = temps.at( level - 1 );

        // Fill just the coarse/fine ghosts on the fine MultiFab using interpolated
        // values from the coarse MultiFab.
        // TODO: Does this also fill interior or exterior fine ghosts?
        amrex::FillPatchTwoLevels(  //
            fine_mf,                // <- Destination
            time_unused,
            { &non_const_coarse_mf },  // <- Source, coarse at multiple times
            { time_unused },
            { &fine_mf },  // <- Source, fine at multiple times
            { time_unused },
            scomp,
            dcomp,
            ncomp,
            coarse_geom,
            fine_geom,
            cphysbc,
            cbccomp,
            fphysbc,
            fbccomp,
            ref_ratio,
            interp,
            bcrecs,
            bcscomp );

    } else {

        // Fill all fine data using interpolated values from coarse MultiFab.  We
        // only need to fill the coarse/fine ghosts, but this will also overwrite
        // the data in the "valid" cells.
        // TODO: Does this also fill interior or exterior fine ghosts?
        amrex::InterpFromCoarseLevel(  //
            fine_mf,
            time_unused,
            coarse_mf,
            scomp,
            dcomp,
            ncomp,
            coarse_geom,
            fine_geom,
            cphysbc,
            cbccomp,
            fphysbc,
            fbccomp,
            ref_ratio,
            interp,
            bcrecs,
            bcscomp );
    }
}


void CalculateBuoyancy (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    const amrex::Vector<amrex::MultiFab>& temps,
    amrex::Real dt,
    const std::array<amrex::Real, 3>& gs,
    amrex::Real beta,
    amrex::Real T_ref )
{
    // For each level, calculate the buoyancy term and apply it to the velocity
    // field.
    for ( int level = 0; level < vels.size(); ++level ) {
        auto& temp = temps.at( level );
        auto& vel = vels.at( level );

        for ( amrex::MFIter mfi( vel.at( U ) ); mfi.isValid(); ++mfi ) {
            const amrex::Box& bx = mfi.validbox();

            const auto temp_arr = temp.const_array( mfi );
            const auto velU_arr = vel.at( U ).array( mfi );
            const amrex::Real coeff = -1.0 * gs.at( U ) * beta;

            amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept {
                const amrex::Real face_temp = 0.5 * ( temp_arr( i, j, k ) + temp_arr( i - 1, j, k ) );
                const amrex::Real buoyancy_accel = coeff * ( face_temp - T_ref );
                velU_arr( i, j, k ) += dt * buoyancy_accel;
            } );
        }

        for ( amrex::MFIter mfi( vel.at( V ) ); mfi.isValid(); ++mfi ) {
            const amrex::Box& bx = mfi.validbox();

            const auto temp_arr = temp.const_array( mfi );
            const auto velV_arr = vel.at( V ).array( mfi );
            const amrex::Real coeff = -1.0 * gs.at( V ) * beta;

            amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept {
                const amrex::Real face_temp = 0.5 * ( temp_arr( i, j, k ) + temp_arr( i, j - 1, k ) );
                const amrex::Real buoyancy_accel = coeff * ( face_temp - T_ref );
                velV_arr( i, j, k ) += dt * buoyancy_accel;
            } );
        }

        for ( amrex::MFIter mfi( vel.at( W ) ); mfi.isValid(); ++mfi ) {
            const amrex::Box& bx = mfi.validbox();

            const auto temp_arr = temp.const_array( mfi );
            const auto velW_arr = vel.at( W ).array( mfi );
            const amrex::Real coeff = -1.0 * gs.at( W ) * beta;

            amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept {
                const amrex::Real face_temp = 0.5 * ( temp_arr( i, j, k ) + temp_arr( i, j, k - 1 ) );
                const amrex::Real buoyancy_accel = coeff * ( face_temp - T_ref );
                velW_arr( i, j, k ) += dt * buoyancy_accel;
            } );
        }
    }

    // Overwrite every "covered" cell with a value averaged down from the level
    // above it.
    const amrex::IntVect ref_ratio( 2 );
    for ( int coarseLevel = vels.size() - 2; coarseLevel >= 0; --coarseLevel ) {
        auto& coarseVel = vels.at( coarseLevel );
        const auto& fineVel = vels.at( coarseLevel + 1 );
        amrex::average_down_faces( fineVel.at( U ), coarseVel.at( U ), ref_ratio );
        amrex::average_down_faces( fineVel.at( V ), coarseVel.at( V ), ref_ratio );
        amrex::average_down_faces( fineVel.at( W ), coarseVel.at( W ), ref_ratio );
    }
} /* end: CalculateBuoyancy() */


void CheckResults (  //
    const amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    const amrex::Vector<amrex::Geometry>& geoms )
{
    // Analytic result defined by
    //    U(x,y,z) = 100.0 (m/s)/m * x
    //
    // However, the smallest/largest temperature values are at the first/last
    // cell centers at
    //    x_min = DX/2
    // and
    //    x_max = 1.0 m - DX/2.
    //
    // Using a reflect_even boundary condition on the temperatures introduces
    // an error, where the interpolated temperatures on the low/high faces are
    // equal to the temperatures at x_min/x_max.
    //
    // Note that the DX above should correspond to the finest level at the given
    // location, which is different on the high and low ends.
    //
    // This definition of the expected results encodes this error in the
    // expectation.  Without this error, the expectation could be defined at
    // all levels as
    //    { {0.0, 0.0}, {1.0, 100.0} };
    //
    const auto& geom_fine = geoms.back();
    const auto& geom_coarse = geoms.front();

    const double dx_fine = geom_fine.data().CellSize( U );
    const double dx_coarse = geom_coarse.data().CellSize( U );

    const double length = geom_coarse.ProbLength( U );
    assert( geom_fine.ProbLength( U ) == geom_coarse.ProbLength( U ) );

    const double firstCenter = dx_fine / 2.0;
    const double lastCenter = length - dx_coarse / 2.0;

    const std::map<double, double> expectedValues{
        { firstCenter, 100.0 * firstCenter },  // i.e. { 0.03125, 3.125 },
        { lastCenter, 100.0 * lastCenter },    // i.e. { 0.875, 87.5 },
    };

    for ( int level = 0; level < vels.size(); ++level ) {
        const auto& UVel = vels.at( level ).at( U );
        const auto& geom = geoms.at( level );

        for ( amrex::MFIter mfi( UVel ); mfi.isValid(); ++mfi ) {
            const amrex::Box& box = mfi.validbox();
            const auto& arr = UVel.array( mfi );

            const auto lo = amrex::lbound( box );
            const auto hi = amrex::ubound( box );
            for ( int i = lo.x; i <= hi.x; ++i ) {
                for ( int j = lo.y; j <= hi.y; ++j ) {
                    for ( int k = lo.z; k <= hi.z; ++k ) {
                        // Check "valid" cells against the expected expression.
                        const amrex::Real x = geom.LoEdge( i, U );
                        const auto diff = arr( i, j, k ) - Interpolate( expectedValues, x );
                        constexpr auto tol_arbitrary = 1.0E-5;
                        AMREX_ALWAYS_ASSERT( std::abs( diff ) < tol_arbitrary );
                    }
                }
            }
        }
    }
}


static double Interpolate ( const std::map<double, double>& data, double x )
{
    assert( !data.empty() );

    // Clamp if x is out of bounds
    if ( x <= data.begin()->first ) {
        return data.begin()->second;
    }
    if ( x >= data.rbegin()->first ) {
        return data.rbegin()->second;
    }

    // Find the first element with key > x
    const auto upper = data.upper_bound( x );
    assert( upper != data.begin() );
    const auto lower = std::prev( upper );

    const auto [x0, y0] = *lower;
    const auto [x1, y1] = *upper;

    // Linear interpolation
    const double t = ( x - x0 ) / ( x1 - x0 );
    return y0 + t * ( y1 - y0 );
};

/*--------------------------------------------------------------------
  End of file
  --------------------------------------------------------------------*/
