/*
 *  Testbed implementing a buoyancy calculation on a MAC grid on a multi-level
 *  mesh without subcycling-in-time.
 */

/*--------------------------------------------------------------------
  standard includes
  --------------------------------------------------------------------*/
#include <AMReX.H>

/*--------------------------------------------------------------------
  non-standard includes
  --------------------------------------------------------------------*/
#include "helpers.H"

/*--------------------------------------------------------------------
  forward declarations
  --------------------------------------------------------------------*/
int MyMain ();

/*--------------------------------------------------------------------
  function definitions
  --------------------------------------------------------------------*/
int main (int argc, char** argv) {
    amrex::Initialize(argc, argv);
    const int ret = MyMain();
    amrex::Finalize();
    return ret;
}

int MyMain ()
{
    //
    // Allocation of data structures
    //

    constexpr int nx = 4;
    constexpr int ny = 1;
    constexpr int nz = 1;
    constexpr double dx = 0.25;  // meters

    constexpr amrex::Real dt = 1.0;                             // Seconds
    constexpr std::array<amrex::Real, 3> gs{ -1.0, 0.0, 0.0 };  // m/s
    constexpr amrex::Real beta = 1.0;                           // 1/Kelvin
    constexpr amrex::Real T_ref = 300.0;                        // Kelvin

    constexpr auto re = amrex::BCType::reflect_even;
    const amrex::BCRec bcrec( re, re, re, re, re, re );

    const amrex::Vector<amrex::Geometry> geoms = DefineGeometry( nx, ny, nz, dx );
    const amrex::Vector<amrex::BoxArray> bas = DefineBoxArrays( nx, ny, nz );
    const amrex::Vector<amrex::DistributionMapping> dms = DefineDMs( bas );

    amrex::Vector<std::array<amrex::MultiFab, 3>> vels;
    amrex::Vector<amrex::MultiFab> temps;
    DefineFABs( vels, temps, bas, dms );


    //
    // Solution
    //

    const int nTimeSteps = 1;
    for ( int timeStep = 0; timeStep < nTimeSteps; ++timeStep ) {

        // Initial velocities and temperatures are assumed to be calculated by some
        // other part of the overall simulation
        CalculateVelocities( vels );
        CalculateTemperatures( temps, geoms, bcrec );

        // At every level, calculate the buoyancy term and modify the velocities
        CalculateBuoyancy( vels, temps, dt, gs, beta, T_ref );
    }

    //
    // Post-process
    //

    CheckResults( vels, geoms );

    return 0;
}

/*--------------------------------------------------------------------
  End of file
  --------------------------------------------------------------------*/
