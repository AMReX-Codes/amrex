/*
 *  Testbed implementing a buoyancy calculation on a MAC grid on a multi-level
 *  mesh without subcycling-in-time.
 *
 *
 *  Physics
 *  =======
 *
 *  Uses the Boussinesq approximation to calculate a change in velocity due to
 *  temperature-induced density differences:
 *
 *    \Delta \vec{V} = -\beta ( T - T_{ref} ) \vec{g} \Delta t
 *
 *  or
 *
 *    U_i = U_i - dt * beta * g_x * ( T_i - T_ref )
 *
 *  Temperatures are stored at cell centers and velocity components are stored
 *  at cell faces, so the temperature at the cell face (T_i above) has to be
 *  interpolated from the neighboring cell-center temperatures.
 *
 *
 *  Geometry and Mesh
 *  =================
 *
 *               ════════════════════════════════╗
 *              ╱       ╱       ╱       ╱       ╱║
 *             ╔═══════╤═══════╤═══════╤═══════╗ ║      \vec{g} ◄────
 *  Level 0:   ║       │       │       │       ║ ║
 *             ║       │       │       │       ║ ║
 *             ║       │       │       │       ║╱
 *             ╚═══════╧═══════╧═══════╧═══════╝
 *
 *               ════════════════╗
 *              ╱   ╱   ╱   ╱   ╱║
 *             ╔═══╤═══╤═══╤═══╗ ║
 *  Level 1:   ║   │   │   │   ║╱║
 *             ╟───┼───┼───┼───╢ ║
 *             ║   │   │   │   ║╱
 *             ╚═══╧═══╧═══╧═══╝
 *
 *               ════════╗
 *              ╱ ╱ ╱ ╱ ╱║
 *             ╔═╤═╤═╤═╗╱║
 *  Level 2:   ╟─┼─┼─┼─╢╱║
 *             ╟─┼─┼─┼─╢╱║
 *             ╟─┼─┼─┼─╢╱
 *             ╚═╧═╧═╧═╝
 *
 *               ════════════════════════════════╗
 *              ╱ ╱ ╱ ╱ ╱   ╱   ╱       ╱       ╱║
 *             ╔═╤═╤═╤═╦═══╤═══╦═══════╤═══════╗ ║
 *  Composite: ╟─┼─┼─┼─╢   │   ║       │       ║ ║
 *             ╟─┼─┼─┼─╫───┼───╢       │       ║ ║
 *             ╟─┼─┼─┼─╢   │   ║       │       ║╱
 *             ╚═╧═╧═╧═╩═══╧═══╩═══════╧═══════╝
 *
 *  Manually allocates the multi-level mesh to minimize dependencies and allow
 *  us to create a very specific multi-level mesh.  The coarse domain is 4x1x1,
 *  and each finer level only spans the "left half" of the "l-1" domain below
 *  it.
 *
 *    L2:  ├───┼───┼───┼───┤   ·   ·   ·   ·   ·   ·   ·   ·   ·   ·   ·   :
 *
 *    L1:  ├───────┼───────┼───────┼───────┤       ·       ·       ·       :
 *
 *    L0:  ├───────────────┼───────────────┼───────────────┼───────────────┤
 *        0.0 m                           0.5 m                           1.0 m
 *
 *  The velocity data is stored at cell faces.  For a given level, we use an
 *  array of three MultiFabs to separately store the U/V/W velocity components.
 *  We have no use for ghosts with the velocity data.
 *
 *  The temperature data is stored at cell centers.  We require one layer of
 *  ghosts to be able to calculate temperatures at cell faces for the velocity
 *  update calculation.
 *
 *  Each MultiFab contains a single FAB.  This mesh has "domain" boundaries
 *  and "coarse/fine" boundaries, but no "interior" boundaries.
 *
 *
 *  Initialization and Boundary Conditions
 *  ======================================
 *
 *  The initial test case specifies temperatures which vary linearly, assigned
 *  by the formula
 *
 *      T(x,y,z) = 300 K + 100 K/m * x.
 *
 *  The data dependencies for calculating the "delta U" term at 0.5 m on
 *  the coarse grid looks like this.  We need a temperature at the U face,
 *  but temperatures are stored at cell centers.  So we interpolate, by
 *  averaging the two neighboring cell-centered temperatures, T_a and T_b.
 *
 *                                         |
 *                                        \|/
 *                                         '
 *                                         U
 *    L0:  ├───────────────┼────── x ──────┼────── x ──────┼───────────────┤
 *        0.0 m                    |      0.5 m    |                      1.0 m
 *                                T_a             T_b
 *
 *  Since the test case temperature profile is linear, this linear
 *  interpolation is exact, generating the same 350 K temperature as the
 *  analytic equation.
 *
 *  On the medium grid, the temperature on the right is outside the mesh
 *  at this level, making it a "coarse/fine ghost."  If we initialize the
 *  coarse/fine ghosts via a linear interpolation from the coarse mesh,
 *  the ghost will have the exact expected value and our interpolated face
 *  temperature will be an exact 350 K, matching the coarse grid.
 *
 *                                         |
 *                                        \|/
 *                                         '     .─── Ghost Temperature
 *                                         U    /
 *    L1:  ├───────┼───────┼───────┼── x ──┤   x                           :
 *        0.0 m                        |       |                          1.0 m
 *                                    T_c     T_d
 *
 *  We also have ghosts at the domain boundaries.  Exploiting knowledge of
 *  our analytic temperature distribution, we could set these ghosts to
 *  the exact expected values via extrapolation from the temperatures
 *  inside the domain.  However, this is probably not a good general solution.
 *  Instead we use an "even reflection" domain boundary condition.  The
 *  temperature immediately inside the domain boundary is copied to the ghost
 *  immediately outside the domain boundary.
 *
 *
 *        Ghost T -.       .- Boundary T
 *                  \     /
 *    L2:            x ├ x ┼───┼───┼───┼   ·   ·   ·   ·   ·
 *
 *    L1:          x   ├── x ──┼───────┼───────┼───────┼   ·
 *
 *    L0:      x       ├────── x ──────┼───────────────┼──────
 *            /       0.0 m     \                     0.5 m
 *  Ghost T -'                   '- Boundary T
 *
 *  The calculated face temperature now corresponds to the temperature half a
 *  cell width inside the domain.  On the x=0 domain boundary, this
 *  overestimates the temperature by 12.5 K on the coarse grid, but by only
 *  3.125 K on the fine grid.
 *
 *  The following graph illustrates the difference between the effective
 *  temperature distribution represented by this choice of boundary conditions,
 *  and the analytic temperature distribution we are attempting to represent.
 *
 *        Temp (K)
 *             ┌───┼──────────┼───┐     ┌────────────────────────────┐
 *             │   :           XX │     │ *  Effective temperatures  │
 *        400  ┼   :         XX   │     │ X  Analytic temperatures   │
 *             │   :       *******│     └────────────────────────────┘
 *             │   :     **       │
 *             │   :   **         │
 *             │*******           │
 *        300  ┼   XX             │
 *             │ XX:              │
 *             └───┼──────────┼───┘ X (m)
 *                0.0        1.0
 *
 *  Zooming in on the left edge shows how this varies at the different levels.
 *
 *        Temp (K)
 *             ┌───┼─┼─┼───┼──────┐     ┌──────────────────────────────┐
 *        ·    │   :           ** │     │ ## Effective temperature L0  │
 *        ·    │   :         **   │     │ == Effective temperature L1  │
 *     312.500 ┼###########**     │     │ -- Effective temperature L2  │
 *     309.375 │   :     **       │     │ XX Analytic temperature      │
 *     306.250 ┼=======**         │     └──────────────────────────────┘
 *     303.125 ┼-----**           │
 *     300.000 ┼···XX·············│
 *        ·    │ XX               │
 *             └───┼─┼─┼───┼──────┘ X (m)
 *                0.0     0.125
 *
 *
 *  Syncronization of Levels
 *  ========================
 *
 *  Note: This describes what is believed to be the desired implementation,
 *  which varies from what is current implemented and may vary from best
 *  practices.
 *
 *  On the Way Up: Syncronizing Temperatures
 *  ========================================
 *
 *  Before calculating the velocity change, it is assumed the temperature for
 *  all valid cells has been updated for the current timestep.  We do that with
 *  `SetTemperatureProfile()`.
 *
 *  Next, all ghosts need to be updated for the given level.
 *
 *  First, external ghosts are updated first via `FillDomainBoundary()`.  This
 *  imposes the "reflect_even" condition.
 *
 *  Next, coarse/fine ghosts are interpolated from the coarse data via
 *  `FillPatchTwoLevels()`.
 *
 *  Internal ghosts are then filled from other FABs at the same level via
 *  `FillBoundary()`.
 *
 *  Now the velocity change can be calculated for the given level, via
 *  `CalculateBuoyancy()`.
 *
 *  On the Way Down: Syncronizing Velocities
 *  ========================================
 *
 *  After updating velocities at the finest level, we need to bring our better
 *  fine level results onto our coarser levels.
 *
 *  TODO: Document how we average down fine results on the coarse level.
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

        // Debug
        WriteAllResults( vels, temps, geoms, "init" );

        // At every level, calculate the buoyancy term and modify the velocities
        CalculateBuoyancy( vels, temps, dt, gs, beta, T_ref );

        // Debug
        WriteAllResults( vels, temps, geoms, "final" );
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
