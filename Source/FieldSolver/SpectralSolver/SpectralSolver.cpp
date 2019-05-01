#include <SpectralKSpace.H>
#include <SpectralSolver.H>
#include <PsatdAlgorithm.H>

/* \brief Initialize the spectral Maxwell solver
 *
 * This function selects the spectral algorithm to be used, allocates the
 * corresponding coefficients for the discretized field update equation,
 * and prepares the structures that store the fields in spectral space.
 */
SpectralSolver::SpectralSolver(
                const amrex::BoxArray& realspace_ba,
                const amrex::DistributionMapping& dm,
                const int norder_x, const int norder_y,
                const int norder_z, const bool nodal,
                const amrex::RealVect dx, const amrex::Real dt ) {

    // Initialize all structures using the same distribution mapping dm

    // - Initialize k space object (Contains info about the size of
    // the spectral space corresponding to each box in `realspace_ba`,
    // as well as the value of the corresponding k coordinates)
    const SpectralKSpace k_space= SpectralKSpace(realspace_ba, dm, dx);

    // - Select the algorithm depending on the input parameters
    //   Initialize the corresponding coefficients over k space
    // TODO: Add more algorithms + selection depending on input parameters
    //       For the moment, this only uses the standard PsatdAlgorithm
    algorithm = std::unique_ptr<PsatdAlgorithm>( new PsatdAlgorithm(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt ) );

    // - Initialize arrays for fields in spectral space + FFT plans
    field_data = SpectralFieldData( realspace_ba, k_space, dm );

};
