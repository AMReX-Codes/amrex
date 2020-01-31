/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <SpectralKSpace.H>
#include <SpectralSolver.H>
#include <PsatdAlgorithm.H>
#include <PMLPsatdAlgorithm.H>

/* \brief Initialize the spectral Maxwell solver
 *
 * This function selects the spectral algorithm to be used, allocates the
 * corresponding coefficients for the discretized field update equation,
 * and prepares the structures that store the fields in spectral space.
 *
 * \param norder_x Order of accuracy of the spatial derivatives along x
 * \param norder_y Order of accuracy of the spatial derivatives along y
 * \param norder_z Order of accuracy of the spatial derivatives along z
 * \param nodal    Whether the solver is applied to a nodal or staggered grid
 * \param dx       Cell size along each dimension
 * \param dt       Time step
 * \param pml      Whether the boxes in which the solver is applied are PML boxes
 */
SpectralSolver::SpectralSolver(
                const amrex::BoxArray& realspace_ba,
                const amrex::DistributionMapping& dm,
                const int norder_x, const int norder_y,
                const int norder_z, const bool nodal,
                const amrex::RealVect dx, const amrex::Real dt,
                const bool pml ) {

    // Initialize all structures using the same distribution mapping dm

    // - Initialize k space object (Contains info about the size of
    // the spectral space corresponding to each box in `realspace_ba`,
    // as well as the value of the corresponding k coordinates)
    const SpectralKSpace k_space= SpectralKSpace(realspace_ba, dm, dx);

    // - Select the algorithm depending on the input parameters
    //   Initialize the corresponding coefficients over k space
    if (pml) {
        algorithm = std::unique_ptr<PMLPsatdAlgorithm>( new PMLPsatdAlgorithm(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt ) );
    } else {
        algorithm = std::unique_ptr<PsatdAlgorithm>( new PsatdAlgorithm(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt ) );
    }

    // - Initialize arrays for fields in spectral space + FFT plans
    field_data = SpectralFieldData( realspace_ba, k_space, dm,
            algorithm->getRequiredNumberOfFields() );

};
