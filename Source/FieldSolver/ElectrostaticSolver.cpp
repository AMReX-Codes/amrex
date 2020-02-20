/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeTensorLaplacian.H>

#include <WarpX.H>

using namespace amrex;

void
WarpX::ComputeSpaceChargeField (bool const reset_fields)
{
    if (reset_fields) {
        // Reset all E and B fields to 0, before calculating space-charge fields
        const int num_levels = max_level + 1;
        for (int lev = 0; lev <= max_level; lev++) {
            for (int comp=0; comp<3; comp++) {
                Efield_fp[lev][comp]->setVal(0);
                Bfield_fp[lev][comp]->setVal(0);
            }
        }
    }

    // Loop over the species and add their space-charge contribution to E and B
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        if (species.initialize_self_fields || do_electrostatic) {
            AddSpaceChargeField(species);
        }
    }
}

void
WarpX::AddSpaceChargeField (WarpXParticleContainer& pc)
{

#ifdef WARPX_DIM_RZ
    amrex::Abort("The initialization of space-charge field has not yet been implemented in RZ geometry.");
#endif

    // Allocate fields for charge and potential
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    const int ng = WarpX::nox;
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rho[lev].reset(new MultiFab(nba, dmap[lev], 1, ng)); // Make ng big enough/use rho from sim
        phi[lev].reset(new MultiFab(nba, dmap[lev], 1, 1));
        phi[lev]->setVal(0.);
    }

    // Deposit particle charge density (source of Poisson solver)
    bool const local = false;
    bool const reset = true;
    bool const do_rz_volume_scaling = true;
    pc.DepositCharge(rho, local, reset, do_rz_volume_scaling);

    // Get the particle beta vector
    bool const local_average = false; // Average across all MPI ranks
    std::array<Real, 3> beta = pc.meanParticleVelocity(local_average);
    for (Real& beta_comp : beta) beta_comp /= PhysConst::c; // Normalize

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi, beta, pc.self_fields_required_precision );

    // Compute the corresponding electric and magnetic field, from the potential phi
    computeE( Efield_fp, phi, beta );
    computeB( Bfield_fp, phi, beta );

}

/* Compute the potential `phi` by solving the Poisson equation with `rho` as
   a source, assuming that the source moves at a constant speed \f$\vec{\beta}\f$.
   This uses the amrex solver.

   More specifically, this solves the equation
   \f[
       \vec{\nabla}^2\phi - (\vec{\beta}\cdot\vec{\nabla})^2\phi = -\frac{\rho}{\epsilon_0}
   \f]

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<Real, 3> const beta,
                   Real const required_precision) const
{
    // Define the boundary conditions
    Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    for (int idim=0; idim<AMREX_SPACEDIM; idim++){
        if ( Geom(0).isPeriodic(idim) ) {
            lobc[idim] = LinOpBCType::Periodic;
            hibc[idim] = LinOpBCType::Periodic;
        } else {
            // Use Dirichlet boundary condition by default.
            // Ideally, we would often want open boundary conditions here.
            lobc[idim] = LinOpBCType::Dirichlet;
            hibc[idim] = LinOpBCType::Dirichlet;
        }
    }

    // Define the linear operator (Poisson operator)
    MLNodeTensorLaplacian linop( Geom(), boxArray(), DistributionMap() );
    linop.setDomainBC( lobc, hibc );
    // Set the value of beta
    amrex::Array<amrex::Real,AMREX_SPACEDIM> beta_solver =
#if (AMREX_SPACEDIM==2)
        {{ beta[0], beta[2] }};  // beta_x and beta_z
#else
        {{ beta[0], beta[1], beta[2] }};
#endif
    linop.setBeta( beta_solver );

    // Solve the Poisson equation
    MLMG mlmg(linop);
    mlmg.setVerbose(2);
    mlmg.solve( GetVecOfPtrs(phi), GetVecOfConstPtrs(rho), required_precision, 0.0);

    // Normalize by the correct physical constant
    for (int lev=0; lev < rho.size(); lev++){
        phi[lev]->mult(-1./PhysConst::ep0);
    }
}

/* \bried Compute the electric field that corresponds to `phi`, and
          add it to the set of MultiFab `E`.

   The electric field is calculated by assuming that the source that
   produces the `phi` potential is moving with a constant speed \f$\vec{\beta}\f$:
   \f[
    \vec{E} = -\vec{\nabla}\phi + (\vec{\beta}\cdot\vec{\beta})\phi \vec{\beta}
   \f]
   (where the second term represent the term \f$\partial_t \vec{A}\f$, in
    the case of a moving source)

   \param[inout] E Electric field on the grid
   \param[in] phi The potential from which to compute the electric field
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computeE (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                 const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                 std::array<amrex::Real, 3> const beta ) const
{
    for (int lev = 0; lev <= max_level; lev++) {

        const Real* dx = Geom(lev).CellSize();

#ifdef _OPENMP
        #pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1./dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1./dx[1];
            const Real inv_dz = 1./dx[2];
#else
            const Real inv_dz = 1./dx[1];
#endif
            const Box& tbx  = mfi.tilebox(Ex_nodal_flag);
            const Box& tby  = mfi.tilebox(Ey_nodal_flag);
            const Box& tbz  = mfi.tilebox(Ez_nodal_flag);

            const auto& phi_arr = phi[lev]->array(mfi);
            const auto& Ex_arr = (*E[lev][0])[mfi].array();
            const auto& Ey_arr = (*E[lev][1])[mfi].array();
            const auto& Ez_arr = (*E[lev][2])[mfi].array();

            Real beta_x = beta[0];
            Real beta_y = beta[1];
            Real beta_z = beta[2];

            // Calculate the electric field
            // Use discretized derivative that matches the staggering of the grid.
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) +=
                        +(beta_x*beta_x-1)*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) )
                        +beta_x*beta_y*0.25*inv_dy*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k))
                        +beta_x*beta_z*0.25*inv_dz*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k-1)
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ey_arr(i,j,k) +=
                        +beta_y*beta_x*0.25*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_y-1)*inv_dy*( phi_arr(i,j+1,k)-phi_arr(i,j,k) )
                        +beta_y*beta_z*0.25*inv_dz*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k-1)
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25*inv_dx*(phi_arr(i+1,j,k  )-phi_arr(i-1,j,k  )
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i-1,j,k+1))
                        +beta_z*beta_y*0.25*inv_dy*(phi_arr(i,j+1,k  )-phi_arr(i,j-1,k  )
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j-1,k+1))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j,k+1)-phi_arr(i,j,k) );
                }
            );
#else
            amrex::ParallelFor( tbx, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) +=
                        +(beta_x*beta_x-1)*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) )
                        +beta_x*beta_z*0.25*inv_dz*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j+1,k)-phi_arr(i,j,k) );
                }
            );
#endif
        }
    }
}


/* \bried Compute the magnetic field that corresponds to `phi`, and
          add it to the set of MultiFab `B`.

   The magnetic field is calculated by assuming that the source that
   produces the `phi` potential is moving with a constant speed \f$\vec{\beta}\f$:
   \f[
    \vec{B} = -\frac{1}{c}\vec{\beta}\times\vec{\nabla}\phi
   \f]
   (this represents the term \f$\vec{\nabla} \times \vec{A}\f$, in the case of a moving source)

   \param[inout] E Electric field on the grid
   \param[in] phi The potential from which to compute the electric field
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computeB (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& B,
                 const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                 std::array<amrex::Real, 3> const beta ) const
{
    for (int lev = 0; lev <= max_level; lev++) {

        const Real* dx = Geom(lev).CellSize();

#ifdef _OPENMP
        #pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1./dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1./dx[1];
            const Real inv_dz = 1./dx[2];
#else
            const Real inv_dz = 1./dx[1];
#endif
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

            const auto& phi_arr = phi[0]->array(mfi);
            const auto& Bx_arr = (*B[lev][0])[mfi].array();
            const auto& By_arr = (*B[lev][1])[mfi].array();
            const auto& Bz_arr = (*B[lev][2])[mfi].array();

            Real beta_x = beta[0];
            Real beta_y = beta[1];
            Real beta_z = beta[2];

            constexpr Real inv_c = 1./PhysConst::c;

            // Calculate the magnetic field
            // Use discretized derivative that matches the staggering of the grid.
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bx_arr(i,j,k) += inv_c * (
                        -beta_y*inv_dz*0.5*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k)
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k))
                        +beta_z*inv_dy*0.5*(phi_arr(i,j+1,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j,k+1)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    By_arr(i,j,k) += inv_c * (
                        -beta_z*inv_dx*0.5*(phi_arr(i+1,j,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i+1,j,k+1)-phi_arr(i,j,k+1))
                        +beta_x*inv_dz*0.5*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bz_arr(i,j,k) += inv_c * (
                        -beta_x*inv_dy*0.5*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i+1,j,k))
                        +beta_y*inv_dx*0.5*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i,j+1,k)));
                }
            );
#else
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bx_arr(i,j,k) += inv_c * (
                        -beta_y*inv_dz*( phi_arr(i,j+1,k)-phi_arr(i,j,k) ));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    By_arr(i,j,k) += inv_c * (
                        -beta_z*inv_dx*0.5*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i,j+1,k))
                        +beta_x*inv_dz*0.5*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i+1,j,k)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bz_arr(i,j,k) += inv_c * (
                        +beta_y*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) ));
                }
            );
#endif
        }
    }
}
