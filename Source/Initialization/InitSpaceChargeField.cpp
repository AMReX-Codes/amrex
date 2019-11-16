
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <WarpX.H>

using namespace amrex;

void
WarpX::InitSpaceChargeField (WarpXParticleContainer& pc)
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
        phi[lev].reset(new MultiFab(nba, dmap[lev], 1, 0));
        phi[lev]->setVal(0.);
    }

    // Deposit particle charge density (source of Poisson solver)
    bool const local = false;
    bool const reset = true;
    bool const do_rz_volume_scaling = true;
    pc.DepositCharge(rho, local, reset, do_rz_volume_scaling);

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi );

    // Compute the corresponding electric field, from the potential phi
    computeE( Efield_fp, phi );

}

/* Compute the potential `phi` by solving the Poisson equation with `rho` as
   a source. This uses the amrex solver.

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
*/
void
WarpX::computePhi(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const
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
    MLNodeLaplacian linop( Geom(), boxArray(), DistributionMap() );
    linop.setDomainBC( lobc, hibc );
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray cba = boxArray(lev);
        cba.enclosedCells();
        MultiFab sigma(cba, DistributionMap(lev), 1, 0);
        sigma.setVal(-PhysConst::ep0);
        linop.setSigma(lev, sigma);
    }

    // Solve the Poisson equation
    MLMG mlmg(linop);
    mlmg.setVerbose(1);
    const Real reltol = 1.e-11;
    mlmg.solve( GetVecOfPtrs(phi), GetVecOfConstPtrs(rho), reltol, 0.0);

}

/* \bried Compute the electric field that corresponds to `phi`, and
          add it to the set of MultiFab `E`.

   \param[inout] E Electric field on the grid
   \param[in] phi The potential from which to compute the electric field
*/
void
WarpX::computeE(amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
          const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const
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

#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) += -inv_dx*( phi_arr(i+1,j,k) - phi_arr(i,j,k) );
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ey_arr(i,j,k) += -inv_dy*( phi_arr(i,j+1,k) - phi_arr(i,j,k) );
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) += -inv_dz*( phi_arr(i,j,k+1) - phi_arr(i,j,k) );
                }
            );
#else
            amrex::ParallelFor( tbx, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) += -inv_dx*( phi_arr(i+1,j,k) - phi_arr(i,j,k) );
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) += -inv_dz*( phi_arr(i,j+1,k) - phi_arr(i,j,k) );
                }
            );
#endif
        }
    }
}
