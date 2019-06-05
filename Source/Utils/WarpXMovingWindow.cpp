
#include <WarpX.H>
#include <WarpXConst.H>

using namespace amrex;

void
WarpX::UpdatePlasmaInjectionPosition (Real a_dt)
{
    int dir = moving_window_dir;
    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::warpx_do_continuous_injection and (WarpX::gamma_boost > 1)){
        // In boosted-frame simulations, the plasma has moved since the last
        // call to this function, and injection position needs to be updated
        current_injection_position -= WarpX::beta_boost *
#if ( AMREX_SPACEDIM == 3 )
            WarpX::boost_direction[dir] * PhysConst::c * a_dt;
#elif ( AMREX_SPACEDIM == 2 )
            // In 2D, dir=0 corresponds to x and dir=1 corresponds to z
            // This needs to be converted in order to index `boost_direction`
            // which has 3 components, for both 2D and 3D simulations.
            WarpX::boost_direction[2*dir] * PhysConst::c * a_dt;
#endif
    }
}

int
WarpX::MoveWindow (bool move_j)
{
    if (do_moving_window == 0) return 0;

    // Update the continuous position of the moving window,
    // and of the plasma injection
    moving_window_x += moving_window_v * dt[0];
    int dir = moving_window_dir;

    // Update warpx.current_injection_position
    // PhysicalParticleContainer uses this injection position
    UpdatePlasmaInjectionPosition( dt[0] );
    if (WarpX::warpx_do_continuous_injection){
        // Update injection position for WarpXParticleContainer in mypc.
        // Nothing to do for PhysicalParticleContainers
        // For LaserParticleContainer, need to update the antenna position.
        mypc->UpdateContinuousInjectionPosition( dt[0] );
    }

    // compute the number of cells to shift on the base level
    Real new_lo[AMREX_SPACEDIM];
    Real new_hi[AMREX_SPACEDIM];
    const Real* current_lo = geom[0].ProbLo();
    const Real* current_hi = geom[0].ProbHi();
    const Real* cdx = geom[0].CellSize();
    int num_shift_base = static_cast<int>((moving_window_x - current_lo[dir]) / cdx[dir]);

    if (num_shift_base == 0) return 0;

    // update the problem domain. Note the we only do this on the base level because
    // amrex::Geometry objects share the same, static RealBox.
    for (int i=0; i<AMREX_SPACEDIM; i++) {
        new_lo[i] = current_lo[i];
        new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift_base * cdx[dir];
    new_hi[dir] = current_hi[dir] + num_shift_base * cdx[dir];

    ResetProbDomain(RealBox(new_lo, new_hi));

    // Moving slice coordinates - lo and hi - with moving window //
    // slice box is modified only if slice diagnostics is initialized in input //
    if ( slice_plot_int > 0 )
    {
       Real new_slice_lo[AMREX_SPACEDIM];
       Real new_slice_hi[AMREX_SPACEDIM];
       const Real* current_slice_lo = slice_realbox.lo();
       const Real* current_slice_hi = slice_realbox.hi();
       for ( int i = 0; i < AMREX_SPACEDIM; i++) {
           new_slice_lo[i] = current_slice_lo[i];
           new_slice_hi[i] = current_slice_hi[i];
       }
       int num_shift_base_slice = static_cast<int> ((moving_window_x - 
                                  current_slice_lo[dir]) / cdx[dir]);
       new_slice_lo[dir] = current_slice_lo[dir] + num_shift_base_slice*cdx[dir];
       new_slice_hi[dir] = current_slice_hi[dir] + num_shift_base_slice*cdx[dir];
       slice_realbox.setLo(new_slice_lo);
       slice_realbox.setHi(new_slice_hi);
    }

    int num_shift      = num_shift_base;
    int num_shift_crse = num_shift;

    // Shift the mesh fields
    for (int lev = 0; lev <= finest_level; ++lev) {

        if (lev > 0) {
            num_shift_crse = num_shift;
            num_shift *= refRatio(lev-1)[dir];
        }

        // Shift each component of vector fields (E, B, j)
        for (int dim = 0; dim < 3; ++dim) {

            // Fine grid
            shiftMF(*Bfield_fp[lev][dim], geom[lev], num_shift, dir);
            shiftMF(*Efield_fp[lev][dim], geom[lev], num_shift, dir);
            if (move_j) {
                shiftMF(*current_fp[lev][dim], geom[lev], num_shift, dir);
            }
            if (do_pml && pml[lev]->ok()) {
                const std::array<MultiFab*, 3>& pml_B = pml[lev]->GetB_fp();
                const std::array<MultiFab*, 3>& pml_E = pml[lev]->GetE_fp();
                shiftMF(*pml_B[dim], geom[lev], num_shift, dir);
                shiftMF(*pml_E[dim], geom[lev], num_shift, dir);
            }

            if (lev > 0) {
                // Coarse grid
                shiftMF(*Bfield_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                shiftMF(*Efield_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                shiftMF(*Bfield_aux[lev][dim], geom[lev], num_shift, dir);
                shiftMF(*Efield_aux[lev][dim], geom[lev], num_shift, dir);
                if (move_j) {
                    shiftMF(*current_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                }
                if (do_pml && pml[lev]->ok()) {
                    const std::array<MultiFab*, 3>& pml_B = pml[lev]->GetB_cp();
                    const std::array<MultiFab*, 3>& pml_E = pml[lev]->GetE_cp();
                    shiftMF(*pml_B[dim], geom[lev-1], num_shift_crse, dir);
                    shiftMF(*pml_E[dim], geom[lev-1], num_shift_crse, dir);
                }
            }
        }

        // Shift scalar component F for dive cleaning
        if (do_dive_cleaning) {
            // Fine grid
            shiftMF(*F_fp[lev],   geom[lev], num_shift, dir);
            if (do_pml && pml[lev]->ok()) {
                MultiFab* pml_F = pml[lev]->GetF_fp();
                shiftMF(*pml_F, geom[lev], num_shift, dir);
            }
            if (lev > 0) {
                // Coarse grid
                shiftMF(*F_cp[lev], geom[lev-1], num_shift_crse, dir);
                if (do_pml && pml[lev]->ok()) {
                    MultiFab* pml_F = pml[lev]->GetF_cp();
                    shiftMF(*pml_F, geom[lev-1], num_shift_crse, dir);
                }
                shiftMF(*rho_cp[lev], geom[lev-1], num_shift_crse, dir);
            }
        }

        // Shift scalar component rho
        if (move_j) {
            if (rho_fp[lev]){
                // Fine grid
                shiftMF(*rho_fp[lev],   geom[lev], num_shift, dir);
                if (lev > 0){
                    // Coarse grid
                    shiftMF(*rho_cp[lev], geom[lev-1], num_shift_crse, dir);
                }
            }
        }
    }

    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::warpx_do_continuous_injection) {

        const int lev = 0;

        // particleBox encloses the cells where we generate particles
        // (only injects particles in an integer number of cells,
        // for correct particle spacing)
        RealBox particleBox = geom[lev].ProbDomain();
        Real new_injection_position;
        if (moving_window_v >= 0){
            // Forward-moving window
            Real dx = geom[lev].CellSize(dir);
            new_injection_position = current_injection_position +
                std::floor( (geom[lev].ProbHi(dir) - current_injection_position)/dx ) * dx;
        } else {
            // Backward-moving window
            Real dx = geom[lev].CellSize(dir);
            new_injection_position = current_injection_position -
                std::floor( (current_injection_position - geom[lev].ProbLo(dir))/dx) * dx;
        }
        // Modify the corresponding bounds of the particleBox
        if (moving_window_v >= 0) {
            particleBox.setLo( dir, current_injection_position );
            particleBox.setHi( dir, new_injection_position );
        } else {
            particleBox.setLo( dir, new_injection_position );
            particleBox.setHi( dir, current_injection_position );
        }

        if (particleBox.ok() and (current_injection_position != new_injection_position)){
            // Performs continuous injection of all WarpXParticleContainer
            // in mypc.
            mypc->ContinuousInjection(particleBox);
            current_injection_position = new_injection_position;
        }
    }

    return num_shift_base;
}

void
WarpX::shiftMF (MultiFab& mf, const Geometry& geom, int num_shift, int dir)
{
    BL_PROFILE("WarpX::shiftMF()");
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    const int nc = mf.nComp();
    const IntVect& ng = mf.nGrowVect();

    AMREX_ALWAYS_ASSERT(ng.min() >= num_shift);

    MultiFab tmpmf(ba, dm, nc, ng, MFInfo().SetDeviceFab(false));
    MultiFab::Copy(tmpmf, mf, 0, 0, nc, ng);
    tmpmf.FillBoundary(geom.periodicity());

    // Make a box that covers the region that the window moved into
    const IndexType& typ = ba.ixType();
    const Box& domainBox = geom.Domain();
    Box adjBox;
    if (num_shift > 0) {
        adjBox = adjCellHi(domainBox, dir, ng[dir]);
    } else {
        adjBox = adjCellLo(domainBox, dir, ng[dir]);
    }
    adjBox = amrex::convert(adjBox, typ);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim == dir and typ.nodeCentered(dir)) {
            if (num_shift > 0) {
                adjBox.growLo(idim, -1);
            } else {
                adjBox.growHi(idim, -1);
            }
        } else if (idim != dir) {
            adjBox.growLo(idim, ng[idim]);
            adjBox.growHi(idim, ng[idim]);
        }
    }

    IntVect shiftiv(0);
    shiftiv[dir] = num_shift;
    Dim3 shift = shiftiv.dim3();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi )
    {
        auto const& dstfab = mf.array(mfi);
        auto const& srcfab = tmpmf.array(mfi);

        const Box& outbox = mfi.fabbox() & adjBox;
        if (outbox.ok()) {
            AMREX_PARALLEL_FOR_4D ( outbox, nc, i, j, k, n,
            {
                srcfab(i,j,k,n) = 0.0;
            });
        }

        Box dstBox = mf[mfi].box();
        if (num_shift > 0) {
            dstBox.growHi(dir, -num_shift);
        } else {
            dstBox.growLo(dir,  num_shift);
        }
        AMREX_PARALLEL_FOR_4D ( dstBox, nc, i, j, k, n,
        {
            dstfab(i,j,k,n) = srcfab(i+shift.x,j+shift.y,k+shift.z,n);
        });
    }
}

void
WarpX::ResetProbDomain (const RealBox& rb)
{
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        Geometry g = Geom(lev);
        g.ProbDomain(rb);
        SetGeometry(lev, g);
    }
}
