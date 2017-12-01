
#include <WarpX.H>
#include <WarpXConst.H>

using namespace amrex;

void
WarpX::MoveWindow (bool move_j)
{
    if (do_moving_window == 0) return;

    // Update the continuous position of the moving window,
    // and of the plasma injection
    moving_window_x += moving_window_v * dt[0];
    int dir = moving_window_dir;
    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::do_plasma_injection and (WarpX::gamma_boost > 1)){
        // In boosted-frame simulations, the plasma has moved since the last
        // call to this function, and injection position needs to be updated
        current_injection_position -= WarpX::beta_boost *
#if ( BL_SPACEDIM == 3 )
            WarpX::boost_direction[dir] * PhysConst::c * dt[0];
#elif ( BL_SPACEDIM == 2 )
            // In 2D, dir=0 corresponds to x and dir=1 corresponds to z
            // This needs to be converted in order to index `boost_direction`
            // which has 3 components, for both 2D and 3D simulations.
            WarpX::boost_direction[2*dir] * PhysConst::c * dt[0];
#endif
    }

    // compute the number of cells to shift on the base level
    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];
    const Real* current_lo = geom[0].ProbLo();
    const Real* current_hi = geom[0].ProbHi();
    const Real* dx = geom[0].CellSize();
    int num_shift_base = static_cast<int>((moving_window_x - current_lo[dir]) / dx[dir]);

    if (num_shift_base == 0) return;

    // update the problem domain. Note the we only do this on the base level because
    // amrex::Geometry objects share the same, static RealBox.
    for (int i=0; i<BL_SPACEDIM; i++) {
        new_lo[i] = current_lo[i];
        new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift_base * dx[dir];
    new_hi[dir] = current_hi[dir] + num_shift_base * dx[dir];
    RealBox new_box(new_lo, new_hi);
    geom[0].ProbDomain(new_box);

    int num_shift      = num_shift_base;
    int num_shift_crse = num_shift;

    // Shift the mesh fields
    for (int lev = 0; lev <= max_level; ++lev) {

        if (lev > 0) {
            num_shift_crse = num_shift;
            num_shift *= refRatio(lev-1)[dir];
        }

        for (int dim = 0; dim < 3; ++dim) {

            shiftMF(*Bfield_fp[lev][dim], geom[lev], num_shift, dir);
            shiftMF(*Efield_fp[lev][dim], geom[lev], num_shift, dir);

            if (move_j) {
                shiftMF(*current_fp[lev][dim], geom[lev], num_shift, dir);
            }

            if (do_dive_cleaning) {
                shiftMF(*F_fp[lev],   geom[lev], num_shift, dir);
                shiftMF(*rho_fp[lev], geom[lev], num_shift, dir);
            }

            if (do_pml && pml[lev]->ok()) {
                const std::array<MultiFab*, 3>& pml_B = pml[lev]->GetB_fp();
                const std::array<MultiFab*, 3>& pml_E = pml[lev]->GetE_fp();
                shiftMF(*pml_B[dim], geom[lev], num_shift, dir);
                shiftMF(*pml_E[dim], geom[lev], num_shift, dir);
                if (do_dive_cleaning) {
                    MultiFab* pml_F = pml[lev]->GetF_fp();
                    shiftMF(*pml_F, geom[lev], num_shift, dir);
                }
            }

            if (lev > 0) {

                shiftMF(*Bfield_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                shiftMF(*Efield_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                shiftMF(*Bfield_aux[lev][dim], geom[lev], num_shift, dir);
                shiftMF(*Efield_aux[lev][dim], geom[lev], num_shift, dir);

                if (move_j) {
                    shiftMF(*current_cp[lev][dim], geom[lev-1], num_shift_crse, dir);
                }

                if (do_dive_cleaning) {
                    shiftMF(*F_cp[lev], geom[lev-1], num_shift_crse, dir);
                    shiftMF(*rho_cp[lev], geom[lev-1], num_shift_crse, dir);
                }

                if (do_pml && pml[lev]->ok()) {
                    const std::array<MultiFab*, 3>& pml_B = pml[lev]->GetB_cp();
                    const std::array<MultiFab*, 3>& pml_E = pml[lev]->GetE_cp();
                    shiftMF(*pml_B[dim], geom[lev-1], num_shift_crse, dir);
                    shiftMF(*pml_E[dim], geom[lev-1], num_shift_crse, dir);
                    if (do_dive_cleaning) {
                        MultiFab* pml_F = pml[lev]->GetF_cp();
                        shiftMF(*pml_F, geom[lev-1], num_shift_crse, dir);
                    }
                }
            }
        }
    }

    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::do_plasma_injection) {
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
        // Perform the injection of new particles in particleBox
        if (particleBox.ok() and (current_injection_position != new_injection_position)){
            for (int i = 0; i < num_injected_species; ++i) {
                int ispecies = injected_plasma_species[i];
                WarpXParticleContainer& pc = mypc->GetParticleContainer(ispecies);
                auto& ppc = dynamic_cast<PhysicalParticleContainer&>(pc);
                ppc.AddPlasma(lev, particleBox);
            }
            // Update the injection position
            current_injection_position = new_injection_position;
        }
    }
}

void
WarpX::shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, int dir)
{
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    const int nc = mf.nComp();
    const int ng = mf.nGrow();

    BL_ASSERT(ng >= num_shift);

    MultiFab tmpmf(ba, dm, nc, ng);
    MultiFab::Copy(tmpmf, mf, 0, 0, nc, ng);
    tmpmf.FillBoundary(geom.periodicity());

    // Make a box that covers the region that the window moved into
    const IndexType& typ = ba.ixType();
    const Box& domainBox = geom.Domain();
    Box adjBox;
    if (num_shift > 0) {
        adjBox = adjCellHi(domainBox, dir, ng);
    } else {
        adjBox = adjCellLo(domainBox, dir, ng);
    }
    adjBox = amrex::convert(adjBox, typ);

    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        if (idim == dir and typ.nodeCentered(dir)) {
            if (num_shift > 0) {
                adjBox.growLo(idim, -1);
            } else {
                adjBox.growHi(idim, -1);
            }
        } else if (idim != dir) {
            adjBox.growLo(idim, ng);
            adjBox.growHi(idim, ng);
        }
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi )
    {
        FArrayBox& srcfab = tmpmf[mfi];

        Box outbox = mfi.fabbox();
        outbox &= adjBox;
        if (outbox.ok()) {  // outbox is the region that the window moved into
            srcfab.setVal(0.0, outbox, 0, nc);
        }

        FArrayBox& dstfab = mf[mfi];
        dstfab.setVal(0.0);

        Box dstBox = dstfab.box();

        if (num_shift > 0) {
            dstBox.growHi(dir, -num_shift);
        } else {
            dstBox.growLo(dir,  num_shift);
        }

        dstfab.copy(srcfab, amrex::shift(dstBox,dir,num_shift), 0, dstBox, 0, nc);
    }
}
