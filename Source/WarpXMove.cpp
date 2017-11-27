
#include <WarpX.H>

using namespace amrex;

void
WarpX::MoveWindow (bool move_j)
{
    if (do_moving_window == 0) return;

    // compute the number of cells to shift on the base level
    int dir = moving_window_dir;
    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];
    const Real* current_lo = geom[0].ProbLo();
    const Real* current_hi = geom[0].ProbHi();
    const Real* dx = geom[0].CellSize();
    moving_window_x += moving_window_v * dt[0];
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
    for (int lev = 0; lev <= max_level; ++lev) {

        if (lev > 0) {
            num_shift_crse = num_shift;
            num_shift *= refRatio(lev-1)[dir];
        }

        // shift the mesh fields
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

    InjectPlasma(num_shift_base, dir);
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
