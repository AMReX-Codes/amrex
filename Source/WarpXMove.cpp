
#include <WarpX.H>

using namespace amrex;

void
WarpX::MoveWindow (bool move_j)
{
    if (do_moving_window == 0) return;

    if (max_level > 0) {
        amrex::Abort("MoveWindow with mesh refinement has not been implemented");
    }

    const int lev = 0;

    // compute the number of cells to shift
    int dir = moving_window_dir;
    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];
    const Real* current_lo = geom[lev].ProbLo();
    const Real* current_hi = geom[lev].ProbHi();
    const Real* dx = geom[lev].CellSize();
    moving_window_x += moving_window_v * dt[lev];
    int num_shift = static_cast<int>((moving_window_x - current_lo[dir]) / dx[dir]);

    if (num_shift == 0) return;

    // update the problem domain
    for (int i=0; i<BL_SPACEDIM; i++) {
        new_lo[i] = current_lo[i];
        new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift * dx[dir];
    new_hi[dir] = current_hi[dir] + num_shift * dx[dir];
    RealBox new_box(new_lo, new_hi);
    geom[lev].ProbDomain(new_box);

    // shift the mesh fields (Note - only on level 0 for now)
    shiftMF(*Bfield_fp[lev][0], geom[lev], num_shift, dir);
    shiftMF(*Bfield_fp[lev][1], geom[lev], num_shift, dir);
    shiftMF(*Bfield_fp[lev][2], geom[lev], num_shift, dir);
    shiftMF(*Efield_fp[lev][0], geom[lev], num_shift, dir);
    shiftMF(*Efield_fp[lev][1], geom[lev], num_shift, dir);
    shiftMF(*Efield_fp[lev][2], geom[lev], num_shift, dir);
    if (move_j) {
        shiftMF(*current_fp[lev][0], geom[lev], num_shift, dir);
        shiftMF(*current_fp[lev][1], geom[lev], num_shift, dir);
        shiftMF(*current_fp[lev][2], geom[lev], num_shift, dir);
    }

    InjectPlasma(num_shift, dir);

    // Redistribute (note - this removes particles that are outside of the box)
    mypc->Redistribute();
}

void
WarpX::shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, int dir)
{
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    const int nc = mf.nComp();
    const int ng = std::max(mf.nGrow(), std::abs(num_shift));
    MultiFab tmpmf(ba, dm, nc, ng);
    MultiFab::Copy(tmpmf, mf, 0, 0, nc, ng);
    tmpmf.FillBoundary(geom.periodicity());

    // Zero out the region that the window moved into
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
