#ifndef _FILL_PATCH_H_
#define _FILL_PATCH_H_

#include "interface.H"
#include "boundary.H"
#include "interpolator.H"
#include "restrictor.H"

//void internal_copy(MultiFab& r, int destgrid, int srcgrid, const Box& b);

Real inner_product(const MultiFab& r, const MultiFab& s);

// returns grid number if patch is contained in a single grid, -1 otherwise.

// int find_patch(const Box& region, const MultiFab& r);

// Fills all ghost cells that can be obtained from other grids or
// from the boundary conditions.  If w is between 0 and border(),
// inclusive, it will be used as the width of the regions to be
// filled---otherwise, the full width will be used.

void sync_borders(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy);

void fill_borders(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy, int w, bool hg_terrain);

void clear_part_interface(MultiFab& r, const level_interface& lev_interface);

void interpolate_patch(MultiFab& dmf, int dgrid, const Box& region,
		       const MultiFab& r, const IntVect& rat,
		       const amr_interpolator_class& interp,
		       const level_interface& lev_interface);

void restrict_level(MultiFab& dest, 
		    MultiFab& r, const IntVect& rat,
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class* bdy);

#endif
