#ifndef _FILL_PATCH_H_
#define _FILL_PATCH_H_

#include "interface.H"
#include "boundary.H"
#include "cache.H"
#include "interpolator.H"
#include "restrictor.H"

void internal_copy(MultiFab& r, int destgrid, int srcgrid, const Box& b);

Real inner_product(const MultiFab& r, const MultiFab& s);

// individual bits of flags are interpreted separately:
// 1-bit:  action when part of patch cannot be filled
//   0 -> ignore (i.e., fill-on-intersect) (default)
//   1 -> return 0 indicating error
// 2-bit:  controls use of ghost cells of source level_real
//   0 -> ignore values in ghost cells (default)
//   1 -> use values from ghost cells
// 4-bit:  what to do with part of patch that can't be filled
//   0 -> ignore (default)
//   1 -> fill this section with zeroes to avoid floating-point errors

// returns grid number if patch is contained in a single grid,
// uses only 2-bit of flags

int find_patch(const Box& region, const MultiFab& r, int flags);

bool fill_patch(FArrayBox& patch, const Box& region,
		const MultiFab& r,
		const level_interface& lev_interface,
		const amr_boundary_class* bdy,
		int flags, int idim = -1, int index = -1);



// Fills all ghost cells that can be obtained from other grids or
// from the boundary conditions.  If w is between 0 and border(),
// inclusive, it will be used as the width of the regions to be
// filled---otherwise, the full width will be used.

void sync_borders(MultiFab& r,
		  const copy_cache* sync_cache,
		  const level_interface& lev_interface, 
		  const amr_boundary_class* bdy);

void fill_borders(MultiFab& r,
		  const copy_cache* border_cache,
		  const level_interface& lev_interface, const amr_boundary_class* bdy, int w);

void clear_part_interface(MultiFab& r, const level_interface& lev_interface);

void interpolate_patch(FArrayBox& patch, const Box& region,
		       const MultiFab& r, const IntVect& rat,
		       const amr_interpolator_class& interp,
		       const level_interface& lev_interface,
		       const amr_boundary_class* bdy);

void restrict_level(MultiFab& dest, bool bflag, MultiFab& r, const IntVect& rat,
		    const copy_cache* border_cache,
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class* bdy);

#endif
