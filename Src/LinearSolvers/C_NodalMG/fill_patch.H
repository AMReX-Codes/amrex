#ifndef _FILL_PATCH_H_
#define _FILL_PATCH_H_

#include "interface.H"
#include "boundary.H"
#include "cache.H"
#include "interpolator.H"
#include "restrictor.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FFCPY     fcpy_
#else
#  define FORT_FFCPY     FCPY
#endif

extern "C" 
{
    void FORT_FFCPY(Real*, intS, intS, const Real*, intS, const int&);
}

inline void internal_copy(MultiFab& r, int destgrid, int srcgrid, const Box& b) 
{
    Real *const dptr = r[destgrid].dataPtr();
    Real *const sptr = r[srcgrid].dataPtr();
    const Box& dbx = r[destgrid].box();
    const Box& sbx = r[srcgrid].box();
    FORT_FFCPY(dptr, DIMLIST(dbx), DIMLIST(b), sptr, DIMLIST(sbx), r.nComp());
}

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

int find_patch(const Box& region, const MultiFab& r, int flags = 0);

bool fill_patch(FArrayBox& patch, const Box& region,
		const MultiFab& r,
		const level_interface& lev_interface,
		const amr_boundary_class& bdy = error_boundary,
		int flags = 0, int idim = -1, int index = -1);

inline bool fill_patch(FArrayBox& patch,
		       const MultiFab& r,
		       const level_interface& lev_interface,
		       const amr_boundary_class& bdy = error_boundary,
		       int flags = 0, int idim = -1, int index = -1)
{
    return fill_patch(patch, patch.box(), r, lev_interface, bdy, flags, idim, index);
}

// Fills all ghost cells that can be obtained from other grids or
// from the boundary conditions.  If w is between 0 and border(),
// inclusive, it will be used as the width of the regions to be
// filled---otherwise, the full width will be used.

void sync_borders(MultiFab& r,
		  const copy_cache* sync_cache,
		  const level_interface& lev_interface, 
		  const amr_boundary_class& bdy = error_boundary);

void fill_borders(MultiFab& r,
		  const copy_cache* border_cache,
		  const level_interface& lev_interface, const amr_boundary_class& bdy = error_boundary, int w = -1);

void clear_part_interface(MultiFab& r, const level_interface& lev_interface);

void interpolate_patch(FArrayBox& patch, const Box& region,
		       const MultiFab& r, const IntVect& rat,
		       const amr_interpolator_class& interp = bilinear_interpolator_class(),
		       const level_interface& lev_interface = null_level_interface,
		       const amr_boundary_class& bdy = error_boundary);
inline void interpolate_patch(FArrayBox& patch,
			      const MultiFab& r, const IntVect& rat,
			      const amr_interpolator_class& interp = bilinear_interpolator_class(),
			      const level_interface& lev_interface = null_level_interface,
			      const amr_boundary_class& bdy = error_boundary)
{
    interpolate_patch(patch, patch.box(), r, rat, interp, lev_interface, bdy);
}

void restrict_patch(FArrayBox& patch, const Box& region,
		    MultiFab& r, const IntVect& rat,
		    const copy_cache* border_cache = 0,
		    const amr_restrictor_class& restric = default_restrictor(),
		    const level_interface& lev_interface = null_level_interface,
		    const amr_boundary_class& bdy = error_boundary);

inline void  restrict_patch(FArrayBox& patch,
			    MultiFab& r, const IntVect& rat,
			    const copy_cache* border_cache = 0,
			    const amr_restrictor_class& restric = default_restrictor(),
			    const level_interface& lev_interface = null_level_interface,
			    const amr_boundary_class& bdy = error_boundary)
{
    restrict_patch(patch, patch.box(), r, rat, border_cache, restric, lev_interface, bdy);
}

void restrict_level(MultiFab& dest, int bflag, MultiFab& r, const IntVect& rat,
		    const copy_cache* border_cache = 0,
		    const amr_restrictor_class& restric = default_restrictor(),
		    const level_interface& lev_interface = null_level_interface,
		    const amr_boundary_class& bdy = error_boundary);

inline
void restrict_level(MultiFab& dest, MultiFab& r, const IntVect& rat,
		    const copy_cache* border_cache = 0,
		    const amr_restrictor_class& restric = default_restrictor(),
		    const level_interface& lev_interface = null_level_interface,
		    const amr_boundary_class& bdy = error_boundary)
{
    restrict_level(dest, 0, r, rat,  border_cache, restric, lev_interface, bdy);
}

#endif
