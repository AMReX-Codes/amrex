#ifndef _FILL_PATCH_H_
#define _FILL_PATCH_H_

#include "interface.H"
#include "boundary.H"
#include "cache.H"
#include "interpolator.H"
#include "restrictor.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FFCPY     fcpy_
#else
#  define FFCPY     FCPY
#endif

extern "C" {
  void FFCPY(Real*, intS, intS, const Real*, intS, const int&);
}

/*
class alias_fab {
 protected:
  FArrayBox* ptr;
  int owned;
 public:
  
  ~alias_fab() {
    if (owned) delete ptr;
  }
  FArrayBox& operator FArrayBox() {
    return *ptr;
  }

 private:
  
};
*/

inline void internal_copy(MultiFab& r,
			  int destgrid, int srcgrid, const Box& b) {
  Real *const dptr = r[destgrid].dataPtr();
  Real *const sptr = r[srcgrid].dataPtr();
  const Box& dbx = r[destgrid].box();
  const Box& sbx = r[srcgrid].box();
  FFCPY(dptr, dimlist(dbx), dimlist(b), sptr, dimlist(sbx), r.nComp());
}

Real inner_product(const MultiFab& r, const MultiFab& s);

/*
const level_real&
    initialize(Real (*f)(const Intvect&,const Intvect&,int,int)) const;
*/

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

/*
  grid_real get_patch(const Box& region,
		      const level_interface& interface,
		      amr_boundary bdy = error_boundary,
		      int flags = 0);
  int get_patch(grid_real& patch, const Box& region,
		const level_interface& interface,
		amr_boundary bdy = error_boundary,
		int flags = 0);
*/

// returns grid number if patch is contained in a single grid,
// uses only 2-bit of flags

int find_patch(const Box& region, const MultiFab& r, int flags = 0);

int fill_patch(FArrayBox& patch, const Box& region,
	       const MultiFab& r,
	       const level_interface& interface,
	       amr_boundary bdy = error_boundary,
	       int flags = 0, int idim = -1, int index = -1);

inline int fill_patch(FArrayBox& patch,
		      const MultiFab& r,
		      const level_interface& interface,
		      amr_boundary bdy = error_boundary,
		      int flags = 0, int idim = -1, int index = -1)
{
  return fill_patch(patch, patch.box(), r, interface, bdy,
		    flags, idim, index);
}

/*
  // fill_patch returns 1 if successful, 0 if not
  // if (flags & 1) is 0, it may return 1 even if patch is not filled
  int fill_patch(FArrayBox& patch, const Box& region,
		 const Box& active,
		 amr_boundary bdy = error_boundary,
		 int flags = 0, int igrid = -1);
  int fill_patch(FArrayBox& patch, const Box& region,
		 amr_boundary bdy = error_boundary,
		 int flags = 0, int igrid = -1)
    { return fill_patch(patch, region, region, bdy, flags, igrid); }
  int fill_patch(FArrayBox& patch,
		 amr_boundary bdy = error_boundary,
		 int flags = 0, int igrid = -1)
    { return fill_patch(patch, patch.box(), bdy, flags, igrid); }
*/

  // fills region of patch wherever it intersects interior grids
  // does not check to see if all of patch is filled, ignores boundaries
int fill_patch_blindly(FArrayBox& patch, const Box& region,
		       MultiFab& r,
		       int flags = 0);

  // fills region of patch wherever it intersects interior grids
  // does not check to see if all of patch is filled, ignores boundaries
int fill_exterior_patch_blindly(FArrayBox& patch, const Box& region,
				MultiFab& r,
				const level_interface& interface,
				amr_boundary bdy,
				int flags = 0);

  // Fills all ghost cells that can be obtained from other grids or
  // from the boundary conditions.  If w is between 0 and border(),
  // inclusive, it will be used as the width of the regions to be
  // filled---otherwise, the full width will be used.

void sync_internal_borders(MultiFab& r, const level_interface& interface);
void fill_internal_borders(MultiFab& r, const level_interface& interface,
			   int w = -1);

inline void sync_borders(MultiFab& r,
#ifdef HG_USE_CACHE
			 const copy_cache* sync_cache,
#endif
			 const level_interface& interface,
			 amr_boundary bdy = error_boundary)
{
#ifdef HG_USE_CACHE
  if (sync_cache) {
    // assumes cache built properly---does not check current bdy
    sync_cache->run();
    return;
  }
#endif
    sync_internal_borders(r, interface);
    bdy.sync_borders(r, interface);
}

inline void fill_borders(MultiFab& r,
#ifdef HG_USE_CACHE
			 const copy_cache* border_cache,
#endif
			 const level_interface& interface,
			 amr_boundary bdy = error_boundary,
			 int w = -1)
{
#ifdef HG_USE_CACHE
  if (border_cache) {
    // assumes cache built properly---does not check current bdy and w
    border_cache->run();
    return;
  }
#endif
    fill_internal_borders(r, interface, w);
    bdy.fill_borders(r, interface, w);
}

void clear_part_interface(MultiFab& r, const level_interface& interface);

void interpolate_patch(FArrayBox& patch, const Box& region,
		       MultiFab& r, const IntVect& rat,
		       const amr_interpolator& interp = bilinear_interpolator(),
		       const level_interface& interface = null_level_interface,
		       amr_boundary bdy = error_boundary);
inline
void interpolate_patch(FArrayBox& patch,
		       MultiFab& r, const IntVect& rat,
		       const amr_interpolator& interp = bilinear_interpolator(),
		       const level_interface& interface = null_level_interface,
		       amr_boundary bdy = error_boundary)
{
  interpolate_patch(patch, patch.box(), r, rat, interp, interface, bdy);
}

void restrict_patch(FArrayBox& patch, const Box& region,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache = 0,
#endif
		    amr_restrictor restric = default_restrictor_class(),
		    const level_interface& interface = null_level_interface,
		    amr_boundary bdy = error_boundary);
inline
void restrict_patch(FArrayBox& patch,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache = 0,
#endif
		    amr_restrictor restric = default_restrictor_class(),
		    const level_interface& interface = null_level_interface,
		    amr_boundary bdy = error_boundary)
{
  restrict_patch(patch, patch.box(), r, rat, 
#ifdef HG_USE_CACHE
      border_cache,
#endif
		 restric, interface, bdy);
}

void restrict_level(MultiFab& dest, int bflag,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache = 0,
#endif
		    amr_restrictor restric = default_restrictor_class(),
		    const level_interface& interface = null_level_interface,
		    amr_boundary bdy = error_boundary);
inline
void restrict_level(MultiFab& dest,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache = 0,
#endif
		    amr_restrictor restric = default_restrictor_class(),
		    const level_interface& interface = null_level_interface,
		    amr_boundary bdy = error_boundary)
{
  restrict_level(dest, 0, r, rat, 
#ifdef HG_USE_CACHE
      border_cache, 
#endif
      restric, interface, bdy);
}

#endif
