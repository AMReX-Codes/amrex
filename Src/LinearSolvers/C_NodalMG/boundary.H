#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "interface.H"

#ifdef HG_USE_CACHE
class copy_cache;
#endif

typedef const class amr_boundary_class& amr_boundary;

// This is the base class for all amr_boundaries.  The default function dir
// should always be acceptable.  The default box and check_against_boundary
// functions are appropriate for wall-type boundaries.  The default
// boundary_mesh function is appropriate for all boundary types I've
// created so far, since it is little more than a front end for the
// internal work routine check_against_boundary.

// The domain (2nd Box) arguments to dir and box should be cell-based
// regardless of the type of the region.

class amr_boundary_class {
 protected:
  virtual void check_against_boundary(BoxList& bl, List<int>& il,
				      const Box& b, int ib,
				      const Box& d, int dim1) const;
 public:
  amr_boundary_class() { }
  virtual ~amr_boundary_class() { }
  int dir(const Box&, const Box&) const;
  virtual int defined() const { return 0; }
  virtual Box box(const Box&, const Box&, int) const;
  virtual void fill(FArrayBox&, const Box&,
		    const FArrayBox&, const Box&, const Box&, int) const
    { BoxLib::Error("amr_boundary_class::fill---boundary fill function not defined."); }
  virtual void fill(FArrayBox&, const Box&,
		    const MultiFab&, int, const Box&) const
    { BoxLib::Error("amr_boundary_class::fill---boundary fill function not defined."); }
  virtual void sync_borders(MultiFab&,
			    const level_interface&) const { }
  virtual void fill_borders(MultiFab& r,
			    const level_interface& interface,
			    int w = -1) const;
#ifdef HG_USE_CACHE
  virtual void set_sync_cache(copy_cache*, int, int&,
			      MultiFab&,
			      const level_interface&) const { }
  virtual void set_border_cache(copy_cache*, int, int&,
				MultiFab&,
				const level_interface&, int) const {
    BoxLib::Error("amr_boundary_class::set_border_cache---boundary condition not defined.");
  }
#endif
  virtual void boundary_mesh(BoxArray& exterior_mesh, int *&grid_ref,
			     const BoxArray& interior_mesh,
			     const Box& domain) const;
  virtual void duplicate(List<Box>&, const Box&) const { }
  virtual int fill_patch_special(FArrayBox&, const Box&,
				 MultiFab&, const Box&, int) const {
    return 0;
  }
  virtual int singular() const {
    BoxLib::Error("amr_boundary_class::singular---boundary condition not defined.");
    return 0;
  }
};

inline int amr_boundary_class::dir(const Box& region, const Box& domain) const
{
#if 1
  for ( int i = 0; i < BL_SPACEDIM; ++i )
  {
      if ( region.bigEnd(i) < domain.smallEnd(i) )
          return -(i+1);
      if ( region.smallEnd(i) > domain.bigEnd(i) )
          return +(i+1);
  }
  for ( int i = 0; i < BL_SPACEDIM; ++i )
  {
      if ( region.bigEnd(i) == domain.smallEnd(i) )
	  return -(i+1);
      if ( region.smallEnd(i) == domain.bigEnd(i) )
	  return (i+1);
  }
  BoxLib::Error("amr_boundary_class::dir---boundary box not outside domain.");
  return 0;

#else
  int retval = 0;
  if      (region.bigEnd(0) < domain.smallEnd(0))
    retval = -1;
  else if (region.smallEnd(0) > domain.bigEnd(0))
    retval = 1;
  else if (region.bigEnd(1) < domain.smallEnd(1))
    retval = -2;
  else if (region.smallEnd(1) > domain.bigEnd(1))
    retval = 2;
#if (BL_SPACEDIM == 3)
  else if (region.bigEnd(2) < domain.smallEnd(2))
    retval = -3;
  else if (region.smallEnd(2) > domain.bigEnd(2))
    retval = 3;
#endif
  else if (region.bigEnd(0) == domain.smallEnd(0))
    retval = -1;
  else if (region.smallEnd(0) == domain.bigEnd(0))
    retval = 1;
  else if (region.bigEnd(1) == domain.smallEnd(1))
    retval = -2;
  else if (region.smallEnd(1) == domain.bigEnd(1))
    retval = 2;
#if (BL_SPACEDIM == 3)
  else if (region.bigEnd(2) == domain.smallEnd(2))
    retval = -3;
  else if (region.smallEnd(2) == domain.bigEnd(2))
    retval = 3;
#endif
  else
    BoxLib::Error("amr_boundary_class::dir---boundary box not outside domain.");
  return retval;
#endif
}

class reflection_boundary_class : public amr_boundary_class {
 public:
  reflection_boundary_class() { }
  int defined() const { return 1; }
  void fill(FArrayBox&, const Box&,
	    const FArrayBox&, const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const MultiFab&, int, const Box&) const
    { BoxLib::Error("reflection_boundary_class::fill---interface version not defined"); }
  int singular() const {
    return 1;
  }
};

class negation_boundary_class : public amr_boundary_class {
 public:
  negation_boundary_class() { }
  int defined() const { return 1; }
  void fill(FArrayBox&, const Box&,
	    const FArrayBox&, const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const MultiFab&, int, const Box&) const
    { BoxLib::Error("negation_boundary_class::fill---interface version not defined"); }
  int singular() const {
    return 0;
  }
};

class periodic_boundary_class : public amr_boundary_class {
 protected:
  void check_against_boundary(BoxList& bl, List<int>& il,
			      const Box& b, int ib,
			      const Box& d, int dim1) const;
 public:
  periodic_boundary_class() { }
  int defined() const { return 1; }
  Box box(const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const FArrayBox&, const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const MultiFab&, int, const Box&) const
    { BoxLib::Error("periodic_boundary_class::fill---interface version not defined"); }
  void duplicate(List<Box>& bl, const Box& domain) const;
  int singular() const {
    return 1;
  }
};

const amr_boundary_class error_boundary;
const reflection_boundary_class reflection_boundary;
const negation_boundary_class negation_boundary;
const periodic_boundary_class periodic_boundary;

#include "RegType.H"

// This is the base class for physical boundary conditions of the
// types used in the AMR code.  Members return appropriate
// amr_boundary objects for various flow variables.  The defaults
// are all error_boundary, so this class should only be used as
// a base for derivation.

class amr_fluid_boundary_class {
 protected:
  const amr_boundary_class *v[BL_SPACEDIM];
  const amr_boundary_class *s;
  const amr_boundary_class *p;
  const amr_boundary_class *ts;
 public:
  amr_fluid_boundary_class();
  virtual ~amr_fluid_boundary_class() { }

  virtual amr_boundary velocity(int i) const {
    return *v[i];
  }
  virtual amr_boundary scalar() const {
    return *s;
  }
  virtual amr_boundary pressure() const {
    return *p;
  }
  virtual amr_boundary terrain_sigma() const {
    return *ts;
  }
};

typedef const amr_fluid_boundary_class& amr_fluid_boundary;

// This is the object returned by the inviscid_fluid_boundary class.
// It chooses among periodic, reflection and negation boundary behavior
// based on the boundary types accessed through the pointer and the
// flow direction flowdim.
// flowdim = -1 indicates an advected scalar
// flowdim = -2 indicates a pressure
// flowdim = -3 indicates a streamfunction
// flowdim = -4 indicates the terrain sigma (all components)

class inviscid_fluid_boundary;

class mixed_boundary_class : public amr_boundary_class {
 protected:
  inviscid_fluid_boundary *ptr;
  int flowdim;
  void check_against_boundary(BoxList& bl, List<int>& il,
			      const Box& b, int ib,
			      const Box& d, int dim1) const;
 public:
  mixed_boundary_class(inviscid_fluid_boundary *Ptr, int idim) {
    ptr = Ptr; flowdim = idim;
  }

  int defined() const { return 1; }
  Box box(const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const FArrayBox&, const Box&, const Box&, int) const;
  void fill(FArrayBox&, const Box&,
	    const MultiFab&, int, const Box&) const;
  void sync_borders(MultiFab&,
		    const level_interface&) const;
  void fill_borders(MultiFab& r,
		    const level_interface& interface,
		    int w = -1) const;
#ifdef HG_USE_CACHE
  void set_sync_cache(copy_cache* cache, int nsets, int& iset,
		      MultiFab& r,
		      const level_interface& interface) const;
  void set_border_cache(copy_cache* cache, int nsets, int& iset,
			MultiFab& r,
			const level_interface& interface,
			int w) const;
#endif
  void duplicate(List<Box>& bl, const Box& domain) const;
  int fill_patch_special(FArrayBox&, const Box&,
			 MultiFab&, const Box&, int flags) const;
  int singular() const;
};

// This is a simple boundary condition with slip walls.
// Takes an array of RegType specifiers, one for each side of
// the domain.  Currently only periodic and refWall type
// boundaries are supported.

class inviscid_fluid_boundary : public amr_fluid_boundary_class {
  friend mixed_boundary_class;
 protected:
  RegType bc[BL_SPACEDIM][2];
 public:
  inviscid_fluid_boundary(RegType Bc[BL_SPACEDIM][2]);
  ~inviscid_fluid_boundary();
};

#endif
