#ifndef _RESTRICTOR_H_
#define _RESTRICTOR_H_

#ifndef _INTERFACE_H_
#  include "interface.H"
#endif
#ifndef _BOUNDARY_H_
#  include "boundary.H"
#endif

typedef const class amr_restrictor_class& amr_restrictor;

// This is the base class for all amr_restrictors.  The basic functions
// box and fill work in the interiors of patches:  Give box the fine
// grid dimensions and it will tell you what coarse patch you can fill
// from it; fill then does the work.  The interface routine is for
// restrictors like node-based full-weighting that do special stuff
// at the edges of grids.

class amr_restrictor_class {
 public:
  amr_restrictor_class() { }
  virtual ~amr_restrictor_class() { }
  virtual Box box(const Box&, const IntVect&) const
    { BoxLib::Error("amr_restrictor_class::box---Restriction function not defined.");
      return Box(); }
  virtual void fill(Fab&, const Box&,
	            Fab&, const IntVect&) const
    { BoxLib::Error("amr_restrictor_class::fill---Restriction function not defined."); }
  virtual void interface(Fab&, const Box&, MultiFab&,
			 const copy_cache*, const level_interface&,
			 const amr_boundary_class&, const IntVect&) const { }
};

class cell_average_restrictor_class : public amr_restrictor_class {
 protected:
  int integrate;
 public:
  cell_average_restrictor_class(int i) { integrate = i; }
  Box box(const Box& fb, const IntVect& rat) const;
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
};

class terrain_velocity_restrictor_class
  : public cell_average_restrictor_class {
 public:
  terrain_velocity_restrictor_class(int i)
    : cell_average_restrictor_class(i) { }
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
};

class injection_restrictor_class : public amr_restrictor_class {
 public:
  injection_restrictor_class() { }
  Box box(const Box& fb, const IntVect& rat) const;
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
};

class default_restrictor_class : public amr_restrictor_class {
 public:
  default_restrictor_class() { }
  Box box(const Box& fb, const IntVect& rat) const;
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
};

// Vanilla bilinear (full-weighting) restrictor class.
// Restricts fine-fine interfaces, even those at the boundary.
// Does nothing at coarse-fine interfaces.
class bilinear_restrictor_class : public amr_restrictor_class {
 protected:
  int integrate;
 public:
  bilinear_restrictor_class(int i) { integrate = i; }
  Box box(const Box& fb, const IntVect& rat) const;
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
  void interface(Fab& patch, const Box& region,
	         MultiFab& fine, const copy_cache*,
		 const level_interface& interface,
		 const amr_boundary_class& bdy, const IntVect& rat) const;
};

// Restricts all interfaces, uses data at coarse points
// along coarse-fine interfaces.
class bilinear_restrictor_coarse_class
 : public bilinear_restrictor_class {
 public:
  bilinear_restrictor_coarse_class(int i)
    : bilinear_restrictor_class(i) { }
  void interface(Fab& patch, const Box& region,
	         MultiFab& fine, const copy_cache*,
		 const level_interface& interface,
		 const amr_boundary_class& bdy, const IntVect& rat) const;
};

const amr_restrictor_class error_restrictor;
const cell_average_restrictor_class cell_average_restrictor(0);
const cell_average_restrictor_class cell_integrator(1);
const injection_restrictor_class injection_restrictor;
const default_restrictor_class default_restrictor;
const bilinear_restrictor_class bilinear_restrictor(0);
const bilinear_restrictor_coarse_class bilinear_restrictor_coarse(0);
const bilinear_restrictor_class bilinear_integrator(1);
const bilinear_restrictor_coarse_class bilinear_integrator_coarse(1);

#endif
