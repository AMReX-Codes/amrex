#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#ifndef _AMR_DEFS_H_
#  include "amr_defs.H"
#endif

typedef const class amr_interpolator_class& amr_interpolator;

// This is the base class for all amr_interpolators.  Currently these
// are self-contained, i.e., they don't need to talk to boundary conditions
// or interface objects or whatever.  Tell box what fine area you want
// filled, and it will tell you how much coarse data it needs.  Just
// fill the requested patch with coarse data and call fill.  How you
// get the coarse data is your problem.

class amr_interpolator_class {
 public:
  amr_interpolator_class() { }
  virtual ~amr_interpolator_class() { }
  virtual Box box(const Box&, int) const
    { BoxLib::Error("amr_interpolator_class::box---Interpolation function not defined.");
      return Box(); }
  virtual void fill(Fab&, const Box&,
		    Fab&, const Box&, int) const
    { BoxLib::Error("amr_interpolator_class::box---Interpolation function not defined."); }
};

class bilinear_interpolator_class : public amr_interpolator_class {
 public:
  bilinear_interpolator_class() { }
  Box box(const Box&, int) const;
  void fill(Fab&, const Box&,
	    Fab&, const Box&, int) const;
};

const amr_interpolator_class error_interpolator;
const bilinear_interpolator_class bilinear_interpolator;

#endif
