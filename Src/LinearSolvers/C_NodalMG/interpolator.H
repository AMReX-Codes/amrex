#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "amr_defs.H"

// This is the base class for all amr_interpolators.  Currently these
// are self-contained, i.e., they don't need to talk to boundary conditions
// or interface objects or whatever.  Tell box what fine area you want
// filled, and it will tell you how much coarse data it needs.  Just
// fill the requested patch with coarse data and call fill.  How you
// get the coarse data is your problem.

class amr_interpolator_class 
{
public:
    amr_interpolator_class() 
    { 
    }
    virtual ~amr_interpolator_class() 
    { 
    }
    virtual Box box(const Box&, const IntVect&) const = 0;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const IntVect&) const = 0;
};

class bilinear_interpolator_class : public amr_interpolator_class 
{
public:
    bilinear_interpolator_class() 
    { 
    }
    Box box(const Box&, const IntVect&) const;
    void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const IntVect&) const;
};

#endif
