
#ifndef BL_IFARRAYBOX_H
#define BL_IFARRAYBOX_H

#include <Pointers.H>
#include <Box.H>
#include <BaseFab.H>
#include <REAL.H>
#include <SPACE.H>

class IArrayBox;

/*
  A Fortran Array of ints

  iFortran Array Box's (generally called iFAB's) are objects constructed
  to emulate the FORTRAN array.  Useful operations can be performed
  upon FAB's in C++, and they provide a convenient interface to
  FORTRAN when it is necessary to retreat into that language.

  IArrayBox is derived from BaseFab<int>.
  IArrayBox adds additional useful capabilities which make sense
  for int types, such as L**p norms.

  The C pre-processor macro BL_SPACEDIM must be defined to use
  this class.  The internal precision of FARRAYBOX objects is
  set by defining either BL_USE_FLOAT or BL_USE_DOUBLE

  This is NOT a polymorphic class.

  This class does NOT provide a copy constructor or assignment operator.
*/

class IArrayBox
    :
    public BaseFab<int>
{
public:
    //
    // Construct an invalid FAB with no memory.
    //
    IArrayBox ();
    //
    // Construct an initial FAB with the data space allocated but
    // not inititialized. ncomp is the number of components
    // (variables) at each data point in the Box.
    //
    explicit IArrayBox (const Box& b,
                        int        ncomp=1,
			bool       alloc=true,
			bool       shared=false);
    //
    // Set the fab to the value r.
    //
    IArrayBox& operator= (const int& r);
    //
    //  The destructor.
    //
    ~IArrayBox ();

    //
    // For debugging purposes we hide BaseFab version and do some extra work.
    //
    void resize (const Box& b,
                 int        N = 1);
    //
    // Compute the Lp-norm of this FAB using components
    // (scomp : scomp+ncomp-1).  p < 0  -> ERROR.
    // p = 0  -> infinity norm (max norm).
    // p = 1  -> sum of ABS(FAB)
    // p > 1  -> Lp-norm
    //
    int norm (int p     = 2,
              int scomp = 0,
              int ncomp = 1) const;
    //
    // Same as above except only on given subbox.
    //
    int norm (const Box& subbox,
              int        p     = 2,
              int        scomp = 0,
              int        ncomp = 1) const;

    // 
    // Initialize from ParmParse with "fab" prefix.
    //
    static void Initialize ();
    static void Finalize ();

private:

    static bool do_initval;

    //
    // These are not allowed.
    //
    IArrayBox (const IArrayBox&);
    IArrayBox& operator= (const IArrayBox&);
};

#endif /*BL_IFARRAYBOX_H*/
