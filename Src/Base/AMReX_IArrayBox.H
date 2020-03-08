
#ifndef BL_IFARRAYBOX_H
#define BL_IFARRAYBOX_H

#include <AMReX_Box.H>
#include <AMReX_BaseFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

namespace amrex {

/**
* \brief  A Fortran Array of ints

*  iFortran Array Box's (generally called iFAB's) are objects constructed
*  to emulate the FORTRAN array.  Useful operations can be performed
*  upon FAB's in C++, and they provide a convenient interface to
*  FORTRAN when it is necessary to retreat into that language.

*  IArrayBox is derived from BaseFab<int>.
*  IArrayBox adds additional useful capabilities which make sense
*  for int types, such as L**p norms.

*  The C pre-processor macro AMREX_SPACEDIM must be defined to use
*  this class.  The internal precision of FARRAYBOX objects is
*  set by defining either BL_USE_FLOAT or BL_USE_DOUBLE

*  This is NOT a polymorphic class.

*  This class does NOT provide a copy constructor or assignment operator.
*/

class IArrayBox
    :
    public BaseFab<int>
{
public:
    //! Construct an invalid FAB with no memory.
    IArrayBox () noexcept;

    explicit IArrayBox (Arena* ar) noexcept;

    IArrayBox (const Box& b, int ncomp, Arena* ar);

    /**
    * \brief Construct an initial FAB with the data space allocated but
    * not inititialized. ncomp is the number of components
    * (variables) at each data point in the Box.
    */
    explicit IArrayBox (const Box& b,
                        int        ncomp=1,
			bool       alloc=true,
			bool       shared=false,
                        Arena*     ar = nullptr);

    IArrayBox (const IArrayBox& rhs, MakeType make_type, int scomp, int ncomp);

    explicit IArrayBox (Array4<int> const& a) noexcept : BaseFab<int>(a) {}

    explicit IArrayBox (Array4<int> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    explicit IArrayBox (Array4<int const> const& a) noexcept : BaseFab<int>(a) {}

    explicit IArrayBox (Array4<int const> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    //!  The destructor.
    ~IArrayBox () noexcept {}

    IArrayBox (IArrayBox&& rhs) noexcept = default;

    IArrayBox (const IArrayBox&) = delete;
    IArrayBox& operator= (const IArrayBox&) = delete;
    IArrayBox& operator= (IArrayBox&&) = delete;

    //! Set the fab to the value r.
    template <RunOn run_on>
    IArrayBox& operator= (int r) noexcept;

    //! For debugging purposes we hide BaseFab version and do some extra work.
    void resize (const Box& b, int N = 1);

    //! Initialize from ParmParse with "fab" prefix.
    static void Initialize ();
    static void Finalize ();

private:

    static bool do_initval;

};

template <RunOn run_on>
IArrayBox&
IArrayBox::operator= (int v) noexcept
{
    BaseFab<int>::operator=<run_on>(v);
    return *this;
}

}

#endif /*BL_IFARRAYBOX_H*/
