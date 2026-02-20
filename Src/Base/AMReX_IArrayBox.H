
#ifndef BL_IFARRAYBOX_H
#define BL_IFARRAYBOX_H
#include <AMReX_Config.H>

/**
 * \file AMReX_IArrayBox.H
 * \brief Integer-valued FAB (Field Array Box) type.
 */

#include <AMReX_Box.H>
#include <AMReX_BaseFab.H>
#include <AMReX_SPACE.H>

#include <iosfwd>

namespace amrex {

class IArrayBox;
class IntDescriptor;

/**
 * \class amrex::IFABio
 * \ingroup amrex_fab_containers
 * \brief Low-level I/O helper for IArrayBox; handles header writing and binary reading.
 */
class IFABio
{
public:
    //! Write the IArrayBox header for \p nvar components to \p os.
    static void write_header (std::ostream& os, const IArrayBox& fab, int nvar) ;
    //! Read IArrayBox data from \p is using the given \p data_descriptor.
    static void read (std::istream& is, IArrayBox& fab, IntDescriptor const& data_descriptor) ;
};

/**
* \class amrex::IArrayBox
* \ingroup amrex_fab_containers
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
    friend class IFABio;

    //! Construct an invalid FAB with no memory.
    IArrayBox () noexcept = default;

    //! Construct an empty FAB that allocates from arena \p ar.
    explicit IArrayBox (Arena* ar) noexcept;

    //! Construct a FAB over \p b with \p ncomp components, allocated from \p ar.
    IArrayBox (const Box& b, int ncomp, Arena* ar);

    /**
    * \brief Construct an initial FAB with the data space allocated but
    * not initialized. ncomp is the number of components
    * (variables) at each data point in the Box.
    */
    explicit IArrayBox (const Box& b,
                        int        ncomp=1,
                        bool       alloc=true,
                        bool       shared=false,
                        Arena*     ar = nullptr);

    //! Construct an alias FAB viewing components [\p scomp, \p scomp+\p ncomp) of \p rhs.
    IArrayBox (const IArrayBox& rhs, MakeType make_type, int scomp, int ncomp);

    explicit IArrayBox (Array4<int> const& a) noexcept : BaseFab<int>(a) {}

    explicit IArrayBox (Array4<int> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    explicit IArrayBox (Array4<int const> const& a) noexcept : BaseFab<int>(a) {}

    explicit IArrayBox (Array4<int const> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    //!  The destructor.
    ~IArrayBox () noexcept override = default;

    IArrayBox (IArrayBox&& rhs) noexcept = default;
    IArrayBox& operator= (IArrayBox&&) = default;

    IArrayBox (const IArrayBox&) = delete;
    IArrayBox& operator= (const IArrayBox&) = delete;

    //! Set the fab to the value r.
    template <RunOn run_on>
    IArrayBox& operator= (int v) noexcept;

    //! For debugging purposes we hide BaseFab version and do some extra work.
    void resize (const Box& b, int N = 1, Arena* ar = nullptr);

    //! Read FAB data from binary stream \p is.
    void readFrom (std::istream& is);

    //! Initialize the IArrayBox subsystem.
    static void Initialize ();
    //! Finalize the IArrayBox subsystem.
    static void Finalize ();

    //! Return the IntDescriptor describing the native integer format.
    static std::unique_ptr<IntDescriptor> getDataDescriptor ();

    //! Return the I/O helper object.
    static IFABio const& getFABio ();

    //! Return the class name string used for I/O.
    static std::string getClassName ();

private:

    static bool do_initval;
    static std::unique_ptr<IFABio> ifabio;

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
