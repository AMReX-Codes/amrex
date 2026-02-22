
#ifndef AMREX_MASK_H_
#define AMREX_MASK_H_
#include <AMReX_Config.H>

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

/**
 * \file AMReX_Mask.H
 * \brief Logical mask type used to classify regions of FABs.
 */

namespace amrex {

/**
        A Mask is a derived Fab class whose elements maintain a logical
        "TRUE" or "FALSE" state.  A Mask can perform several elementary
        logical operations on its data members.

        Implement a logical class Mask to maintain the TRUE or FALSE state
        of elements in a box.  A Mask can perform element-wise logical
        operations, but is otherwise very similar to an FArrayBox object.
        A Mask object can set its components based on comparing (a subregion
        of a) real-valued FArrayBox element magnitudes to a scalar value, or
        by "and", "or" operations with (the a subregion of) other Masks.

        This class does NOT provide a copy constructor or assignment operator.
*/
class Mask final
    :
    public BaseFab<int>
{
public:

    //! Construct an empty mask with no allocation.
    Mask () noexcept = default;

    //! Construct a mask that allocates from Arena \p ar.
    explicit Mask (Arena* ar) noexcept;

    //! Construct a mask over Box \p bx with \p nc components using Arena \p ar.
    Mask (const Box& bx, int nc, Arena* ar);

    /**
     * \brief Construct a mask with fine-grained control over allocation.
     *
     * \param bx     Region covered by the mask.
     * \param nc     Number of components (default 1).
     * \param alloc  Whether to allocate storage immediately.
     * \param shared Whether storage is shared by reference.
     * \param ar     Arena used for allocations (or default if nullptr).
     */
    explicit Mask (const Box& bx,
                   int        nc = 1,
                   bool       alloc=true,
                   bool       shared=false,
                   Arena*     ar = nullptr);

    //! Construct from an input stream \p is.
    explicit Mask (std::istream& is);

    //! Wrap an existing Array4 \p a without taking ownership.
    explicit Mask (Array4<int> const& a) noexcept : BaseFab<int>(a) {}

    //! Wrap an existing Array4 \p a with IndexType \p t.
    explicit Mask (Array4<int> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    //! Wrap a const Array4 \p a without taking ownership.
    explicit Mask (Array4<int const> const& a) noexcept : BaseFab<int>(a) {}

    //! Wrap a const Array4 \p a with IndexType \p t.
    explicit Mask (Array4<int const> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    ~Mask () noexcept override = default;

    Mask (Mask&& rhs) noexcept = default;

    //! Construct by copying \p ncomp components starting at \p scomp from \p rhs using \p make_type semantics.
    Mask (Mask const& rhs, MakeType make_type, int scomp, int ncomp);

    Mask (const Mask&) = delete;
    Mask& operator= (const Mask&) = delete;
    Mask& operator= (Mask&&) = delete;


    /**
     * \brief Initialize from stream in ASCII FAB format.
     *
     * \param is Input stream.
     * \param m  Mask populated from \p is.
     * \return The input stream for chaining.
     */
    friend std::istream& operator>> (std::istream& is, Mask& m);

    /**
     * \brief Initialize from stream, FAB-style.
     *
     * \param is Input stream.
     */
    void readFrom (std::istream& is);

    /**
     * \brief Output to stream in ASCII format.
     *
     * \param os Output stream receiving the data.
     * \param m  Mask whose contents will be written.
     */
    friend std::ostream& operator<< (std::ostream& os, const Mask& m);

    /**
     * \brief Output to stream in FAB binary format.
     *
     * \param os Output stream.
     */
    void writeOn (std::ostream& os) const;
    //
    //! In-place logical AND with mask \p src.
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& operator&= (const Mask& src) noexcept { return And<run_on>(src); }

    /**
    * \brief In-place logical AND over every component.
    *
    * \param src Source mask supplying the AND operand.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& And (const Mask& src) noexcept;

    /**
    * \brief As above, but specify source/destination component ranges.
    *
    * \param src      Source mask.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components to combine.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& And (const Mask& src,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;
    /**
    * \brief As above, but restrict operation to a subbox.
    *
    * \param src      Source mask.
    * \param subbox   Region where the operation is performed.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& And (const Mask& src,
               const Box&  subbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;
    /**
    * \brief In-place AND on potentially different source/destination boxes.
    *
    * \param src      Source mask.
    * \param srcbox   Source region.
    * \param destbox  Destination region.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& And (const Mask& src,
               const Box&  srcbox,
               const Box&  destbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;

    //! In-place logical OR with mask \p src.
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& operator|= (const Mask& src) noexcept { return Or<run_on>(src); }

    /**
    * \brief In-place logical OR over every component.
    *
    * \param src Source mask supplying the OR operand.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& Or (const Mask& src) noexcept;

    /**
    * \brief As above, but specify source/destination component ranges.
    *
    * \param src      Source mask.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& Or (const Mask& src,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1) noexcept;
    /**
    * \brief As above, and specify a subregion.
    *
    * \param src      Source mask.
    * \param subbox   Region where to apply the OR.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& Or (const Mask& src,
              const Box&  subbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1) noexcept;
    /**
    * \brief In-place OR allowing separate source and destination boxes.
    *
    * \param src      Source mask.
    * \param srcbox   Source region.
    * \param destbox  Destination region.
    * \param srccomp  Starting component in \p src.
    * \param destcomp Starting component in this mask.
    * \param numcomp  Number of components.
    */
    template <RunOn run_on AMREX_DEFAULT_RUNON>
    Mask& Or (const Mask& src,
              const Box&  srcbox,
              const Box&  destbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1) noexcept;
};

template <RunOn run_on>
Mask&
Mask::And (const Mask& src) noexcept
{
    return this->And<run_on>(src,domain,domain,0,0,nvar);
}

template <RunOn run_on>
Mask&
Mask::And (const Mask& src,
           int         srccomp,
           int         destcomp,
           int         numcomp) noexcept
{
    return this->And<run_on>(src,domain,domain,srccomp,destcomp,numcomp);
}

template <RunOn run_on>
Mask&
Mask::And (const Mask& src,
           const Box&  subbox,
           int         srccomp,
           int         destcomp,
           int         numcomp) noexcept
{
    return this->And<run_on>(src,subbox,subbox,srccomp,destcomp,numcomp);
}

template <RunOn run_on>
Mask&
Mask::And (const Mask& src,
           const Box&  srcbox,
           const Box&  destbox,
           int         srccomp,
           int         destcomp,
           int         numcomp) noexcept
{
    auto const& d = this->array();
    auto const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) = d(i,j,k,n+destcomp) ? s(i+offset.x,j+offset.y,k+offset.z,n+srccomp) : 0;
    });
    return *this;
}

template <RunOn run_on>
Mask&
Mask::Or (const Mask& src) noexcept
{
    return this->Or<run_on>(src,domain,domain,0,0,nvar);
}

template <RunOn run_on>
Mask&
Mask::Or (const Mask& src,
          int         srccomp,
          int         destcomp,
          int         numcomp) noexcept
{
    return this->Or<run_on>(src,domain,domain,srccomp,destcomp,numcomp);
}

template <RunOn run_on>
Mask&
Mask::Or (const Mask& src,
          const Box&  subbox,
          int         srccomp,
          int         destcomp,
          int         numcomp) noexcept
{
    return this->Or<run_on>(src,subbox,subbox,srccomp,destcomp,numcomp);
}

template <RunOn run_on>
Mask&
Mask::Or (const Mask& src,
          const Box&  srcbox,
          const Box&  destbox,
          int         srccomp,
          int         destcomp,
          int         numcomp) noexcept
{
    auto const& d = this->array();
    auto const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) = d(i,j,k,n+destcomp) ? 1: s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });
    return *this;
}

}

#endif /*_MASK_H_*/
