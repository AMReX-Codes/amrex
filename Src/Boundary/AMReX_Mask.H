
#ifndef AMREX_MASK_H_
#define AMREX_MASK_H_

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

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

class Mask
    :
    public BaseFab<int>
{
public:

    Mask () noexcept;

    explicit Mask (Arena* ar) noexcept;

    Mask (const Box& bx, int nc, Arena* ar);

    explicit Mask (const Box& bx,
		   int        nc = 1,
		   bool       alloc=true,
		   bool       shared=false,
                   Arena*     ar = nullptr);

    /**
    * \brief construct from input stream
    */
    explicit Mask (std::istream& is);

    explicit Mask (Array4<int> const& a) noexcept : BaseFab<int>(a) {}

    explicit Mask (Array4<int> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    explicit Mask (Array4<int const> const& a) noexcept : BaseFab<int>(a) {}

    explicit Mask (Array4<int const> const& a, IndexType t) noexcept : BaseFab<int>(a,t) {}

    ~Mask () noexcept {}

    Mask (Mask&& rhs) noexcept = default;

    Mask (Mask const& rhs, MakeType make_type, int scomp, int ncomp);

    Mask (const Mask&) = delete;
    Mask& operator= (const Mask&) = delete;
    Mask& operator= (Mask&&) = delete;


    /**
    * \brief Initialize from stream.
    *
    */
    friend std::istream& operator>> (std::istream&, Mask&);

    /**
    * \brief initialize from stream, FAB-style
    *
    */
    void readFrom (std::istream&);

    /**
    * \brief output to stream
    *
    */
    friend std::ostream& operator<< (std::ostream&, const Mask&);

    /**
    * \brief output to stream, FAB-style
    *
    */
    void writeOn (std::ostream&) const;
    //
    //! in-place And operator
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& operator&= (const Mask& src) noexcept { return And<run_on>(src); }

    /**
    * \brief in-place And
    *
    * \param src
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& And (const Mask& src) noexcept;

    /**
    * \brief as above, but specify source/destination/number of components
    *
    * \param src
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& And (const Mask& src,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;
    /**
    * \brief as above, and specify subregion
    *
    * \param src
    * \param subbox
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& And (const Mask& src,
               const Box&  subbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;
    /**
    * \brief in-pace And, over source components to destination components,
    * and on intersection of source and destination boxes
    *
    * \param src
    * \param srcbox
    * \param destbox
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& And (const Mask& src,
               const Box&  srcbox,
               const Box&  destbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1) noexcept;

    //! in-place Or operator
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& operator|= (const Mask& src) noexcept { return Or<run_on>(src); }

    /**
    * \brief in-place Or
    *
    * \param src
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& Or (const Mask& src) noexcept;

    /**
    * \brief as above, but specify source/destination/number of components
    *
    * \param src
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& Or (const Mask& src,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1) noexcept;
    /**
    * \brief as above, and specify subregion
    *
    * \param src
    * \param subbox
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Mask& Or (const Mask& src,
              const Box&  subbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1) noexcept;
    /**
    * \brief in-pace Or, over source components to destination components,
    * and on intersection of source and destination boxes
    *
    * \param src
    * \param srcbox
    * \param destbox
    * \param srccomp
    * \param destcomp
    * \param numcomp
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
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
