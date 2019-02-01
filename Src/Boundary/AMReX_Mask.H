
#ifndef AMREX_MASK_H_
#define AMREX_MASK_H_

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

namespace amrex {

/*
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
    //
    // default constructor
    //
    Mask ();
    //
    // allocating constructor (nc=number of components)
    //
    explicit Mask (const Box& bx,
		   int        nc = 1,
		   bool       alloc=true,
		   bool       shared=false);
    //
    // construct from input stream
    //
    explicit Mask (std::istream& is);

    ~Mask () = default;

    Mask (Mask&& rhs) noexcept = default;

#ifdef AMREX_USE_GPU
    Mask (Mask const& rhs, MakeType make_type);
#endif

    Mask (const Mask&) = delete;
    Mask& operator= (const Mask&) = delete;
    Mask& operator= (Mask&&) = delete;

    //
    // Initialize from stream.
    //
    friend std::istream& operator>> (std::istream&, Mask&);
    //
    // initialize from stream, FAB-style
    //
    void readFrom (std::istream&);
    //
    // output to stream
    //
    friend std::ostream& operator<< (std::ostream&, const Mask&);
    //
    // output to stream, FAB-style
    //
    void writeOn (std::ostream&) const;
    //
    // in-place And operator
    //
    AMREX_GPU_HOST_DEVICE
    Mask& operator&= (const Mask& src) { return And(src); }
    //
    // in-place And
    //
    AMREX_GPU_HOST_DEVICE
    Mask& And (const Mask& src);
    //
    // as above, but specify source/destination/number of components
    //
    AMREX_GPU_HOST_DEVICE
    Mask& And (const Mask& src,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // as above, and specify subregion
    //
    AMREX_GPU_HOST_DEVICE
    Mask& And (const Mask& src,
               const Box&  subbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // in-pace And, over source components to destination components,
    //  and on intersection of source and destination boxes
    //
    AMREX_GPU_HOST_DEVICE
    Mask& And (const Mask& src,
               const Box&  srcbox,
               const Box&  destbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // in-place Or operator
    //
    AMREX_GPU_HOST_DEVICE
    Mask& operator|= (const Mask& src) { return Or(src); }
    //
    // in-place Or
    //
    AMREX_GPU_HOST_DEVICE
    Mask& Or (const Mask& src);
    //
    // as above, but specify source/destination/number of components
    //
    AMREX_GPU_HOST_DEVICE
    Mask& Or (const Mask& src,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
    //
    // as above, and specify subregion
    //
    AMREX_GPU_HOST_DEVICE
    Mask& Or (const Mask& src,
              const Box&  subbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
    //
    // in-pace Or, over source components to destination components,
    //  and on intersection of source and destination boxes
    //
    AMREX_GPU_HOST_DEVICE
    Mask& Or (const Mask& src,
              const Box&  srcbox,
              const Box&  destbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::And (const Mask& src)
{
    return this->And(src,domain,domain,0,0,nvar);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::And (const Mask& src,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    return this->And(src,domain,domain,srccomp,destcomp,numcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::And (const Mask& src,
           const Box&  subbox,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    return this->And(src,subbox,subbox,srccomp,destcomp,numcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::And (const Mask& src,
           const Box&  srcbox,
           const Box&  destbox,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = dp(i,j,k,n) ? sp(i,j,k,n) : 0;
                }
            }
        }
    }

    return *this;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::Or (const Mask& src)
{
    return this->Or(src,domain,domain,0,0,nvar);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::Or (const Mask& src,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    return this->Or(src,domain,domain,srccomp,destcomp,numcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::Or (const Mask& src,
          const Box&  subbox,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    return this->Or(src,subbox,subbox,srccomp,destcomp,numcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Mask&
Mask::Or (const Mask& src,
          const Box&  srcbox,
          const Box&  destbox,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = dp(i,j,k,n) ? 1: sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

}

#endif /*_MASK_H_*/
