
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
    Mask& operator&= (const Mask& src) { return And(src); }
    //
    // in-place And
    //
    Mask& And (const Mask& src);
    //
    // as above, but specify source/destination/number of components
    //
    Mask& And (const Mask& src,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // as above, and specify subregion
    //
    Mask& And (const Mask& src,
               const Box&  subbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // in-pace And, over source components to destination components,
    //  and on intersection of source and destination boxes
    //
    Mask& And (const Mask& src,
               const Box&  srcbox,
               const Box&  destbox,
               int         srccomp,
               int         destcomp,
               int         numcomp = 1);
    //
    // in-place Or operator
    //
    Mask& operator|= (const Mask& src) { return Or(src); }
    //
    // in-place Or
    //
    Mask& Or (const Mask& src);
    //
    // as above, but specify source/destination/number of components
    //
    Mask& Or (const Mask& src,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
    //
    // as above, and specify subregion
    //
    Mask& Or (const Mask& src,
              const Box&  subbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
    //
    // in-pace Or, over source components to destination components,
    //  and on intersection of source and destination boxes
    //
    Mask& Or (const Mask& src,
              const Box&  srcbox,
              const Box&  destbox,
              int         srccomp,
              int         destcomp,
              int         numcomp = 1);
};

}

#endif /*_MASK_H_*/
