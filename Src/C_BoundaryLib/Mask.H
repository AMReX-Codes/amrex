#ifndef _MASK_H_
#define _MASK_H_
//
// $Id: Mask.H,v 1.1 1998-03-24 07:06:35 almgren Exp $
//

#include <BaseFab.H>
#include <FArrayBox.H>

//@Man:
/*@Memo:
        A Mask is a derived Fab class whose elements maintain a logical
        "TRUE" or "FALSE" state.  A Mask can perform several elementary
        logical operations on its data members.
*/        
/*@Doc:
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
    public ArithFab<int>,
    public OrderedFab<int>
{
public:

//@ManMemo: administrative functions
//@ManDoc: default constructor
    Mask()
        : BaseFab<int>() {};
//@ManDoc: allocating constructor (nc=number of components)
    Mask(const Box& bx, int nc = 1)
        : BaseFab<int>(bx,nc) {};
//@ManDoc: construct from input stream
    Mask(istream &is);
//@ManDoc: destructor
    virtual ~Mask() {};

//@ManMemo: I/O functions
//@ManDoc: initialize from stream
    friend istream& operator >> (istream&, Mask&);
//@ManDoc: initialize from stream, FAB-style
    void readFrom(istream&);
//@ManDoc: output to stream
    friend ostream& operator << (ostream&, const Mask&);
//@ManDoc: output to stream, FAB-style
    void writeOn(ostream&) const;


//@ManMemo: logical And functions
//@ManDoc: in-place And operator
    Mask& operator &= (const Mask& src) { return And(src); };
//@ManDoc: in-place And
    Mask& And(const Mask& src);
//@ManDoc: as above, but specify source/destination/number of components
    Mask& And(const Mask& src,
              int srccomp, int destcomp, int numcomp=1);
//@ManDoc: as above, and specify subregion
    Mask& And(const Mask& src, const Box& subbox,
              int srccomp, int destcomp, int numcomp=1);
//@ManDoc: in-pace And, over source components to destination components, and on intersection of source and destination boxes
    Mask& And(const Mask& src, const Box& srcbox,
              const Box& destbox,
              int srccomp, int destcomp, int numcomp=1);

//@ManMemo: logical Or functions
//@ManDoc: in-place Or operator
    Mask& operator |= (const Mask& src) { return Or(src); };
//@ManDoc: in-place Or
    Mask& Or(const Mask& src);
//@ManDoc: as above, but specify source/destination/number of components
    Mask& Or(const Mask& src,
             int srccomp, int destcomp, int numcomp=1);
//@ManDoc: as above, and specify subregion
    Mask& Or(const Mask& src, const Box& subbox,
             int srccomp, int destcomp, int numcomp=1);
//@ManDoc: in-pace Or, over source components to destination components, and on intersection of source and destination boxes
    Mask& Or(const Mask& src, const Box& srcbox,
             const Box& destbox,
             int srccomp, int destcomp, int numcomp=1);
private:
// Disallow copy operator
    Mask& operator=(const Mask&);
// Disallow copy constructor
    Mask(const Mask&);
};

#endif /*_MASK_H_*/
