//BL_COPYRIGHT_NOTICE

#ifndef BL_BOXARRAY_H
#define BL_BOXARRAY_H

//
// $Id: BoxArray.H,v 1.14 1999-09-30 15:55:12 sstanley Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstddef>
#else
#include <stddef.h>
#endif

#include <IndexType.H>
#include <BoxList.H>
#include <Array.H>
#include <Pointers.H>

//
// Forward declaration.
//
class BoxList;

class BoxArray;
BoxArray boxComplement (const Box& b1in, const Box& b2);
BoxArray complementIn (const Box& b, const BoxArray& ba);
BoxArray intersect (const BoxArray& ba, const Box& b);

//
//@Man:
//@Memo: An Array of Boxes
/*@Doc:

  A BoxArray is a collection of Boxes stored in an Array.  It is a
  reference-counted concrete class, not a polymorphic one; i.e. you
  cannot use any of the List member functions with a BoxList.
*/

class BoxArray
{
public:
    //
    //@ManDoc: Write a BoxArray to an ostream in ASCII format.
    //
    friend ostream& operator<< (ostream&        os,
                                const BoxArray& ba);
    //
    //@ManDoc: Construct an empty BoxArray.
    //
    BoxArray ();
    //
    //@ManDoc: Construct a BoxArray of the specified size.
    //
    explicit BoxArray (size_t size);
    //
    //@ManDoc: Construct a BoxArray from an array of Boxes of size nbox.
    //
    BoxArray (const Box* bxvec,
              int        nbox);
    //
    //@ManDoc: Construct a BoxArray from a BoxList.
    //
    explicit BoxArray (const BoxList& bl);
    //
    //@ManDoc: The copy constructor.
    //
    BoxArray (const BoxArray& bs);
    //
    //@ManDoc: Construct a BoxArray from a checkpoint file.
    //
    explicit BoxArray (istream& is);
    //
    //@ManDoc: The assignment operator.
    //
    BoxArray& operator= (const BoxArray& rhs);
    //
    //@ManDoc: The destructor.
    //
    ~BoxArray ();

    /*@ManDoc: Initialize the BoxArray from the supplied BoxList.
               It is an error if the BoxArray has already been initialized.
    */
    void define (const BoxList& bl);

    /*@ManDoc: Initialize the BoxArray from the supplied BoxArray.
               It is an error if the BoxArray has already been initialized.
    */
    void define (const BoxArray& bs);

    /*@ManDoc: Initialize the BoxArray from the supplied istream.
               It is an error if the BoxArray has already been initialized.
               Note that the BoxArray in the istream must have been written
               using writeOn().
    */
    void define (istream& is);
    //
    //@ManDoc: Remove all Boxes from the BoxArray.
    //
    void clear ();
    //
    //@ManDoc: Resizes the array.  See Array<T>::resize() for the gory details.
    //
    void resize (int len);
    //
    //@ManDoc: Returns the number of elements in the BoxArray.
    //
    int length () const;
    //
    //@ManDoc: Returns true if the BoxArray is not empty.
    //
    bool ready () const;
    //
    //@ManDoc: Output this BoxArray to a checkpoint file.
    //
    ostream& writeOn (ostream&) const;
    //
    //@ManDoc: Are the BoxArrays equal?
    //
    bool operator== (const BoxArray& rhs) const;
    //
    //@ManDoc: Are the BoxArrays not equal?
    //
    bool operator!= (const BoxArray& rhs) const;
    //
    //@ManDoc: Forces each Box in BoxArray to have sides <= block\_size.
    //
    BoxArray& maxSize (int block_size);
    //
    //@ManDoc: Refine each Box in the BoxArray to the specified ratio.
    //
    BoxArray& refine (int refinement_ratio);
    //
    //@ManDoc: Refine each Box in the BoxArray to the specified ratio.
    //
    BoxArray& refine (const IntVect& refinement_ratio);
    //
    //@ManDoc: Coarsen each Box in the BoxArray to the specified ratio.
    //
    BoxArray& coarsen (int refinement_ratio);
    //
    //@ManDoc: Coarsen each Box in the BoxArray to the specified ratio.
    //
    BoxArray& coarsen (const IntVect& refinement_ratio);
    //
    //@ManDoc: Grow each Box in the BoxArray by the specified amount.
    //
    BoxArray& grow (int n);
    //
    //@ManDoc: Grow each Box in the BoxArray by the specified amount.
    //
    BoxArray& grow (const IntVect& iv);

    /*@ManDoc: Grow each Box in the BoxArray on the low and high end
               by n\_cell cells in idir direction.
    */
    BoxArray& grow (int idir,
                    int n_cell);

    /*@ManDoc: Applies surroundingNodes(Box) to each Box in BoxArray.
               See the documentation of Box for details.
    */
    BoxArray& surroundingNodes ();

    /*@ManDoc: Applies surroundingNodes(Box,int) to each Box in
               BoxArray.  See the documentation of Box for details.
    */
    BoxArray& surroundingNodes (int dir);
    //
    //@ManDoc: Applies Box::enclosedCells() to each Box in the BoxArray.
    //
    BoxArray& enclosedCells ();
    //
    //@ManDoc: Applies Box::enclosedCells(int) to each Box in the BoxArray.
    //
    BoxArray& enclosedCells  (int dir);
    //
    //@ManDoc: Applies Box::convert(IndexType) to each Box in the BoxArray.
    //
    BoxArray& convert (IndexType typ);
    //
    //@ManDoc: Applies function (*fp)(Box) to each Box in the BoxArray.
    //
    BoxArray& convert (Box (*fp)(const Box&));
    //
    //@ManDoc: Applies Box::shift(int,int) to each Box in the BoxArray.
    //
    BoxArray& shift (int dir,
                     int nzones);
    //
    //@ManDoc: Applies Box::shiftHalf(int,int) to each Box in the BoxArray.
    //
    BoxArray& shiftHalf (int dir,
                         int num_halfs);
    //
    //@ManDoc: Applies Box::shiftHalf(IntVect) to each Box in BoxArray.
    //
    BoxArray& shiftHalf (const IntVect &iv);
    //
    //@ManDoc: Set element i in this BoxArray to Box ibox.
    //
    void set (int        i,
              const Box& ibox);
    //
    //@ManDoc: Returns element index of this BoxArray.
    //
    const Box& operator[] (int index) const;
    //
    //@ManDoc: Returns element index of this BoxArray.
    //
    const Box& get (int index) const;

    /*@ManDoc: Returns true if Box is valid and they all have the same
               IndexType.  Is true by default if the BoxArray is empty.
    */
    bool ok () const;
    //
    //@ManDoc: Returns true if set of intersecting Boxes in BoxArray is null.
    //
    bool isDisjoint () const;
    //
    //@ManDoc: Create a BoxList from this BoxArray.
    //
    BoxList boxList () const;
    //
    //@ManDoc: True if the IntVect is within any of the Boxes in this BoxArray.
    //
    bool contains (const IntVect& v) const;

    /*@ManDoc: True if the Box is within any of the Boxes in the List.
               The Box must also have the same IndexType as those in this
               BoxArray.
    */
    bool contains (const Box& b) const;
    //
    //@ManDoc: True if all Boxes in bl are contained in this BoxArray.
    //
    bool contains (const BoxArray& bl) const;
    //
    //@ManDoc: Returns smallest Box that contains all Boxes in this BoxArray.
    //
    Box minimalBox () const;
    //
    //@ManDoc: Makes a BoxArray from the the complement of b2 in b1in.
    //
    friend BoxArray boxComplement (const Box& b1in,
                                   const Box& b2);
    //
    //@ManDoc: Makes a BoxArray from the complement of BoxArray ba in Box b.
    //
    friend BoxArray complementIn (const Box&      b,
                                  const BoxArray& ba);
    //
    //@ManDoc: Makes a BoxArray from the intersection of Box b and BoxArray ba.
    //
    friend BoxArray intersect (const BoxArray& ba,
                               const Box&      b);

    void reserve (long _truesize);

protected:
    //
    // Helper for ASCII stream output.
    //
    ostream& print (ostream& os) const;

    class Ref
    {
        friend class BoxArray;
        //
        // Constructors to match those in BoxArray ....
        //
        Ref ();

        Ref (size_t size);

        Ref (const Box* bxvec,
             int        nbox);

        Ref (const BoxList& bl);

        Ref (istream& is);

        Ref (const Ref& rhs);
        //
        // Some defines()s to match those in BoxArray.
        //
        void define (const BoxList& bl);

        void define (istream& is);
        //
        // The data -- an Array of Boxes.
        //
        Array<Box> m_abox;   
    };
    //
    // Make `this' unique.
    //
    void uniqify ();
    //
    // The data -- a reference-counted pointer to a Ref.
    //
    LnClassPtr<Ref> m_ref;
};

//
// Inlines.
//

inline
int
BoxArray::length () const
{
    return m_ref->m_abox.length();
}

inline
bool
BoxArray::ready () const
{
    return m_ref->m_abox.ready();
}

inline
const Box&
BoxArray::operator[] (int index) const
{
    return m_ref->m_abox.get(index);
}

inline
const Box&
BoxArray::get (int index) const
{
    return m_ref->m_abox.get(index);
}

inline
void
BoxArray::reserve (long _truesize)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.reserve(_truesize);
}

#endif /*BL_BOXARRAY_H*/
