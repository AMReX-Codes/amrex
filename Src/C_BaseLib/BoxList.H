//BL_COPYRIGHT_NOTICE

#ifndef BL_BOXLIST_H
#define BL_BOXLIST_H

//
// $Id: BoxList.H,v 1.3 1997-09-18 20:12:46 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#include <Boolean.H>
#include <IntVect.H>
#include <IndexType.H>
#include <Box.H>
#include <BoxArray.H>
#include <List.H>

class BoxArray;
class BoxListIterator;

//
//@Man:
//@Memo: A List of Boxes
/*@Doc:

  A BoxList is a class for managing a List of Boxes that share a common
  IndexType.  This class implements operations for sets of Boxes.

  This is a concrete class, not a polymorphic one.
*/

class BoxList
{
    friend class BoxListIterator;
    friend class BoxDomain;

public:
    //
    //@ManDoc: Construct an empty BoxList with IndexType::TheCellType().
    //
    BoxList ();
    //
    //@ManDoc: Construct an empty BoxList with IndexType btyp
    //
    BoxList (IndexType btyp);
    //
    //@ManDoc:  Construct a BoxList from a BoxArray.
    //
    BoxList (const BoxArray& ba);
    //
    //@ManDoc: The copy constructor.
    //
    BoxList (const BoxList& _blst);
    //
    //@ManDoc: The assignment operator.
    //
    BoxList& operator= (const BoxList& rhs);
    //
    //@ManDoc: The destructor.
    //
    ~BoxList ();
    //
    //@ManDoc: Append a Box to this BoxList.
    //
    void append (const Box& bn);
    //
    //@ManDoc: Add a Box to this BoxList.
    //
    void add (const Box& bn);
    //
    //@ManDoc:  Prepend a Box to this BoxList.
    //
    void prepend (const Box& bn);
    //
    //@ManDoc: Join the BoxList to ourselves.
    //
    void join (const BoxList& blist);
    //
    //@ManDoc: Add the Box after where the iterator points.
    //
    void addAfter (BoxListIterator& bli,
                   const Box&       bn);
    //
    //@ManDoc: Add the Box before where the iterator points.
    //
    void addBefore (BoxListIterator& bli,
                    const Box&       bn);
    //
    //@ManDoc: Remove all Boxes from this BoxList.
    //
    void clear ();
    //
    //@ManDoc: The number of Boxes in this BoxList.
    //
    int length () const;

    /*@ManDoc: True if this BoxList is valid; i.e. all the Boxes are
               valid and they all have the same IndexType.  Also returns
               true if the BoxList is empty.
    */
    bool ok () const;
    //
    //@ManDoc: Is this BoxList equal to rhs?
    //
    bool operator== (const BoxList& rhs) const;
    //
    //@ManDoc: Is this BoxList notequal to rhs?
    //
    bool operator!= (const BoxList& rhs) const;
    //
    //@ManDoc: Is this BoxList empty?
    //
    bool isEmpty () const;
    //
    //@ManDoc: Is this BoxList not empty?
    //
    bool isNotEmpty () const;
    //
    //@ManDoc: Access the Box to which the iterator points.
    //
    const Box& operator[] (const BoxListIterator& bli) const;
    //
    //@ManDoc: Access the Box to which the iterator points.
    //
    Box& operator[] (const BoxListIterator& bli);
    //
    //@ManDoc: True if the set of intersecting Boxes is empty.
    //
    bool isDisjoint () const;
    //
    //@ManDoc: True if the IntVect is within any of the Boxes in this BoxList.
    //
    bool contains (const IntVect& v) const;

    /*@ManDoc: True if the Box is within any of the Boxes in the List.
               The Box must also have the same IndexType as those in this
               BoxList.
    */
    bool contains (const Box& b) const;
    //
    //@ManDoc: True if all Boxes in bl are contained in this BoxList.
    //
    bool contains (const BoxList& bl) const;

    /*@ManDoc: Modify this BoxList to contain only its intersection with
               Box b.
    */
    BoxList& intersect (const Box& b);

    /*@ManDoc: Modify this BoxList to contain only its intersection with
               BoxList bl.
    */
    BoxList& intersect (const BoxList &bl);
    //
    //@ManDoc: Remove the Box from this BoxList.
    //
    BoxList& remove (const Box& bx);
    //
    //@ManDoc: Remove the pointed to Box from this BoxList.
    //
    BoxList& remove (BoxListIterator& bli);
    //
    //@ManDoc: Creates the complement of BoxList bl in Box b.
    //
    BoxList& complementIn (const Box&     b,
                           const BoxList& bl);
    //
    //@ManDoc: Refine each Box in the BoxList by the ratio.
    //
    BoxList& refine (int ratio);
    //
    //@ManDoc: Refine each Box in the BoxList by the ratio.
    //
    BoxList& refine (const IntVect& ratio);
    //
    //@ManDoc: Coarsen each Box in the BoxList by the ratio.
    //
    BoxList& coarsen (int ratio);
    //
    //@ManDoc: Coarsen each Box in the BoxList by the ratio.
    //
    BoxList& coarsen (const IntVect& ratio);
    //
    //@ManDoc: Grow each Box in the BoxList by size sz.
    //
    BoxList& accrete  (int sz);
    //
    //@ManDoc: Applies Box::shift(int,int) to each Box in the BoxList.
    //
    BoxList& shift (int dir,
                    int nzones);
    //
    //@ManDoc: Applies Box::shiftHalf(int,int) to each Box in the BoxList.
    //
    BoxList& shiftHalf (int dir,
                        int num_halfs);
    //
    //@ManDoc: Applies Box::shiftHalf(IntVect) to each Box in BoxList.
    //
    BoxList& shiftHalf (const IntVect& iv);

    /*@ManDoc: Merge adjacent Boxes in this BoxList. Return the number
               of Boxes merged.
    */
    int simplify ();

    /*@ManDoc: Merges all possible adjacent Boxes.  Returns the total
               number of Boxes merged.
    */
    int minimize ();

    /*@ManDoc: Forces each Box in the BoxList to have sides of 
               length <= chunk.
    */
    BoxList& maxSize (int chunk);
    //
    //@ManDoc: Returns smallest Box that contains all Boxes in this BoxList.
    //
    Box minimalBox () const;
    //
    //@ManDoc: Returns the IndexType of Boxes in this BoxList.
    //
    IndexType ixType () const;

    /*@ManDoc: Applies surroundingNodes(Box) to each Box in BoxArray.
               See the documentation of Box for details.
    */
    BoxList& surroundingNodes ();

    /*@ManDoc: Applies surroundingNodes(Box,int) to each Box in
               BoxList.  See the documentation of Box for details.
    */
    BoxList& surroundingNodes (int dir);
    //
    //@ManDoc: Applies Box::enclosedCells() to each Box in the BoxList.
    //
    BoxList& enclosedCells ();
    //
    //@ManDoc: Applies Box::enclosedCells(int) to each Box in the BoxList.
    //
    BoxList& enclosedCells  (int dir);
    //
    //@ManDoc: Applies Box::convert(IndexType) to each Box in the BoxList.
    //
    BoxList& convert (IndexType typ);
    //
    //@ManDoc: Output a BoxList to an ostream in ASCII format.
    //
    friend ostream& operator<< (ostream&       os,
                                const BoxList& blist);
    //
    //@ManDoc: Returns a BoxList defining the intersection of bl with b.
    //
    friend inline BoxList intersect (const BoxList& bl,
                                     const Box&     b);

    /*@ManDoc: Returns a new BoxList in which each Box is refined by the
               given ratio.
    */
    friend inline BoxList refine (const BoxList&  bl,
                                  int             ratio);

    /*@ManDoc: Returns a new BoxList in which each Box is coarsened by
               the given ratio.
    */
    friend inline BoxList coarsen (const BoxList& bl,
                                   int            ratio);

    /*@ManDoc: Returns a new BoxList in which each Box is grown by the
               given size.
    */
    friend inline BoxList accrete (const BoxList& bl,
                                   int            sz);
    //
    //@ManDoc: Returns BoxList defining the compliment of b2 in b1in.
    //
    friend BoxList boxDiff (const Box& b1in,
                            const Box& b2);

    /*@ManDoc: Returns a BoxList defining the complement of BoxList bl
               in Box b.
    */
    friend BoxList complementIn  (const Box&     b,
                                  const BoxList& bl);
protected:
    //
    // The list of Boxes.
    //
    List<Box> lbox;
    //
    // Returns a reference to the List<Box>.
    //
    List<Box>& listBox();
    //
    // Returns a constant reference to the List<Box>.
    //
    const List<Box>& listBox() const;

private:
    //
    // The IndexType of Boxes in the BoxList.
    //
    IndexType btype;
};

//
//@Man:
//@Memo: An Iterator for BoxList
/*@Doc:

  An iterator class allowing traversal over the Boxes in a BoxList.

  This class does support copying -- the compiler will generate
  a copy constructor on the fly if you need it.

  This class does NOT provide a default constructor or an assignment operator.
*/

class BoxListIterator
    : public ListIterator<Box>
{
public:
    //
    //@ManDoc: Construct a BoxListIterator on the supplied BoxList.
    //
    BoxListIterator (const BoxList& bl)
        : ListIterator<Box>(bl.lbox) {}
};

//
// Inlines.
//

inline
BoxList::~BoxList()
{}

inline
IndexType
BoxList::ixType () const
{
    return btype;
}

inline
void
BoxList::append (const Box& bn)
{
    assert(ixType() == bn.ixType());
    lbox.append(bn);
}

inline
void
BoxList::add (const Box& bn)
{
    append(bn);
}

inline
void
BoxList::prepend (const Box& bn)
{
    assert(ixType() == bn.ixType());
    lbox.prepend(bn);
}

inline
void
BoxList::join (const BoxList& blist)
{
    assert(ixType() == blist.ixType());
    lbox.join(blist.lbox);
}

inline
void
BoxList::clear ()
{
    lbox.clear();
}

inline
void
BoxList::addAfter (BoxListIterator& bli,
                   const Box&       bn)
{
    assert(ixType() == bn.ixType());
    lbox.addAfter(bli, bn);
}

inline
int
BoxList::length () const
{
    return lbox.length();
}

inline
bool
BoxList::isEmpty () const
{
    return lbox.isEmpty();
}

inline
bool
BoxList::isNotEmpty () const
{
    return lbox.isNotEmpty();
}

inline
bool
BoxList::contains (const Box& b) const
{
    BoxList bnew = ::complementIn(b,*this);
    return bnew.isEmpty();
}

inline
BoxList&
BoxList::remove (const Box& bx)
{
    assert(ixType() == bx.ixType());
    lbox.remove(bx);
    return *this;
}

inline
void
BoxList::addBefore (BoxListIterator& bli,
                    const Box&       bn)
{
    assert(ixType() == bn.ixType());
    lbox.addBefore(bli, bn);
}

inline
const Box&
BoxList::operator[] (const BoxListIterator& bli) const
{
    return lbox[bli];
}

inline
Box&
BoxList::operator[] (const BoxListIterator& bli)
{
    return lbox[bli];
}

inline
BoxList&
BoxList::remove (BoxListIterator& bli)
{
    assert(ixType() == bli().ixType());
    lbox.remove(bli);
    return *this;
}

inline
BoxList
intersect (const BoxList& bl,
           const Box&     b)
{
    BoxList newbl(bl);
    return newbl.intersect(b);
}

inline
BoxList
refine (const BoxList& bl,
        int            ratio)
{
    BoxList nbl(bl);
    return nbl.refine(ratio);
}

inline
BoxList
coarsen (const BoxList& bl,
         int            ratio)
{
    BoxList nbl(bl);
    return nbl.coarsen(ratio);
}

inline
BoxList
accrete (const BoxList& bl,
         int            sz)
{
    BoxList nbl(bl);
    return nbl.accrete(sz);
}

inline
List<Box>&
BoxList::listBox()
{
    return lbox;
}

inline
const List<Box>&
BoxList::listBox() const
{
    return lbox;
}

inline
bool
BoxList::operator!= (const BoxList& rhs) const
{
    return !operator==(rhs);
}

#endif /*BL_BOXLIST_H*/
