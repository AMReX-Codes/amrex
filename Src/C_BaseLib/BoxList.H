#ifndef BL_BOXLIST_H
#define BL_BOXLIST_H

//
// $Id: BoxList.H,v 1.24 2001-07-24 18:16:52 lijewski Exp $
//

#include <iostream>
#include <list>

#include <IntVect.H>
#include <IndexType.H>
#include <Box.H>
//
// Forward declarations.
//
class BoxArray;
class BoxList;

namespace BoxLib
{
    /*@ManDoc: Returns a BoxList defining the complement of BoxList bl
      in Box b.
    */
    BoxList complementIn (const Box& b, const BoxList& bl);

    //
    //@ManDoc: Returns BoxList defining the compliment of b2 in b1in.
    //
    BoxList boxDiff (const Box& b1in, const Box& b2);

    /*@ManDoc: Returns a new BoxList in which each Box is refined by the
      given ratio.
    */
    BoxList refine (const BoxList&  bl, int ratio);

    /*@ManDoc: Returns a new BoxList in which each Box is coarsened by
      the given ratio.
    */
    BoxList coarsen (const BoxList& bl, int ratio);
    //
    //@ManDoc: Returns a BoxList defining the intersection of bl with b.
    //
    BoxList intersect (const BoxList& bl, const Box& b);
    //
    //@ManDoc: Returns a BoxList defining the intersection of bl with br.
    //
    BoxList intersect (const BoxList& bl, const BoxList& br);

    /*@ManDoc: Returns a new BoxList in which each Box is grown by the
      given size.
    */
    BoxList accrete (const BoxList& bl, int sz);
}

//
//@ManDoc: Output a BoxList to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream&  os, const BoxList& blist);

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
public:

    friend class BoxDomain;

    typedef std::list<Box>::iterator       Iterator;
    typedef std::list<Box>::const_iterator ConstIterator;
    //
    //@ManDoc: Construct an empty BoxList with IndexType::TheCellType().
    //
    BoxList ();
    //
    //@ManDoc: Construct an empty BoxList with IndexType btyp
    //
    explicit BoxList (IndexType btyp);
    //
    //@ManDoc:  Construct a BoxList from a BoxArray.
    //
    explicit BoxList (const BoxArray& ba);
    //
    //@ManDoc: Append a Box to this BoxList.
    //
    void append (const Box& bn);
    //
    //@ManDoc:  Prepend a Box to this BoxList.
    //
    void prepend (const Box& bn);
    //
    //@ManDoc: Join the BoxList to ourselves.
    //
    void join (const BoxList& blist);
    //
    //@ManDoc: Catenate the BoxList to ourselves. Removes entries from blist.
    //
    void catenate (BoxList& blist);
    //
    //@ManDoc: Remove all Boxes from this BoxList.
    //
    void clear ();
    //
    //@ManDoc: The number of Boxes in this BoxList.
    //
    size_t size () const;

    Iterator begin ();
    ConstIterator begin () const;

    Iterator end ();
    ConstIterator end () const;

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
    //
    //@ManDoc: True if all Boxes in ba are contained in this BoxList.
    //
    bool contains (const BoxArray& ba) const;

    /*@ManDoc: Modify this BoxList to contain only its intersection with
               Box b.
    */
    BoxList& intersect (const Box& b);

    /*@ManDoc: Modify this BoxList to contain only its intersection with
               BoxList b.
    */
    BoxList& intersect (const BoxList& b);
    //
    //@ManDoc: Remove the Box from this BoxList.
    //
    BoxList& remove (const Box& bx);
    //
    //@ManDoc: Remove the pointed to Box from this BoxList.
    //
    BoxList& remove (Iterator bli);
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
    //
    //@ManDoc: Forces each Box in the BoxList to have sides of length <= chunk.
    //
    BoxList& maxSize (int chunk);

    /*@ManDoc: Forces each Box in the BoxList to have dimth side of 
               length <= chunk[dim].
    */
    BoxList& maxSize (const IntVect& chunk);
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

protected:
    //
    // The list of Boxes.
    //
    std::list<Box> lbox;
    //
    // Returns a reference to the List<Box>.
    //
    std::list<Box>& listBox();
    //
    // Returns a constant reference to the List<Box>.
    //
    const std::list<Box>& listBox() const;

private:
    //
    // The IndexType of Boxes in the BoxList.
    //
    IndexType btype;
};

#endif /*BL_BOXLIST_H*/
