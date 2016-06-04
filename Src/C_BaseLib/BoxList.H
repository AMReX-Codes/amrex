
#ifndef BL_BOXLIST_H
#define BL_BOXLIST_H

#include <iosfwd>
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
    //
    // Returns a BoxList defining the complement of BoxList bl in Box b.
    //
    BoxList complementIn (const Box& b, const BoxList& bl);
    //
    // Returns BoxList defining the compliment of b2 in b1in.
    //
    BoxList boxDiff (const Box& b1in, const Box& b2);
    //
    // Returns a new BoxList in which each Box is refined by the given ratio.
    //
    BoxList refine (const BoxList&  bl, int ratio);
    //
    // Returns a new BoxList in which each Box is coarsened by the given ratio.
    //
    BoxList coarsen (const BoxList& bl, int ratio);
    //
    // Returns a BoxList defining the intersection of bl with b.
    //
    BoxList intersect (const BoxList& bl, const Box& b);
    //
    // Returns a BoxList defining the intersection of bl with br.
    //
    BoxList intersect (const BoxList& bl, const BoxList& br);
    //
    // Returns a new BoxList in which each Box is grown by the given size.
    //
    BoxList accrete (const BoxList& bl, int sz);
    //
    // Return BoxList which covers the same area but has no overlapping boxes.
    //
    BoxList removeOverlap (const BoxList& bl);
}

//
// Output a BoxList to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream&  os, const BoxList& blist);

//
// A List of Boxes
//
// A BoxList is a class for managing a List of Boxes that share a common
// IndexType.  This class implements operations for sets of Boxes.
//
// This is a concrete class, not a polymorphic one.
//

class BoxList
{
public:

    friend class BoxDomain;

    typedef std::list<Box>::iterator       iterator;
    typedef std::list<Box>::const_iterator const_iterator;
    //
    // Construct an empty BoxList with IndexType::TheCellType().
    //
    BoxList ();
    //
    /// Make a boxList consisting of a single Box
    //
    BoxList (const Box& bx);
    //
    // Construct an empty BoxList with IndexType btyp
    //
    explicit BoxList (IndexType btyp);
    //
    //  Construct a BoxList from a BoxArray.
    //
    explicit BoxList (const BoxArray& ba);
    //
    // Make a boxList consisting of non-overlapping tile box from a single Box
    //
    BoxList (const Box& bx, const IntVect& tilesize);
    //
    // Append a Box to this BoxList.
    //
    void push_back (const Box& bn) { BL_ASSERT(ixType() == bn.ixType()); lbox.push_back(bn); }
    //
    //  Prepend a Box to this BoxList.
    //
    void push_front (const Box& bn) { BL_ASSERT(ixType() == bn.ixType()); lbox.push_front(bn); }

    Box& front () { BL_ASSERT(!lbox.empty()); return lbox.front(); }

    const Box& front () const { BL_ASSERT(!lbox.empty()); return lbox.front(); }

    void pop_front () { BL_ASSERT(!lbox.empty()); lbox.pop_front(); }
    //
    // Join the BoxList to ourselves.
    //
    void join (const BoxList& blist);
    //
    // Catenate the BoxList to ourselves. Removes entries from blist.
    //
    void catenate (BoxList& blist);
    // 
    // Splic the BoxList to the front
    //
    void splice_front (BoxList& blist);
    //
    // Remove all Boxes from this BoxList.
    //
    void clear ();
    //
    // The number of Boxes in this BoxList.
    //
    int size () const { return lbox.size(); }

    iterator begin () { return lbox.begin(); }
    const_iterator begin () const { return lbox.begin(); }

    iterator end () { return lbox.end(); }
    const_iterator end () const { return lbox.end(); }
    //
    // True if this BoxList is valid; i.e. all the Boxes are
    // valid and they all have the same IndexType.  Also returns
    // true if the BoxList is empty.
    //
    bool ok () const;
    //
    // Is this BoxList equal to rhs?
    //
    bool operator== (const BoxList& rhs) const;
    //
    // Is this BoxList notequal to rhs?
    //
    bool operator!= (const BoxList& rhs) const;
    //
    // Is this BoxList empty?
    //
    bool isEmpty () const { return lbox.empty(); }
    //
    // Is this BoxList not empty?
    //
    bool isNotEmpty () const { return !lbox.empty(); }
    //
    // True if the set of intersecting Boxes is empty.
    //
    bool isDisjoint () const;
    //
    // True if the IntVect is within any of the Boxes in this BoxList.
    //
    bool contains (const IntVect& v) const;

    /* True if the Box is within any of the Boxes in the List.
               The Box must also have the same IndexType as those in this
               BoxList.
    */
    bool contains (const Box& b) const;
    //
    // True if all Boxes in bl are contained in this BoxList.
    //
    bool contains (const BoxList& bl) const;
    //
    // True if all Boxes in ba are contained in this BoxList.
    //
    bool contains (const BoxArray& ba) const;
    //
    // Modify this BoxList to contain only its intersection with Box b.
    //
    BoxList& intersect (const Box& b);
    //
    // Modify this BoxList to contain only its intersection with BoxList b.
    //
    BoxList& intersect (const BoxList& b);
    //
    // Remove the Box from this BoxList.
    //
    BoxList& remove (const Box& bx);
    //
    // Remove the pointed to Box from this BoxList.
    //
    BoxList& remove (iterator bli);
    //
    // Creates the complement of BoxList bl in Box b.
    //
    BoxList& complementIn (const Box&     b,
                           const BoxList& bl);
    //
    // Refine each Box in the BoxList by the ratio.
    //
    BoxList& refine (int ratio);
    //
    // Refine each Box in the BoxList by the ratio.
    //
    BoxList& refine (const IntVect& ratio);
    //
    // Coarsen each Box in the BoxList by the ratio.
    //
    BoxList& coarsen (int ratio);
    //
    // Coarsen each Box in the BoxList by the ratio.
    //
    BoxList& coarsen (const IntVect& ratio);
    //
    // Grow each Box in the BoxList by size sz.
    //
    BoxList& accrete  (int sz);

    BoxList& accrete  (IntVect sz);
    //
    // Applies Box::shift(int,int) to each Box in the BoxList.
    //
    BoxList& shift (int dir,
                    int nzones);
    //
    // Applies Box::shiftHalf(int,int) to each Box in the BoxList.
    //
    BoxList& shiftHalf (int dir,
                        int num_halfs);
    //
    // Applies Box::shiftHalf(IntVect) to each Box in BoxList.
    //
    BoxList& shiftHalf (const IntVect& iv);
    //
    // Merge adjacent Boxes in this BoxList. Return the number
    // of Boxes merged.  If "best" is specified we do a single
    // bruteforce pass over the list checking each Box against
    // all Boxes after it in the list to see if they can be
    // merged.  If "best" is not specified we limit how fair
    // afield we look for possible matches.  The "best" algorithm
    // is O(N-squared) while the other algorithm is O(N).
    //
    int simplify (bool best = false);
    //
    // Merges all possible adjacent Boxes.  Returns the total number of Boxes merged.
    //
    int minimize ();
    //
    // Forces each Box in the BoxList to have sides of length <= chunk.
    //
    BoxList& maxSize (int chunk);
    //
    // Forces each Box in the BoxList to have dimth side of length <= chunk[dim].
    //
    BoxList& maxSize (const IntVect& chunk);
    //
    // Returns smallest Box that contains all Boxes in this BoxList.
    //
    Box minimalBox () const;
    //
    // Returns the IndexType of Boxes in this BoxList.
    //
    IndexType ixType () const { return btype; }
    //
    // Set the type of the BoxList.  It's an error if the BoxList isn't empty.
    //
    void set (IndexType ixtyp) { BL_ASSERT(lbox.empty()); btype = ixtyp; }
    //
    // Applies surroundingNodes(Box) to each Box in BoxArray.
    // See the documentation of Box for details.
    //
    BoxList& surroundingNodes ();
    //
    // Applies surroundingNodes(Box,int) to each Box in
    // BoxList.  See the documentation of Box for details.
    //
    BoxList& surroundingNodes (int dir);
    //
    // Applies Box::enclosedCells() to each Box in the BoxList.
    //
    BoxList& enclosedCells ();
    //
    // Applies Box::enclosedCells(int) to each Box in the BoxList.
    //
    BoxList& enclosedCells  (int dir);
    //
    // Applies Box::convert(IndexType) to each Box in the BoxList.
    //
    BoxList& convert (IndexType typ);

private:
    //
    // This is a helper function for the real complementIn().
    //
    BoxList& complementIn_base (const Box&     b,
                                const BoxList& bl);
    //
    // Core simplify routine.
    //
    int simplify_doit (bool best);
    //
    // The list of Boxes.
    //
    std::list<Box> lbox;
    //
    // Returns a reference to the List<Box>.
    //
    std::list<Box>& listBox() { return lbox; }
    //
    // Returns a constant reference to the List<Box>.
    //
    const std::list<Box>& listBox() const { return lbox; }
    //
    // The IndexType of Boxes in the BoxList.
    //
    IndexType btype;

};

#endif /*BL_BOXLIST_H*/
