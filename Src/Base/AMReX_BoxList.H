
#ifndef BL_BOXLIST_H
#define BL_BOXLIST_H

#include <iosfwd>

#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

namespace amrex
{
    class BoxArray;
    class BoxList;

    //! Returns a BoxList defining the complement of BoxList bl in Box b.
    BoxList complementIn (const Box& b, const BoxList& bl);

    //! Returns BoxList defining the compliment of b2 in b1in.
    BoxList boxDiff (const Box& b1in, const Box& b2);

    void boxDiff (BoxList& bl_diff, const Box& b1in, const Box& b2);

    //! Returns a new BoxList in which each Box is refined by the given ratio.
    BoxList refine (const BoxList&  bl, int ratio);

    //! Returns a new BoxList in which each Box is coarsened by the given ratio.
    BoxList coarsen (const BoxList& bl, int ratio);

    //! Returns a BoxList defining the intersection of bl with b.
    BoxList intersect (const BoxList& bl, const Box& b);

    //! Returns a new BoxList in which each Box is grown by the given size.
    BoxList accrete (const BoxList& bl, int sz);

    //! Return BoxList which covers the same area but has no overlapping boxes.
    BoxList removeOverlap (const BoxList& bl);


    //! Output a BoxList to an ostream in ASCII format.
    std::ostream& operator<< (std::ostream&  os, const BoxList& blist);


/** \brief A class for managing a List of Boxes that share a common
* IndexType.  This class implements operations for sets of Boxes.
* This is a concrete class, not a polymorphic one.
*/

class BoxList
{
public:

    friend class BoxDomain;

    typedef Vector<Box>::iterator       iterator;
    typedef Vector<Box>::const_iterator const_iterator;

    //! Construct an empty BoxList with IndexType::TheCellType().
    BoxList ();

    //! Make a boxList consisting of a single Box
    explicit BoxList (const Box& bx);

    //! Construct an empty BoxList with IndexType btyp
    explicit BoxList (IndexType btyp);

    //!  Construct a BoxList from a BoxArray.
    explicit BoxList (const BoxArray& ba);

    //! Construct a boxlist from a rvalue Vector<Box>.
    explicit BoxList (Vector<Box>&& bxs);

    //! Make a boxList consisting of non-overlapping tile box from a single Box
    BoxList (const Box& bx, const IntVect& tilesize);

    //! Make a BoxList consisting of nboxes Boxes from a single Box
    BoxList (const Box& bx, int nboxes);

    //! Make a BoxList consisting of nboxes Boxes chopped in dir-direction from a single Box
    BoxList (const Box& bx, int nboxes, Direction dir);

    BoxList (const BoxList& rhs) = default;
    BoxList (BoxList&& rhs) = default;
    BoxList& operator= (const BoxList& rhs) = default;
    BoxList& operator= (BoxList&& rhs) = default;

    void reserve (std::size_t n) { m_lbox.reserve(n); }

    //! Append a Box to this BoxList.
    void push_back (const Box& bn) { BL_ASSERT(ixType() == bn.ixType()); m_lbox.push_back(bn); }

    Box& front () noexcept { BL_ASSERT(!m_lbox.empty()); return m_lbox.front(); }

    const Box& front () const noexcept { BL_ASSERT(!m_lbox.empty()); return m_lbox.front(); }

    //! Join the BoxList to ourselves.
    void join (const BoxList& blist);
    //! Join the Array of Boxes to ourselves.
    void join (const Vector<Box>& barr);
    //! Catenate the BoxList to ourselves. Removes entries from blist.
    void catenate (BoxList& blist);

    //! Remove all Boxes from this BoxList.
    void clear ();
    //! The number of Boxes in this BoxList.
    long size () const noexcept { return m_lbox.size(); }

    std::size_t capacity () const noexcept { return m_lbox.capacity(); }

    iterator begin () noexcept { return m_lbox.begin(); }
    const_iterator begin () const noexcept { return m_lbox.begin(); }

    iterator end () noexcept { return m_lbox.end(); }
    const_iterator end () const noexcept { return m_lbox.end(); }
    /**
    * \brief True if this BoxList is valid; i.e. all the Boxes are
    * valid and they all have the same IndexType.  Also returns
    * true if the BoxList is empty.
    */
    bool ok () const noexcept;
    //! Is this BoxList equal to rhs?
    bool operator== (const BoxList& rhs) const;
    //! Is this BoxList notequal to rhs?
    bool operator!= (const BoxList& rhs) const;
    //! Is this BoxList empty?
    bool isEmpty () const noexcept { return m_lbox.empty(); }
    //! Is this BoxList not empty?
    bool isNotEmpty () const noexcept { return !m_lbox.empty(); }
    //! True if the set of intersecting Boxes is empty.
    bool isDisjoint () const;
    //! True if all Boxes in bl are contained in this BoxList.
    bool contains (const BoxList& bl) const;
    //! Modify this BoxList to contain only its intersection with Box b.
    BoxList& intersect (const Box& b);
    //! Modify this BoxList to contain only its intersection with BoxList bl.
    BoxList& intersect (const BoxList& bl);
    //! Remove empty Boxes from this BoxList.
    BoxList& removeEmpty();

    BoxList& complementIn (const Box&     b,
                           const BoxList& bl);
    BoxList& complementIn (const Box&     b,
                            BoxList&&     bl);
    BoxList& complementIn (const Box& b, const BoxArray& ba);
    //! Refine each Box in the BoxList by the ratio.
    BoxList& refine (int ratio);
    //! Refine each Box in the BoxList by the ratio.
    BoxList& refine (const IntVect& ratio);
    //! Coarsen each Box in the BoxList by the ratio.
    BoxList& coarsen (int ratio);
    //! Coarsen each Box in the BoxList by the ratio.
    BoxList& coarsen (const IntVect& ratio);
    //! Grow each Box in the BoxList by size sz.
    BoxList& accrete  (int sz);
    BoxList& accrete  (const IntVect& sz);
    //! Applies Box::shift(int,int) to each Box in the BoxList.
    BoxList& shift (int dir,
                    int nzones);
    //! Applies Box::shiftHalf(int,int) to each Box in the BoxList.
    BoxList& shiftHalf (int dir,
                        int num_halfs);
    //! Applies Box::shiftHalf(IntVect) to each Box in BoxList.
    BoxList& shiftHalf (const IntVect& iv);
    /**
    * \brief Merge adjacent Boxes in this BoxList. Return the number
    * of Boxes merged.  If "best" is specified we do a single
    * bruteforce pass over the list checking each Box against
    * all Boxes after it in the list to see if they can be
    * merged.  If "best" is not specified we limit how fair
    * afield we look for possible matches.  The "best" algorithm
    * is O(N-squared) while the other algorithm is O(N).
    */
    int simplify (bool best = false);
    //! Forces each Box in the BoxList to have sides of length <= chunk.
    BoxList& maxSize (int chunk);
    //! Forces each Box in the BoxList to have dimth side of length <= chunk[dim].
    BoxList& maxSize (const IntVect& chunk);
    //! Returns smallest Box that contains all Boxes in this BoxList.
    Box minimalBox () const;
    //! Returns the IndexType of Boxes in this BoxList.
    IndexType ixType () const noexcept { return btype; }
    //! Set the type of the BoxList.  It's an error if the BoxList isn't empty.
    void set (IndexType ixtyp) noexcept { BL_ASSERT(m_lbox.empty()); btype = ixtyp; }
    /**
    * \brief Applies surroundingNodes(Box) to each Box in BoxArray.
    * See the documentation of Box for details.
    */
    BoxList& surroundingNodes () noexcept;
    /**
    * \brief Applies surroundingNodes(Box,int) to each Box in
    * BoxList.  See the documentation of Box for details.
    */
    BoxList& surroundingNodes (int dir) noexcept;
    //! Applies Box::enclosedCells() to each Box in the BoxList.
    BoxList& enclosedCells () noexcept;
    //! Applies Box::enclosedCells(int) to each Box in the BoxList.
    BoxList& enclosedCells  (int dir) noexcept;
    //! Applies Box::convert(IndexType) to each Box in the BoxList.
    BoxList& convert (IndexType typ) noexcept;

    //! Returns a reference to the Vector<Box>.
    Vector<Box>& data() noexcept { return m_lbox; }
    //! Returns a constant reference to the Vector<Box>.
    const Vector<Box>& data() const noexcept { return m_lbox; }

    void swap (BoxList& rhs) {
	std::swap(m_lbox, rhs.m_lbox);
	std::swap(btype, rhs.btype);
    }

private:
    //! Core simplify routine.
    int simplify_doit (bool best);

    //! The list of Boxes.
    Vector<Box> m_lbox;
    //! The IndexType of Boxes in the BoxList.
    IndexType btype;

};

}

#endif /*BL_BOXLIST_H*/
