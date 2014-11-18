
#ifndef BL_BOXARRAY_H
#define BL_BOXARRAY_H

#include <cstddef>
#include <map>

#include <IndexType.H>
#include <BoxList.H>
#include <Array.H>
#include <Pointers.H>
#include <BaseFab.H>

class BoxArray;

namespace BoxLib
{
    //
    // Makes a BoxArray from the the complement of b2 in b1in.
    //
    BoxArray boxComplement (const Box& b1in, const Box& b2);
    //
    // Makes a BoxArray from the complement of BoxArray ba in Box b.
    //
    BoxArray complementIn (const Box& b, const BoxArray& ba);
    //
    // Makes a BoxArray from the intersection of Box b and BoxArray ba.
    //
    BoxArray intersect (const BoxArray& ba, const Box& b);
    //
    // Makes a BoxArray from the intersection of two BoxArrays.
    //
    BoxArray intersect (const BoxArray& lhs, const BoxArray& rhs);
    //
    // Find the ghost cells of a given BoxArray.
    //
    BoxList GetBndryCells (const BoxArray& ba, int ngrow);
}

//
// Write a BoxArray to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream& os, const BoxArray& ba);

//
// An Array of Boxes
//
// A BoxArray is a collection of Boxes stored in an Array.  It is a
// reference-counted concrete class, not a polymorphic one; i.e. you
// cannot use any of the List member functions with a BoxList.
//

class BoxArray
{
public:
    //
    // Construct an empty BoxArray.
    //
    BoxArray ();
    //
    /// Make a boxarray out of a single box
    //
    BoxArray (const Box& bx);
    //
    // Construct a BoxArray of the specified size.
    //
    explicit BoxArray (size_t size);
    //
    // Construct a BoxArray from an array of Boxes of size nbox.
    //
    BoxArray (const Box* bxvec,
              int        nbox);
    //
    // Construct a BoxArray from a BoxList.
    //
    explicit BoxArray (const BoxList& bl);
    //
    // The copy constructor.
    //
    BoxArray (const BoxArray& rhs);
    //
    // The assignment operator.
    //
    BoxArray& operator= (const BoxArray& rhs);
    //
    // The destructor.
    //
    ~BoxArray ();
    //
    // BoxArray can be used like an const STL container.
    // We ONLY supply const_iterators.  We do this because
    // BoxArrays are reference counted and we don't want folks
    // to modify without going through the non-STL BoxArray interface.
    //
    typedef Array<Box>::const_iterator const_iterator;

    const_iterator begin () const { return m_ref->m_abox.begin(); }
    const_iterator   end () const { return m_ref->m_abox.end();   }
    //
    // Initialize the BoxArray from the supplied BoxList.
    // It is an error if the BoxArray has already been initialized.
    //
    void define (const BoxList& bl);
    //
    // Initialize the BoxArray from the supplied BoxArray.
    // It is an error if the BoxArray has already been initialized.
    //
    void define (const BoxArray& bs);
    //
    // Remove all Boxes from the BoxArray.
    //
    void clear ();
    //
    // Resizes the array.  See Array<T>::resize() for the gory details.
    //
    void resize (int len);
    //
    // Returns the number of elements in the BoxArray.
    //
    int size () const;

    bool empty () const;

    long numPts() const;
    long long ll_numPts() const;
    //
    // Returns the number of points contained in the BoxArray.
    // This is inteded for use only in diagnostic messages.
    //
    double d_numPts () const;
    //
    // Initialize the BoxArray from the supplied istream.
    // It is an error if the BoxArray has already been initialized.
    // Note that the BoxArray in the istream must have been written
    // using writeOn().
    //
    void readFrom (std::istream& is);
    //
    // Output this BoxArray to a checkpoint file.
    //
    std::ostream& writeOn (std::ostream&) const;
    //
    // Are the BoxArrays equal?
    //
    bool operator== (const BoxArray& rhs) const;
    //
    // Are the BoxArrays not equal?
    //
    bool operator!= (const BoxArray& rhs) const;
    //
    // Forces each Box in BoxArray to have sides <= block_size.
    //
    BoxArray& maxSize (int block_size);

    BoxArray& maxSize (const IntVect& block_size);
    //
    // Refine each Box in the BoxArray to the specified ratio.
    //
    BoxArray& refine (int refinement_ratio);
    //
    // Refine each Box in the BoxArray to the specified ratio.
    //
    BoxArray& refine (const IntVect& refinement_ratio);
    //
    // Coarsen each Box in the BoxArray to the specified ratio.
    //
    BoxArray& coarsen (int refinement_ratio);
    //
    // Coarsen each Box in the BoxArray to the specified ratio.
    //
    BoxArray& coarsen (const IntVect& refinement_ratio);
    //
    // Grow each Box in the BoxArray by the specified amount.
    //
    BoxArray& grow (int n);
    //
    // Grow each Box in the BoxArray by the specified amount.
    //
    BoxArray& grow (const IntVect& iv);
    //
    // Grow each Box in the BoxArray on the low and high end
    // by n_cell cells in idir direction.
    //
    BoxArray& grow (int idir,
                    int n_cell);
    //
    // Applies surroundingNodes(Box) to each Box in BoxArray.
    // See the documentation of Box for details.
    //
    BoxArray& surroundingNodes ();
    //
    // Applies surroundingNodes(Box,int) to each Box in
    // BoxArray.  See the documentation of Box for details.
    //
    BoxArray& surroundingNodes (int dir);
    //
    // Applies Box::enclosedCells() to each Box in the BoxArray.
    //
    BoxArray& enclosedCells ();
    //
    // Applies Box::enclosedCells(int) to each Box in the BoxArray.
    //
    BoxArray& enclosedCells  (int dir);
    //
    // Applies Box::convert(IndexType) to each Box in the BoxArray.
    //
    BoxArray& convert (IndexType typ);
    //
    // Applies function (*fp)(Box) to each Box in the BoxArray.
    //
    BoxArray& convert (Box (*fp)(const Box&));
    //
    // Applies Box::shift(const IntVect &iv) to each Box in the BoxArray.
    //
    BoxArray& shift (const IntVect &iv);
    //
    // Applies Box::shift(int,int) to each Box in the BoxArray.
    //
    BoxArray& shift (int dir,
                     int nzones);
    //
    // Applies Box::shiftHalf(int,int) to each Box in the BoxArray.
    //
    BoxArray& shiftHalf (int dir,
                         int num_halfs);
    //
    // Applies Box::shiftHalf(IntVect) to each Box in BoxArray.
    //
    BoxArray& shiftHalf (const IntVect &iv);
    //
    // Set element i in this BoxArray to Box ibox.
    //
    void set (int        i,
              const Box& ibox);
    //
    // Returns element index of this BoxArray.
    //
    const Box& operator[] (int index) const { return m_ref->m_abox.get(index); }
    //
    // Returns element index of this BoxArray.
    //
    const Box& get (int index) const { return m_ref->m_abox.get(index); }
    //
    // Returns true if Box is valid and they all have the same
    // IndexType.  Is true by default if the BoxArray is empty.
    //
    bool ok () const;
    //
    // Returns true if set of intersecting Boxes in BoxArray is null.
    //
    bool isDisjoint () const;
    //
    // Create a BoxList from this BoxArray.
    //
    BoxList boxList () const;
    //
    // True if the Box intersects with this BoxArray.
    // The Box must have the same IndexType as those in this
    // BoxArray.
    //
    bool intersects (const Box& b) const;
    //
    // True if the IntVect is within any of the Boxes in this BoxArray.
    //
    bool contains (const IntVect& v) const;
    //
    // True if the Box is within any of the Boxes in the List.
    // The Box must also have the same IndexType as those in this
    // BoxArray.
    //
    bool contains (const Box& b) const;
    //
    // True if all Boxes in bl are contained in this BoxArray.
    //
    bool contains (const BoxArray& bl) const;
    //
    // Returns smallest Box that contains all Boxes in this BoxArray.
    //
    Box minimalBox () const;

    std::vector< std::pair<int,Box> > intersections (const Box& bx, bool first_only=false) const;

    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects, bool first_only=false) const;
    //
    // Clear out the internal hash table used by intersections.
    //
    void clear_hash_bin () const;
    //
    // Changes the BoxArray to one with no overlap.
    // Also tries to "simplify" the BoxArray within reason.
    //
    void removeOverlap ();

    static bool SameRefs (const BoxArray& lhs, const BoxArray& rhs) { return lhs.m_ref == rhs.m_ref; }

private:

    void decrementCounters () const;

    class Ref
    {
        friend class BoxArray;
        //
        // Constructors to match those in BoxArray ....
        //
        Ref () {}

        Ref (size_t size) : m_abox(size) {}

        Ref (const BoxList& bl) { define(bl); }

        Ref (std::istream& is) { define(is); }

        Ref (const Ref& rhs) : m_abox(rhs.m_abox) {}
        //
        // Some defines()s to match those in BoxArray.
        //
        void define (const BoxList& bl);

        void define (std::istream& is);
        //
        // The data.
        //
        Array<Box> m_abox;
        //
        // Box hash stuff.
        //
        mutable Box bbox;

        mutable IntVect crsn;

        mutable std::map< IntVect,std::vector<int>,IntVect::Compare > hash;
    };
    //
    // Make ourselves unique.
    //
    void uniqify ();
    //
    // The data -- a reference-counted pointer to a Ref.
    //
    LnClassPtr<Ref> m_ref;
};

#endif /*BL_BOXARRAY_H*/
