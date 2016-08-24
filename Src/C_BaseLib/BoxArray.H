
#ifndef BL_BOXARRAY_H
#define BL_BOXARRAY_H

#include <iostream>
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
    // Makes a BoxArray from the intersection of Box b and BoxArray(+ghostcells).
    //
    BoxArray intersect (const BoxArray& ba, const Box& b, int ng = 0);
    //
    // Makes a BoxArray from the intersection of two BoxArrays.
    //
    BoxArray intersect (const BoxArray& lhs, const BoxArray& rhs);
    //
    BoxArray convert (const BoxArray& ba, IndexType typ);
    //
    // Find the ghost cells of a given BoxArray.
    //
    BoxList GetBndryCells (const BoxArray& ba, int ngrow);
    //
    // Read a BoxArray from a stream.  If b is true,
    //   read in a special way
    //
    void readBoxArray (BoxArray& ba, std::istream& s, bool b);

    void clearCoarseBoxArrayCache (ptrdiff_t key);

    Array<int> SerializeBoxArray(const BoxArray &ba);
    BoxArray UnSerializeBoxArray(const Array<int> &serarray);
}

class BARef
{
public:
    BARef ();
    explicit BARef (size_t size);
    BARef (const BoxList& bl);
    BARef (std::istream& is);
    BARef (const BARef& rhs);

    ~BARef ();

    // This gives a unique ID of the reference
    ptrdiff_t getRefID () const;

    void define (const Box& bx);
    void define (const BoxList& bl);
    void define (std::istream& is);
    //
    void resize (long n);
    //
#ifdef BL_MEM_PROFILING
    void updateMemoryUsage_box (int s);
    void updateMemoryUsage_hash (int s);
#endif
    //
    // The data.
    //
    Array<Box> m_abox;
    //
    // Box hash stuff.
    //
    mutable Box bbox;
    
    mutable IntVect crsn;
    
    typedef std::map< IntVect,std::vector<int>,IntVect::Compare > HashType;

    mutable HashType hash;
    
    static int  numboxarrays;
    static int  numboxarrays_hwm;
    static long total_box_bytes;
    static long total_box_bytes_hwm;
    static long total_hash_bytes;
    static long total_hash_bytes_hwm;
    
private:
    BARef& operator= (const BARef&); // disallowed
    
#ifdef BL_MEM_PROFILING
    static void Initialize ();
    static bool initialized;
#endif
};


template<>
struct LnClassPtrCallBack<BARef> 
{
    // This will be called when a unique ref is about to be deleted.
    void operator() (BARef* ref) {
	BoxLib::clearCoarseBoxArrayCache(ref->getRefID());
    }
};


//
// Thess classes are used to transform from cell-centered boxes stored in BARef
// (and hidden from the users) to the public boxes BoxArray is supposed to hold.
//
struct BATransformer
{
    virtual ~BATransformer () {}

    virtual BATransformer* clone () const = 0;
    virtual bool equal (const BATransformer& rhs) const = 0;

    virtual IndexType ixType () const = 0;
    virtual void setIxType (IndexType typ) = 0;
    
    virtual IntVect doiLo () const = 0; // doi: Domain of Influence
    virtual IntVect doiHi () const = 0; // doi: Domain of Influence

    virtual Box operator() (const Box& bx) const = 0;
};

template <class T>
struct BATBase
    : public BATransformer
{
    explicit BATBase (IndexType typ = IndexType())
	: m_typ(typ) {}
    virtual ~BATBase () {}

    virtual BATransformer* clone () const BL_OVERRIDE {
	return new T(dynamic_cast<const T&>(*this));
    }

    virtual bool equal (const BATransformer& rhs) const BL_OVERRIDE {
	const T* p = dynamic_cast<const T*>(&rhs);
	return p && (*p) == (dynamic_cast<const T&>(*this));
    }

    virtual IndexType ixType () const BL_OVERRIDE {
	return m_typ;
    }

    virtual void setIxType (IndexType typ) BL_OVERRIDE {
	m_typ = typ;
    }

protected:
    IndexType m_typ;
};

struct BATypeTransformer
    : public BATBase<BATypeTransformer>
{
    explicit BATypeTransformer (IndexType typ = IndexType())
	: BATBase<BATypeTransformer>(typ) {}

    virtual ~BATypeTransformer () BL_OVERRIDE {}

    virtual IntVect doiLo () const BL_OVERRIDE 
        { return IntVect::TheZeroVector(); }

    virtual IntVect doiHi () const BL_OVERRIDE 
        { return m_typ.ixType(); }

    virtual Box operator() (const Box& bx) const BL_OVERRIDE
        { return BoxLib::convert(bx, m_typ); }

    bool operator== (const BATypeTransformer& rhs) const
	{ return m_typ == rhs.m_typ; }
};

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
    explicit BoxArray (const Box& bx);
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
    BoxArray (const BoxArray& rhs, const BATransformer& trans);
    //
    // The assignment operator.
    //
    BoxArray& operator= (const BoxArray& rhs);
    
    //
    // The destructor.
    //
    ~BoxArray ();
    //
    // Initialize the BoxArray from a single box.
    // It is an error if the BoxArray has already been initialized.
    //
    void define (const Box& bx);
    //
    // Initialize the BoxArray from the supplied BoxList.
    // It is an error if the BoxArray has already been initialized.
    //
    void define (const BoxList& bl);
    //
    // Remove all Boxes from the BoxArray.
    //
    void clear ();
    //
    // Resizes the array.  See Array<T>::resize() for the gory details.
    //
    void resize (long len);
    //
    // Returns the number of elements in the BoxArray.
    //
    long size () const;

    long capacity () const;

    bool empty () const;

    long numPts() const;
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
    //  Are the BoxArrays equal after conversion to cell-centered
    //
    bool CellEqual (const BoxArray& rhs) const;
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
    // Grow and then coarsen each Box in the BoxArray.
    //
    BoxArray& growcoarsen (int n, const IntVect& refinement_ratio);
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
    // Applies Box::shift(int,int) to each Box in the BoxArray.
    //
    BoxArray& shift (int dir,
                     int nzones);
    //
    // Applies Box::shift(const IntVect &iv) to each Box in the BoxArray.
    //
    BoxArray& shift (const IntVect &iv);
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
    Box operator[] (int index) const 
	{ return (*m_transformer)(m_ref->m_abox[index]); }
    Box get        (int index) const 
	{ return (*m_transformer)(m_ref->m_abox[index]); }
    //
    // Returns cell-centered box at element index of this BoxArray.
    //
    const Box& getCellCenteredBox (int index) const
	{ return m_ref->m_abox[index]; }    
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
    //
    // True if the Box intersects with this BoxArray(+ghostcells).
    // The Box must have the same IndexType as those in this
    // BoxArray.
    //
    bool intersects (const Box& b, int ng = 0) const;
    //
    // Returns intersections of Box and BoxArray(+ghsotcells).
    //
    std::vector< std::pair<int,Box> > intersections (const Box& bx) const;
    std::vector< std::pair<int,Box> > intersections (const Box& bx, bool first_only, int ng) const;
    //
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects) const; 
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects, 
			bool first_only, int ng) const;
    // Return box - boxarray
    BoxList complement (const Box& b) const;
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

    //
    // This gives a unique ID of the reference
    //
    ptrdiff_t getRefID () const { return m_ref->getRefID(); }

    //
    // Returns index type of this BoxArray
    //
    IndexType ixType () const { return m_transformer->ixType(); }

#ifdef BL_USE_MPI
    static void SendBoxArray(const BoxArray &ba, int whichSidecar);
    static void RecvBoxArray(BoxArray &ba, int whichSidecar);
#endif

    typedef std::map<IntVect, BoxArray*, IntVect::Compare> CBAMap;
    typedef std::map<ptrdiff_t, CBAMap> CBACache;

    static CBACache m_CoarseBoxArrayCache;

private:
    //
    //  Update BoxArray index type according the box type, and then convert boxes to cell-centered.
    //
    void type_update ();

    BARef::HashType& getHashMap () const;

    //
    // Make ourselves unique.
    //
    void uniqify ();

    // When we move to C++11, we will switch to unique pointer
    BATransformer* m_transformer;

    //
    // The data -- a reference-counted pointer to a Ref.
    //
    LnClassPtr<BARef> m_ref;    
};

//
// Write a BoxArray to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream& os, const BoxArray& ba);

#endif /*BL_BOXARRAY_H*/
