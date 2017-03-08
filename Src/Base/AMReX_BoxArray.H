
#ifndef BL_BOXARRAY_H
#define BL_BOXARRAY_H

#include <iostream>
#include <cstddef>
#include <map>
#include <unordered_map>

#include <AMReX_IndexType.H>
#include <AMReX_BoxList.H>
#include <AMReX_Array.H>
#include <AMReX_BaseFab.H>

namespace amrex
{
    class BoxArray;

    //! Make a BoxArray from the the complement of b2 in b1in.
    BoxArray boxComplement (const Box& b1in, const Box& b2);

    //! Make a BoxArray from the complement of BoxArray ba in Box b.
    BoxArray complementIn (const Box& b, const BoxArray& ba);

    //! Make a BoxArray from the intersection of Box b and BoxArray(+ghostcells).
    BoxArray intersect (const BoxArray& ba, const Box& b, int ng = 0);

    //! Make a BoxArray from the intersection of two BoxArrays.
    BoxArray intersect (const BoxArray& lhs, const BoxArray& rhs);

    BoxArray convert (const BoxArray& ba, IndexType typ);

    //! Find the ghost cells of a given BoxArray.
    BoxList GetBndryCells (const BoxArray& ba, int ngrow);

    //! Read a BoxArray from a stream.  If b is true, read in a special way
    void readBoxArray (BoxArray& ba, std::istream& s, bool b);

    //! delete entries in the coarseBoxArrayCache, which have key "key"
    void clearCoarseBoxArrayCache (ptrdiff_t key);

    //! Create a 1D Array of cells from all boxes in a BoxArray 
    Array<int> SerializeBoxArray(const BoxArray &ba);

    //! Construct a BoxArray (in BL_SPACEDIM dimensions) from a 1D Array of cells
    BoxArray UnSerializeBoxArray(const Array<int> &serarray);

    //! Note that two BoxArrays that match are not necessarily equal.
    bool match (const BoxArray& x, const BoxArray& y);

// \cond CODEGEN
struct BARef
{
    BARef ();
    explicit BARef (size_t size);
    explicit BARef (const Box& b);
    explicit BARef (const BoxList& bl);
    explicit BARef (std::istream& is);
    BARef (const BARef& rhs);
    BARef& operator= (const BARef& rhs) = delete;
    
    ~BARef ();

    //! This gives a unique ID of the reference
    ptrdiff_t getRefID () const;

    void define (const Box& bx);
    void define (const BoxList& bl);
    void define (std::istream& is);
    //!
    void resize (long n);
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
    
    typedef std::unordered_map< IntVect, std::vector<int>, IntVect::shift_hasher > HashType;
    //using HashType = std::map< IntVect,std::vector<int> >;

    mutable HashType hash;
    
    static int  numboxarrays;
    static int  numboxarrays_hwm;
    static long total_box_bytes;
    static long total_box_bytes_hwm;
    static long total_hash_bytes;
    static long total_hash_bytes_hwm;
        
    static void Initialize ();
    static bool initialized;
};

struct BATransformer
{
    virtual ~BATransformer () {}

    virtual BATransformer* clone () const = 0;
    virtual bool equal (const BATransformer& rhs) const = 0;

    virtual IndexType ixType () const = 0;
    virtual void setIxType (IndexType typ) = 0;
    
    virtual IntVect doiLo () const = 0; // doi: Domain of Influence
    virtual IntVect doiHi () const = 0; // doi: Domain of Influence

    virtual bool simple () const = 0; // Is this simply a type transformer

    virtual Box operator() (const Box& bx) const = 0;
};

template <class T>
struct BATBase
    : public BATransformer
{
    explicit BATBase (IndexType typ = IndexType())
	: m_typ(typ) {}
    virtual ~BATBase () {}

    virtual BATransformer* clone () const override {
	return new T(dynamic_cast<const T&>(*this));
    }

    virtual bool equal (const BATransformer& rhs) const override {
	const T* p = dynamic_cast<const T*>(&rhs);
	return p && (*p) == (dynamic_cast<const T&>(*this));
    }

    virtual IndexType ixType () const override {
	return m_typ;
    }

    virtual void setIxType (IndexType typ) override {
	m_typ = typ;
    }

    virtual bool simple () const override { return false; }

protected:
    IndexType m_typ;
};

struct BATypeTransformer
    : public BATBase<BATypeTransformer>
{
    explicit BATypeTransformer (IndexType typ = IndexType())
	: BATBase<BATypeTransformer>(typ) {}

    virtual ~BATypeTransformer () override {}

    virtual IntVect doiLo () const override 
        { return IntVect::TheZeroVector(); }

    virtual IntVect doiHi () const override 
        { return m_typ.ixType(); }

    virtual bool simple () const override { return true; }

    virtual Box operator() (const Box& bx) const override
        { return amrex::convert(bx, m_typ); }

    bool operator== (const BATypeTransformer& rhs) const
	{ return m_typ == rhs.m_typ; }
};
// \endcond

/**
* \brief A collection of Boxes stored in an Array.  It is a
* reference-counted concrete class, not a polymorphic one; i.e. you
* cannot use any of the List member functions with a BoxList.
*/

class BoxArray
{
public:
    //! Construct an empty BoxArray.
    BoxArray ();

    //! Make a boxarray out of a single box
    explicit BoxArray (const Box& bx);

    //! Construct a BoxArray of the specified size.
    explicit BoxArray (size_t size);

    //! Construct a BoxArray from an array of Boxes of size nbox.
    BoxArray (const Box* bxvec,
              int        nbox);

    //! Construct a BoxArray from a BoxList.
    explicit BoxArray (const BoxList& bl);

    BoxArray (const BoxArray& rhs, const BATransformer& trans);

    //! The copy constructor.
    BoxArray (const BoxArray& rhs);

    //! The assignment operator.
    BoxArray& operator= (const BoxArray& rhs);

    //! The move constructor
    BoxArray (BoxArray&& rhs) noexcept;
    
    //! The destructor.
    ~BoxArray ();
    /**
    * \brief Initialize the BoxArray from a single box.
    * It is an error if the BoxArray has already been initialized.
    */
    void define (const Box& bx);
    /**
    * \brief Initialize the BoxArray from the supplied BoxList.
    * It is an error if the BoxArray has already been initialized.
    */
    void define (const BoxList& bl);

    //! Remove all Boxes from the BoxArray.
    void clear ();

    //! Resize the BoxArray.  See Array<T>::resize() for the gory details.
    void resize (long len);

    //! Return the number of boxes in the BoxArray.
    long size () const;

    //! Return the number of boxes that can be held in the current allocated storage
    long capacity () const;

    //! Return whether the BoxArray is empty
    bool empty () const;

    //! Returns the total number of cells contained in all boxes in the BoxArray.
    long numPts() const;

    //! Returns the total number of cells (in double type) contained in all boxes in the BoxArray.
    double d_numPts () const;
    /**
    * \brief Initialize the BoxArray from the supplied istream.
    * It is an error if the BoxArray has already been initialized.
    * Note that the BoxArray in the istream must have been written
    * using writeOn().
    */
    void readFrom (std::istream& is);

    //! Output this BoxArray to a checkpoint file.
    std::ostream& writeOn (std::ostream&) const;

    //! Are the BoxArrays equal?
    bool operator== (const BoxArray& rhs) const;

    //! Are the BoxArrays not equal?
    bool operator!= (const BoxArray& rhs) const;

    //!  Are the BoxArrays equal after conversion to cell-centered
    bool CellEqual (const BoxArray& rhs) const;

    //! Forces each Box in BoxArray to have sides <= block_size.
    BoxArray& maxSize (int block_size);

    BoxArray& maxSize (const IntVect& block_size);

    //! Refine each Box in the BoxArray to the specified ratio.
    BoxArray& refine (int refinement_ratio);

    //! Refine each Box in the BoxArray to the specified ratio.
    BoxArray& refine (const IntVect& refinement_ratio);

    //! Coarsen each Box in the BoxArray to the specified ratio.
    BoxArray& coarsen (int refinement_ratio);

    //! Coarsen each Box in the BoxArray to the specified ratio.
    BoxArray& coarsen (const IntVect& refinement_ratio);

    //! Grow and then coarsen each Box in the BoxArray.
    BoxArray& growcoarsen (int n, const IntVect& refinement_ratio);

    //! Grow each Box in the BoxArray by the specified amount.
    BoxArray& grow (int n);

    //! Grow each Box in the BoxArray by the specified amount.
    BoxArray& grow (const IntVect& iv);
    /**
    * \brief Grow each Box in the BoxArray on the low and high end
    * by n_cell cells in idir direction.
    */
    BoxArray& grow (int idir,
                    int n_cell);
    /**
    * \brief Apply surroundingNodes(Box) to each Box in BoxArray.
    * See the documentation of Box for details.
    */
    BoxArray& surroundingNodes ();
    /**
    * \brief Apply surroundingNodes(Box,int) to each Box in
    * BoxArray.  See the documentation of Box for details.
    */
    BoxArray& surroundingNodes (int dir);

    //! Apply Box::enclosedCells() to each Box in the BoxArray.
    BoxArray& enclosedCells ();

    //! Apply Box::enclosedCells(int) to each Box in the BoxArray.
    BoxArray& enclosedCells  (int dir);

    //! Apply Box::convert(IndexType) to each Box in the BoxArray.
    BoxArray& convert (IndexType typ);

    //! Apply function (*fp)(Box) to each Box in the BoxArray.
    BoxArray& convert (Box (*fp)(const Box&));

    //! Apply Box::shift(int,int) to each Box in the BoxArray.
    BoxArray& shift (int dir,
                     int nzones);

    //! Apply Box::shift(const IntVect &iv) to each Box in the BoxArray.
    BoxArray& shift (const IntVect &iv);

    //! Apply Box::shiftHalf(int,int) to each Box in the BoxArray.
    BoxArray& shiftHalf (int dir,
                         int num_halfs);
    //! Apply Box::shiftHalf(IntVect) to each Box in BoxArray.
    BoxArray& shiftHalf (const IntVect &iv);

    //! Set element i in this BoxArray to Box ibox.
    void set (int        i,
              const Box& ibox);

    //! Return element index of this BoxArray.
    Box operator[] (int index) const {
        if (m_simple) {
            return amrex::convert(m_ref->m_abox[index], m_typ);
        } else {
            return (*m_transformer)(m_ref->m_abox[index]); 
        }
    }
    Box get (int index) const { return operator[](index); }

    //! Return cell-centered box at element index of this BoxArray.
    const Box& getCellCenteredBox (int index) const
	{ return m_ref->m_abox[index]; }    

    /**
    * \brief Return true if Box is valid and they all have the same
    * IndexType.  Is true by default if the BoxArray is empty.
    */
    bool ok () const;

    //! Return true if set of intersecting Boxes in BoxArray is null.
    bool isDisjoint () const;

    //! Create a BoxList from this BoxArray.
    BoxList boxList () const;

    //! True if the IntVect is within any of the Boxes in this BoxArray.
    bool contains (const IntVect& v) const;

    /**
    * \brief True if the Box is within any of the Boxes in the List.
    * The Box must also have the same IndexType as those in this BoxArray.
    */
    bool contains (const Box& b, bool assume_disjoint_ba = false) const;

    //! True if all Boxes in bl are contained in this BoxArray.
    bool contains (const BoxArray& bl, bool assume_disjoint_ba = false) const;

    //! Return smallest Box that contains all Boxes in this BoxArray.
    Box minimalBox () const;

    /**
    * \brief True if the Box intersects with this BoxArray(+ghostcells).
    * The Box must have the same IndexType as those in this BoxArray.
    */
    bool intersects (const Box& b, int ng = 0) const;

    //! Return intersections of Box and BoxArray
    std::vector< std::pair<int,Box> > intersections (const Box& bx) const;

    //! Return intersections of Box and BoxArray(+ghostcells).
    std::vector< std::pair<int,Box> > intersections (const Box& bx, bool first_only, int ng) const;

    //! intersect Box and BoxArray, then store the result in isects
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects) const; 

    //! intersect Box and BoxArray(+ghostcells), then store the result in isects
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects, 
			bool first_only, int ng) const;

    //! Return box - boxarray
    BoxList complement (const Box& b) const;

    //! Clear out the internal hash table used by intersections.
    void clear_hash_bin () const;

    //! Change the BoxArray to one with no overlap and then simplify it (see the simplify function in BoxList).
    void removeOverlap ();

    //! whether two BoxArrays share the same data
    static bool SameRefs (const BoxArray& lhs, const BoxArray& rhs) { return lhs.m_ref == rhs.m_ref; }

    //! Return a unique ID of the reference
    ptrdiff_t getRefID () const { return m_ref->getRefID(); }

    //! Return index type of this BoxArray
    IndexType ixType () const { return m_typ; }

#ifdef BL_USE_MPI
    //! Broadcast the BoxArray
    static void SendBoxArray(const BoxArray &ba, int whichSidecar);
    //! Receive the BoxArray from root
    static void RecvBoxArray(BoxArray &ba, int whichSidecar);
#endif

    using CBAMap   = std::map<IntVect, BoxArray*>;
    using CBACache = std::map<ptrdiff_t, CBAMap>;

    static CBACache m_CoarseBoxArrayCache;

    static void Initialize ();
    static bool initialized;

private:
    //!  Update BoxArray index type according the box type, and then convert boxes to cell-centered.
    void type_update ();

    BARef::HashType& getHashMap () const;

    //! Make ourselves unique.
    void uniqify ();

    IntVect getDoiLo () const
        { return m_simple ? IntVect::TheZeroVector() : m_transformer->doiLo(); }

    IntVect getDoiHi () const
        { return m_simple ?           m_typ.ixType() : m_transformer->doiHi(); }

    std::unique_ptr<BATransformer> m_transformer;
    IndexType m_typ;
    bool m_simple;

    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<BARef> m_ref;    
};

//! Write a BoxArray to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const BoxArray& ba);

}

#endif /*BL_BOXARRAY_H*/
