
#ifndef BL_BOXARRAY_H
#define BL_BOXARRAY_H

#include <iostream>
#include <cstddef>
#include <map>
#include <unordered_map>

#include <AMReX_IndexType.H>
#include <AMReX_BoxList.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

namespace amrex
{
    class BoxArray;

    //! Make a BoxArray from the the complement of b2 in b1in.
    BoxArray boxComplement (const Box& b1in, const Box& b2);

    //! Make a BoxArray from the complement of BoxArray ba in Box b.
    BoxArray complementIn (const Box& b, const BoxArray& ba);

    //! Make a BoxArray from the intersection of Box b and BoxArray(+ghostcells).
    BoxArray intersect (const BoxArray& ba, const Box& b, int ng = 0);

    BoxArray intersect (const BoxArray& ba, const Box& b, const IntVect& ng);

    //! Make a BoxArray from the intersection of two BoxArrays.
    BoxArray intersect (const BoxArray& lhs, const BoxArray& rhs);

    //! Make a BoxList from the intersection of BoxArray and BoxList.
    BoxList intersect (const BoxArray& ba, const BoxList& bl);

    BoxArray convert (const BoxArray& ba, IndexType typ);
    BoxArray convert (const BoxArray& ba, const IntVect& typ);

    BoxArray coarsen (const BoxArray& ba, int ratio);
    BoxArray coarsen (const BoxArray& ba, const IntVect& ratio);

    //! Find the ghost cells of a given BoxArray.
    BoxList GetBndryCells (const BoxArray& ba, int ngrow);

    //! Read a BoxArray from a stream.  If b is true, read in a special way
    void readBoxArray (BoxArray& ba, std::istream& s, bool b);

    //! Note that two BoxArrays that match are not necessarily equal.
    bool match (const BoxArray& x, const BoxArray& y);

// \cond CODEGEN
struct BARef
{
    BARef ();
    explicit BARef (size_t size);
    explicit BARef (const Box& b);
    explicit BARef (const BoxList& bl);
    explicit BARef (BoxList&& bl) noexcept;
    explicit BARef (std::istream& is);
    BARef (const BARef& rhs);
    BARef& operator= (const BARef& rhs) = delete;

    ~BARef ();

    void define (const Box& bx);
    void define (const BoxList& bl);
    void define (BoxList&& bl) noexcept;
    void define (std::istream& is, int& ndims);
    //!
    void resize (long n);
#ifdef BL_MEM_PROFILING
    void updateMemoryUsage_box (int s);
    void updateMemoryUsage_hash (int s);
#endif

    inline bool HasHashMap () const {
        bool r;
#ifdef _OPENMP
#pragma omp atomic read
#endif
        r = has_hashmap;
        return r;
    }

    //
    //! The data.
    Vector<Box> m_abox;
    //
    //! Box hash stuff.
    mutable Box bbox;

    mutable IntVect crsn;

    typedef std::unordered_map< IntVect, std::vector<int>, IntVect::shift_hasher > HashType;
    //using HashType = std::map< IntVect,std::vector<int> >;

    mutable HashType hash;

    mutable bool has_hashmap = false;

    static int  numboxarrays;
    static int  numboxarrays_hwm;
    static long total_box_bytes;
    static long total_box_bytes_hwm;
    static long total_hash_bytes;
    static long total_hash_bytes_hwm;

    static void Initialize ();
    static void Finalize ();
    static bool initialized;
};

struct BATransformer
{
    virtual ~BATransformer () {}

    virtual BATransformer* clone () const = 0;
    virtual bool equal (const BATransformer& rhs) const = 0;

    virtual IndexType ixType () const = 0;
    virtual void setIxType (IndexType typ) = 0;

    virtual IntVect crseRatio () const = 0;
    virtual void setCrseRatio (const IntVect& crse_ratio) = 0;

    virtual IntVect doiLo () const = 0; //!< doi: Domain of Influence
    virtual IntVect doiHi () const = 0; //!< doi: Domain of Influence

    virtual bool simple () const = 0; //!< Is this simply a type transformer

    virtual Box operator() (const Box& bx) const = 0;
};

template <class T>
struct BATBase
    : public BATransformer
{
    explicit BATBase (IndexType typ = IndexType(),
                      const IntVect& crse_ratio = IntVect::TheUnitVector())
	: m_typ(typ), m_crse_ratio(crse_ratio) {}
    virtual ~BATBase () {}

    virtual BATransformer* clone () const final {
	return new T(dynamic_cast<const T&>(*this));
    }

    virtual bool equal (const BATransformer& rhs) const final {
	const T* p = dynamic_cast<const T*>(&rhs);
	return p && (*p) == (dynamic_cast<const T&>(*this));
    }

    virtual IndexType ixType () const final {
	return m_typ;
    }

    virtual void setIxType (IndexType typ) final {
	m_typ = typ;
    }

    virtual IntVect crseRatio () const final {
        return m_crse_ratio;
    }

    virtual void setCrseRatio (const IntVect& crse_ratio) final {
        m_crse_ratio = crse_ratio;
    }

    virtual bool simple () const override { return false; }

protected:
    IndexType m_typ;
    IntVect   m_crse_ratio;
};

struct DefaultBATransformer final
    : public BATBase<DefaultBATransformer>
{
    explicit DefaultBATransformer (IndexType typ = IndexType(),
                                const IntVect& crse_ratio = IntVect::TheUnitVector())
	: BATBase<DefaultBATransformer>(typ, crse_ratio) {}

    virtual ~DefaultBATransformer () final {}

    virtual IntVect doiLo () const final
        { return IntVect::TheZeroVector(); }

    virtual IntVect doiHi () const final
        { return m_typ.ixType(); }

    virtual bool simple () const final { return true; }

    virtual Box operator() (const Box& bx) const final
        { return amrex::convert(amrex::coarsen(bx,m_crse_ratio), m_typ); }

    bool operator== (const DefaultBATransformer& rhs) const
	{ return m_typ == rhs.m_typ and m_crse_ratio == rhs.m_crse_ratio; }
};
// \endcond

/**
* \brief A collection of Boxes stored in an Array.  It is a
* reference-counted concrete class, not a polymorphic one; i.e. you
* cannot use any of the List member functions with a BoxList.
*/

class MFIter;

class BoxArray
{
public:
    //! Construct an empty BoxArray.
    BoxArray ();

    //! The copy constructor.
    BoxArray (const BoxArray& rhs);

    //! The assignment operator.
    BoxArray& operator= (const BoxArray& rhs);

    //! The move constructor.
    BoxArray (BoxArray&& rhs) noexcept = default;

    //! The move assignment operator.
    BoxArray& operator=(BoxArray&& rhs) noexcept = default;

    //! The destructor.
    ~BoxArray() noexcept = default;

    //! Make a boxarray out of a single box
    explicit BoxArray (const Box& bx);

    //! Construct a BoxArray of the specified size.
    explicit BoxArray (size_t size);

    //! Construct a BoxArray from an array of Boxes of size nbox.
    BoxArray (const Box* bxvec,
              int        nbox);

    //! Construct a BoxArray from a BoxList.
    explicit BoxArray (const BoxList& bl);
    explicit BoxArray (BoxList&& bl) noexcept;

    BoxArray (const BoxArray& rhs, const BATransformer& trans);

    
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
    void define (BoxList&& bl) noexcept;

    //! Remove all Boxes from the BoxArray.
    void clear ();

    //! Resize the BoxArray.  See Vector<T>::resize() for the gory details.
    void resize (long len);

    //! Return the number of boxes in the BoxArray.
    long size () const noexcept { return m_ref->m_abox.size(); }

    //! Return the number of boxes that can be held in the current allocated storage
    long capacity () const noexcept { return m_ref->m_abox.capacity(); }

    //! Return whether the BoxArray is empty
    bool empty () const noexcept { return m_ref->m_abox.empty(); }

    //! Returns the total number of cells contained in all boxes in the BoxArray.
    long numPts() const noexcept;

    //! Returns the total number of cells (in double type) contained in all boxes in the BoxArray.
    double d_numPts () const noexcept;
    /**
    * \brief Initialize the BoxArray from the supplied istream.
    * It is an error if the BoxArray has already been initialized.
    * Note that the BoxArray in the istream must have been written
    * using writeOn().
    */
    int readFrom (std::istream& is);

    //! Output this BoxArray to a checkpoint file.
    std::ostream& writeOn (std::ostream&) const;

    //! Are the BoxArrays equal?
    bool operator== (const BoxArray& rhs) const noexcept;

    //! Are the BoxArrays not equal?
    bool operator!= (const BoxArray& rhs) const noexcept;

    bool operator== (const Vector<Box>& bv) const noexcept;
    bool operator!= (const Vector<Box>& bv) const noexcept;

    //!  Are the BoxArrays equal after conversion to cell-centered
    bool CellEqual (const BoxArray& rhs) const noexcept;

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
    bool coarsenable(int refinement_ratio, int min_width=1) const;
    bool coarsenable(const IntVect& refinement_ratio, int min_width=1) const;

    //! Coarsen each Box in the BoxArray to the specified ratio.
    BoxArray& coarsen (const IntVect& refinement_ratio);

    //! Grow and then coarsen each Box in the BoxArray.
    BoxArray& growcoarsen (int n, const IntVect& refinement_ratio);
    BoxArray& growcoarsen (IntVect const& ngrow, const IntVect& refinement_ratio);

    //! Grow each Box in the BoxArray by the specified amount.
    BoxArray& grow (int n);

    //! Grow each Box in the BoxArray by the specified amount.
    BoxArray& grow (const IntVect& iv);
    /**
    * \brief Grow each Box in the BoxArray on the low and high ends
    * by n_cell cells in the idir direction.
    */
    BoxArray& grow (int idir, int n_cell);
    /**
    * \brief Grow each Box in the BoxArray on the low end
    * by n_cell cells in the idir direction.
    */
    BoxArray& growLo (int idir, int n_cell);
    /**
    * \brief Grow each Box in the BoxArray on the high end
    * by n_cell cells in the idir direction.
    */
    BoxArray& growHi (int idir, int n_cell);
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

    BoxArray& convert (const IntVect& typ);

    //! Apply function (*fp)(Box) to each Box in the BoxArray.
    BoxArray& convert (Box (*fp)(const Box&));

    //! Apply Box::shift(int,int) to each Box in the BoxArray.
    BoxArray& shift (int dir, int nzones);

    //! Apply Box::shift(const IntVect &iv) to each Box in the BoxArray.
    BoxArray& shift (const IntVect &iv);

    //! Set element i in this BoxArray to Box ibox.
    void set (int i, const Box& ibox);

    //! Return element index of this BoxArray.
    Box operator[] (int index) const noexcept {
        Box r = m_ref->m_abox[index];
        if (m_simple) {
            r.coarsen(m_crse_ratio).convert(m_typ);
        } else {
            r = (*m_transformer)(m_ref->m_abox[index]);
        }
        return r;
    }

    //! Return element index of this BoxArray.
    Box operator[] (const MFIter& mfi) const noexcept;

    //! Return element index of this BoxArray.
    Box get (int index) const noexcept { return operator[](index); }

    //! Return cell-centered box at element index of this BoxArray.
    Box getCellCenteredBox (int index) const noexcept {
        return amrex::coarsen(m_ref->m_abox[index],m_crse_ratio);
    }

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
    Box minimalBox (long& npts_avg_box) const;

    /**
    * \brief True if the Box intersects with this BoxArray(+ghostcells).
    * The Box must have the same IndexType as those in this BoxArray.
    */
    bool intersects (const Box& b, int ng = 0) const;

    bool intersects (const Box& b, const IntVect& ng) const;

    //! Return intersections of Box and BoxArray
    std::vector< std::pair<int,Box> > intersections (const Box& bx) const;

    //! Return intersections of Box and BoxArray(+ghostcells).
    std::vector< std::pair<int,Box> > intersections (const Box& bx, bool first_only, int ng) const;

    std::vector< std::pair<int,Box> > intersections (const Box& bx, bool first_only, const IntVect& ng) const;

    //! intersect Box and BoxArray, then store the result in isects
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects) const;

    //! intersect Box and BoxArray(+ghostcells), then store the result in isects
    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects,
			bool first_only, int ng) const;

    void intersections (const Box& bx, std::vector< std::pair<int,Box> >& isects,
			bool first_only, const IntVect& ng) const;

    //! Return box - boxarray
    BoxList complementIn (const Box& b) const;
    void complementIn (BoxList& bl, const Box& b) const;

    //! Clear out the internal hash table used by intersections.
    void clear_hash_bin () const;

    //! Change the BoxArray to one with no overlap and then simplify it (see the simplify function in BoxList).
    void removeOverlap (bool simplify=true);

    //! whether two BoxArrays share the same data
    static bool SameRefs (const BoxArray& lhs, const BoxArray& rhs) { return lhs.m_ref == rhs.m_ref; }

    struct RefID {
        RefID () noexcept : data(nullptr) {}
        explicit RefID (BARef* data_) noexcept : data(data_) {}
        bool operator<  (const RefID& rhs) const noexcept { return std::less<BARef*>()(data,rhs.data); }
        bool operator== (const RefID& rhs) const noexcept { return data == rhs.data; }
        bool operator!= (const RefID& rhs) const noexcept { return data != rhs.data; }
        friend std::ostream& operator<< (std::ostream& os, const RefID& id);
    private:
        BARef* data;
    };

    //! Return a unique ID of the reference
    RefID getRefID () const noexcept { return RefID { m_ref.get() }; }

    //! Return index type of this BoxArray
    IndexType ixType () const noexcept { return m_typ; }

    //! Return crse ratio of this BoxArray
    IntVect crseRatio () const noexcept { return m_crse_ratio; }

    static void Initialize ();
    static void Finalize ();
    static bool initialized;

    //! Make ourselves unique.
    void uniqify ();

private:
    //!  Update BoxArray index type according the box type, and then convert boxes to cell-centered.
    void type_update ();

    BARef::HashType& getHashMap () const;

    IntVect getDoiLo () const noexcept;
    IntVect getDoiHi () const noexcept;

    std::unique_ptr<BATransformer> m_transformer;
    IndexType m_typ;
    IntVect   m_crse_ratio;
    bool      m_simple;

    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<BARef> m_ref;
};

//! Write a BoxArray to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const BoxArray& ba);

std::ostream& operator<< (std::ostream& os, const BoxArray::RefID& id);

}

#endif /*BL_BOXARRAY_H*/
