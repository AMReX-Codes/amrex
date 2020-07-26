
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
    void readBoxArray (BoxArray& ba, std::istream& s, bool b = false);

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
    void resize (Long n);
#ifdef AMREX_MEM_PROFILING
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
    static Long total_box_bytes;
    static Long total_box_bytes_hwm;
    static Long total_hash_bytes;
    static Long total_hash_bytes_hwm;

    static void Initialize ();
    static void Finalize ();
    static bool initialized;
};

struct BATnull
{
    constexpr Box operator() (const Box& bx) const noexcept { return bx; }
    constexpr Box coarsen (Box const& a_box) const noexcept { return a_box; }
    constexpr IntVect doiLo () const noexcept { return IntVect::TheZeroVector(); }
    constexpr IntVect doiHi () const noexcept { return IntVect::TheZeroVector(); }
    constexpr IndexType index_type () const noexcept { return IndexType(); }
    constexpr IntVect coarsen_ratio () const noexcept { return IntVect::TheUnitVector(); }
};

struct BATindexType
{
    explicit BATindexType (IndexType a_typ) : m_typ(a_typ) {}
    Box operator() (const Box& bx) const noexcept { return amrex::convert(bx,m_typ); }
    Box coarsen (Box const& a_box) const noexcept { return a_box; }
    IntVect doiLo () const noexcept { return IntVect::TheZeroVector(); }
    IntVect doiHi () const noexcept { return m_typ.ixType(); }
    IndexType index_type () const noexcept { return m_typ; }
    IntVect coarsen_ratio () const noexcept { return IntVect::TheUnitVector(); }
    IndexType m_typ;
};

struct BATcoarsenRatio
{
    explicit BATcoarsenRatio (IntVect const& a_crse_ratio) : m_crse_ratio(a_crse_ratio) {}
    Box operator() (const Box& bx) const noexcept { return amrex::coarsen(bx,m_crse_ratio); }
    Box coarsen (Box const& a_box) const noexcept { return amrex::coarsen(a_box,m_crse_ratio); }
    IntVect doiLo () const noexcept { return IntVect::TheZeroVector(); }
    IntVect doiHi () const noexcept { return IntVect::TheZeroVector(); }
    IndexType index_type () const noexcept { return IndexType(); }
    IntVect coarsen_ratio () const noexcept { return m_crse_ratio; }
    IntVect m_crse_ratio;
};

struct BATindexType_coarsenRatio
{
    BATindexType_coarsenRatio (IndexType a_typ, IntVect const& a_crse_ratio)
        : m_typ(a_typ), m_crse_ratio(a_crse_ratio) {}

    Box operator() (const Box& bx) const noexcept {
        return amrex::convert(amrex::coarsen(bx,m_crse_ratio),m_typ);
    }

    Box coarsen (Box const& a_box) const noexcept { return amrex::coarsen(a_box,m_crse_ratio); }

    IntVect doiLo () const noexcept { return IntVect::TheZeroVector(); }
    IntVect doiHi () const noexcept { return m_typ.ixType(); }

    IndexType index_type () const noexcept { return m_typ; }
    IntVect coarsen_ratio () const noexcept { return m_crse_ratio; }

    IndexType m_typ;
    IntVect m_crse_ratio;
};

struct BATbndryReg
{
    BATbndryReg (Orientation a_face, IndexType a_typ,
                 int a_in_rad, int a_out_rad, int a_extent_rad)
        : m_face(a_face), m_typ(a_typ), m_crse_ratio(1)
    {
        m_loshft = IntVect(-a_extent_rad);
        m_hishft = IntVect( a_extent_rad);
        IntVect nodal = a_typ.ixType();
        m_hishft += nodal;
        const int d = a_face.coordDir();
        if (m_face.isLow()) {
            m_loshft[d] = nodal[d] - a_out_rad;
            m_hishft[d] = nodal[d] + a_in_rad - 1;
        } else {
            m_loshft[d] = 1 - a_in_rad;
            m_hishft[d] = a_out_rad;
        }
        m_doilo = IntVect(a_extent_rad);
        m_doihi = IntVect(a_extent_rad);
        m_doihi += nodal;
        if (m_face.isLow()) {  // domain of influence in index space
            m_doilo[d] = std::max(0, a_out_rad - nodal[d]);
            m_doihi[d] = 0;
        } else {
            m_doilo[d] = 0;
            m_doihi[d] = a_out_rad;
        }
    }

    Box operator() (const Box& a_bx) const noexcept {
        IntVect lo = amrex::coarsen(a_bx.smallEnd(), m_crse_ratio);
        IntVect hi = amrex::coarsen(a_bx.bigEnd(), m_crse_ratio);
        const int d = m_face.coordDir();
        if (m_face.isLow()) {
            hi[d] = lo[d];
        } else {
            lo[d] = hi[d];
        }
        lo += m_loshft;
        hi += m_hishft;
        return Box(lo,hi,m_typ);
    }

    Box coarsen (Box const& a_box) const noexcept { return amrex::coarsen(a_box,m_crse_ratio); }

    IntVect doiLo () const noexcept { return m_doilo; }
    IntVect doiHi () const noexcept { return m_doihi; }

    IndexType index_type () const noexcept { return m_typ; }
    IntVect coarsen_ratio () const noexcept { return m_crse_ratio; }

    friend bool operator== (BATbndryReg const& a, BATbndryReg const& b) noexcept {
        return a.m_face == b.m_face and a.m_typ == b.m_typ and a.m_crse_ratio == b.m_crse_ratio
            and a.m_loshft == b.m_loshft and a.m_hishft == b.m_hishft
            and a.m_doilo == b.m_doilo and a.m_doihi == b.m_doihi;
    }

    Orientation m_face;
    IndexType m_typ;
    IntVect m_crse_ratio;
    IntVect m_loshft;
    IntVect m_hishft;
    IntVect m_doilo;
    IntVect m_doihi;
};

struct BATransformer
{
    enum struct BATType { null, indexType, coarsenRatio, indexType_coarsenRatio, bndryReg};

    union BATOp {
        BATOp () noexcept
            : m_null() {}
        BATOp (IndexType t) noexcept
            : m_indexType(t) {}
        BATOp (IntVect const& r) noexcept
            : m_coarsenRatio(r) {}
        BATOp (IndexType t, IntVect const& r) noexcept
            : m_indexType_coarsenRatio(t,r) {}
        BATOp (Orientation f, IndexType t, int in_rad, int out_rad, int extent_rad) noexcept
            : m_bndryReg(f,t,in_rad,out_rad,extent_rad) {}
        BATnull  m_null;
        BATindexType m_indexType;
        BATcoarsenRatio m_coarsenRatio;
        BATindexType_coarsenRatio m_indexType_coarsenRatio;
        BATbndryReg m_bndryReg;
    };

    BATransformer ()
        : m_bat_type(BATType::null), m_op() {}

    BATransformer (IndexType t)
        : m_bat_type(t.cellCentered() ? BATType::null : BATType::indexType),
          m_op      (t.cellCentered() ? BATOp()       : BATOp(t)) {}

    BATransformer (Orientation f, IndexType t, int in_rad, int out_rad, int extent_rad)
        : m_bat_type(BATType::bndryReg),
          m_op(f,t,in_rad,out_rad,extent_rad) {}

    Box operator() (Box const& ab) const noexcept {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null(ab);
        case     BATType::indexType:
            return m_op.m_indexType(ab);
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio(ab);
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio(ab);
        default:
            return m_op.m_bndryReg(ab);
        }
    }

    Box coarsen (Box const& a_box) const noexcept  {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null.coarsen(a_box);
        case     BATType::indexType:
            return m_op.m_indexType.coarsen(a_box);
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio.coarsen(a_box);
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio.coarsen(a_box);
        default:
            return m_op.m_bndryReg.coarsen(a_box);
        }
    }

    IntVect doiLo () const noexcept  {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null.doiLo();
        case     BATType::indexType:
            return m_op.m_indexType.doiLo();
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio.doiLo();
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio.doiLo();
        default:
            return m_op.m_bndryReg.doiLo();
        }
    }

    IntVect doiHi () const noexcept  {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null.doiHi();
        case     BATType::indexType:
            return m_op.m_indexType.doiHi();
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio.doiHi();
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio.doiHi();
        default:
            return m_op.m_bndryReg.doiHi();
        }
    }

    IndexType index_type () const noexcept  {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null.index_type();
        case     BATType::indexType:
            return m_op.m_indexType.index_type();
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio.index_type();
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio.index_type();
        default:
            return m_op.m_bndryReg.index_type();
        }
    }

    IntVect coarsen_ratio () const noexcept  {
        switch (m_bat_type)
        {
        case     BATType::null:
            return m_op.m_null.coarsen_ratio();
        case     BATType::indexType:
            return m_op.m_indexType.coarsen_ratio();
        case     BATType::coarsenRatio:
            return m_op.m_coarsenRatio.coarsen_ratio();
        case     BATType::indexType_coarsenRatio:
            return m_op.m_indexType_coarsenRatio.coarsen_ratio();
        default:
            return m_op.m_bndryReg.coarsen_ratio();
        }
    }

    bool is_null () const noexcept {
        return m_bat_type == BATType::null;
    }

    bool is_simple () const noexcept {
        return m_bat_type != BATType::bndryReg;
    }

    void set_coarsen_ratio (IntVect const& a_ratio) noexcept {
        switch (m_bat_type)
        {
        case BATType::null:
        {
            if (a_ratio == IntVect::TheUnitVector()) {
                return;
            } else {
                m_bat_type = BATType::coarsenRatio;
                m_op.m_coarsenRatio.m_crse_ratio = a_ratio;
                return;
            }
        }
        case BATType::indexType:
        {
            if (a_ratio == IntVect::TheUnitVector()) {
                return;
            } else {
                m_bat_type = BATType::indexType_coarsenRatio;
                auto t = m_op.m_indexType.m_typ;
                m_op.m_indexType_coarsenRatio.m_typ = t;
                m_op.m_indexType_coarsenRatio.m_crse_ratio = a_ratio;
                return;
            }
        }
        case BATType::coarsenRatio:
        {
            if (a_ratio == IntVect::TheUnitVector()) {
                m_bat_type = BATType::null;
                return;
            } else {
                m_op.m_coarsenRatio.m_crse_ratio = a_ratio;
                return;
            }
        }
        case BATType::indexType_coarsenRatio:
        {
            if (a_ratio == IntVect::TheUnitVector()) {
                m_bat_type = BATType::indexType;
                auto t = m_op.m_indexType_coarsenRatio.m_typ;
                m_op.m_indexType.m_typ = t;
                return;
            } else {
                m_op.m_indexType_coarsenRatio.m_crse_ratio = a_ratio;
                return;
            }
        }
        default:
        {
            m_op.m_bndryReg.m_crse_ratio = a_ratio;
            return;
        }
        }
    }

    void set_index_type (IndexType typ) noexcept {
        switch (m_bat_type)
        {
        case BATType::null:
        {
            if (typ.cellCentered()) {
                return;
            } else {
                m_bat_type = BATType::indexType;
                m_op.m_indexType.m_typ = typ;
                return;
            }
        }
        case BATType::indexType:
        {
            if (typ.cellCentered()) {
                m_bat_type = BATType::null;
                return;
            } else {
                m_op.m_indexType.m_typ = typ;
                return;
            }
        }
        case BATType::coarsenRatio:
        {
            if (typ.cellCentered()) {
                return;
            } else {
                m_bat_type = BATType::indexType_coarsenRatio;
                auto r = m_op.m_coarsenRatio.m_crse_ratio;
                m_op.m_indexType_coarsenRatio.m_typ = typ;
                m_op.m_indexType_coarsenRatio.m_crse_ratio = r;
                return;
            }
        }
        case BATType::indexType_coarsenRatio:
        {
            if (typ.cellCentered()) {
                m_bat_type = BATType::coarsenRatio;
                auto r = m_op.m_indexType_coarsenRatio.m_crse_ratio;
                m_op.m_coarsenRatio.m_crse_ratio = r;
                return;
            } else {
                m_op.m_indexType_coarsenRatio.m_typ = typ;
                return;
            }
        }
        default:
        {
            m_op.m_bndryReg.m_typ = typ;
            return;
        }
        }
    }

    friend bool operator== (BATransformer const& a, BATransformer const& b) noexcept {
        switch (a.m_bat_type)
        {
        case BATType::bndryReg:
            switch (b.m_bat_type)
            {
            case BATType::bndryReg:
                return a.m_op.m_bndryReg == b.m_op.m_bndryReg;
            default:
                return false;
            }
        default:
            switch (b.m_bat_type)
            {
            case BATType::bndryReg:
                return false;
            default:
                return a.index_type() == b.index_type()
                    and a.coarsen_ratio() == b.coarsen_ratio();
            }
        }
    }

    BATType m_bat_type;
    BATOp m_op;
};

// for backward compatibility
using BndryBATransformer = BATransformer;

/**
* \brief A collection of Boxes stored in an Array.  It is a
* reference-counted concrete class, not a polymorphic one; i.e. you
* cannot use any of the List member functions with a BoxList.
*/

class MFIter;
class AmrMesh;
class FabArrayBase;

class BoxArray
{
public:
    //! Construct an empty BoxArray.
    BoxArray ();

    //! The copy constructor.
    BoxArray (const BoxArray& rhs);

    //! The move constructor.
    BoxArray (BoxArray&& rhs) noexcept = default;

    //! The move assignment operator.
    BoxArray& operator= (BoxArray const& rhs) noexcept = default;
    BoxArray& operator= (BoxArray&& rhs) noexcept = default;

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

    BoxArray (BoxList&& bl, IntVect const& max_grid_size);
    
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
    void resize (Long len);

    //! Return the number of boxes in the BoxArray.
    Long size () const noexcept { return m_ref->m_abox.size(); }

    //! Return the number of boxes that can be held in the current allocated storage
    Long capacity () const noexcept { return m_ref->m_abox.capacity(); }

    //! Return whether the BoxArray is empty
    bool empty () const noexcept { return m_ref->m_abox.empty(); }

    //! Returns the total number of cells contained in all boxes in the BoxArray.
    Long numPts() const noexcept;

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
        return m_bat(m_ref->m_abox[index]);
    }

    //! Return element index of this BoxArray.
    Box operator[] (const MFIter& mfi) const noexcept;

    //! Return element index of this BoxArray.
    Box get (int index) const noexcept { return operator[](index); }

    //! Return cell-centered box at element index of this BoxArray.
    Box getCellCenteredBox (int index) const noexcept {
        return m_bat.coarsen(m_ref->m_abox[index]);
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
    Box minimalBox (Long& npts_avg_box) const;

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
    IndexType ixType () const noexcept { return m_bat.index_type(); }

    //! Return crse ratio of this BoxArray
    IntVect crseRatio () const noexcept { return m_bat.coarsen_ratio(); }

    static void Initialize ();
    static void Finalize ();
    static bool initialized;

    //! Make ourselves unique.
    void uniqify ();

    friend class AmrMesh;
    friend class FabArrayBase;

private:
    //!  Update BoxArray index type according the box type, and then convert boxes to cell-centered.
    void type_update ();

    BoxList const& simplified_list () const; // For regular AMR grids only
    BoxArray simplified () const;

    BARef::HashType& getHashMap () const;

    IntVect getDoiLo () const noexcept;
    IntVect getDoiHi () const noexcept;

    BATransformer m_bat;
    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<BARef> m_ref;
    mutable std::shared_ptr<BoxList> m_simplified_list;
};

//! Write a BoxArray to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const BoxArray& ba);

std::ostream& operator<< (std::ostream& os, const BoxArray::RefID& id);

}

#endif /*BL_BOXARRAY_H*/
