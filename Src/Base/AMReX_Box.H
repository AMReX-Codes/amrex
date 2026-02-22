
#ifndef AMREX_BOX_H_
#define AMREX_BOX_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Orientation.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Math.H>

#include <iosfwd>

/**
 * \file AMReX_Box.H
 * \brief Integer-lattice boxes and helpers for defining index-space regions.
 */

namespace amrex
{
template<int dim>
class BoxND;
//! \brief Box is an alias for \ref amrex::BoxND instantiated with <tt>AMREX_SPACEDIM</tt>.
using Box = BoxND<AMREX_SPACEDIM>;
class BoxCommHelper;
template<int dim>
class BoxIteratorND;

/**
* \class amrex::BoxND
* \ingroup amrex_index_space
* \brief A Rectangular Domain on an Integer Lattice
*
* A BoxND is an abstraction for defining discrete regions of
* dim indexing space. Boxes have an IndexType, which defines
* IndexType::CELL or IndexType::NODE based points for each direction
* and a low and high INTVECT which defines the lower and upper corners
* of the BoxND. Boxes can exist in positive and negative indexing space.
*
* Box is a dimension dependent alias to BoxND<AMREX_SPACEDIM>, so AMREX_SPACEDIM must be
* defined as either 1, 2, or 3 when compiling.
*/
template<int dim>
class BoxND
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<Box>::type();
    friend class BoxCommHelper;

public:
    /**
     * \brief Default constructor creates an invalid box (use ok() to test).
     *
     * DO NOT CHANGE THIS BEHAVIOR!
     */
    AMREX_GPU_HOST_DEVICE
    constexpr BoxND () noexcept
        : smallend(1),
          bigend(0)
        {}

    //! Construct a cell-centered box spanning \p small through \p big.
    AMREX_GPU_HOST_DEVICE
    constexpr BoxND (const IntVectND<dim>& small, const IntVectND<dim>& big) noexcept
        : smallend(small),
          bigend(big)
        {}

    /**
     * \brief Construct a box with lower corner \p small and lengths \p vec_len.
     * \param small Lower corner.
     * \param vec_len Array of lengths in each direction.
     */
    AMREX_GPU_HOST_DEVICE
    BoxND (const IntVectND<dim>& small, const int* vec_len) noexcept
        : smallend(small),
          bigend(small + IntVectND<dim>(vec_len) - 1)
        {}

    /**
     * \brief Construct a box with explicit nodal/cell flags per direction.
     *
     * \param small Lower corner.
     * \param big Upper corner.
     * \param typ IntVect of 0 (cell) / 1 (node) selectors; bounds must match.
     */
    AMREX_GPU_HOST_DEVICE
    BoxND (const IntVectND<dim>& small, const IntVectND<dim>& big, const IntVectND<dim>& typ) noexcept
        : smallend(small),
          bigend(big),
          btype(typ)
        {
            BL_ASSERT(typ.allGE(0) && typ.allLE(1));
        }

    //! Construct a box with the given index type descriptor.
    AMREX_GPU_HOST_DEVICE
    BoxND (const IntVectND<dim>& small, const IntVectND<dim>& big, IndexTypeND<dim> t) noexcept
        : smallend(small),
          bigend(big),
          btype(t)
        {}

    //! Construct a box from an Array4 view (bounds are inferred from begin/end).
    template <typename T, int Tdim=dim, std::enable_if_t<( 1<=Tdim && Tdim<=3 ), int> = 0>
    AMREX_GPU_HOST_DEVICE
    explicit BoxND (Array4<T> const& a) noexcept
        : smallend(IntVectShrink<dim>(a.begin)),
          bigend(IntVectShrink<dim>(a.end) - 1)
        {}

    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    //! Return the inclusive lower bound of the box.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const IntVectND<dim>& smallEnd () const& noexcept { return smallend; }

    /// \cond DOXYGEN_IGNORE
    [[nodiscard]] const IntVectND<dim>& smallEnd () && = delete;
    /// \endcond

    //! Return the inclusive lower bound in the given direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int smallEnd (int dir) const& noexcept { return smallend[dir]; }

    //! Return the inclusive upper bound of the box.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const IntVectND<dim>& bigEnd () const& noexcept { return bigend; }

    /// \cond DOXYGEN_IGNORE
    [[nodiscard]] const IntVectND<dim>& bigEnd () && = delete;
    /// \endcond

    //! Return the inclusive upper bound in the given direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int bigEnd (int dir) const noexcept { return bigend[dir]; }

    //! Return the indexing type.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IndexTypeND<dim> ixType () const noexcept { return btype; }

    //! Return the indexing type.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVectND<dim> type () const noexcept { return btype.ixType(); }

    //! Return the indexing type in the specified direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IndexType::CellIndex type (int dir) const noexcept { return btype.ixType(dir); }

    //! Return the length of the BoxND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVectND<dim> size () const noexcept
    {
        return bigend - smallend + 1;
    }

    //! Return the length of the BoxND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVectND<dim> length () const noexcept
    {
        return bigend - smallend + 1;
    }

    //! Return the length of the BoxND in given direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int length (int dir) const noexcept { return bigend[dir] - smallend[dir] + 1; }

    //! Return lengths padded to 3D (unused components set to 1).
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> length3d () const noexcept {
        Dim3 len3d = length().dim3(1);
        return {{len3d.x, len3d.y, len3d.z}};
    }

    //! Return low-end coordinates padded to 3D (missing dims default to 0).
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> loVect3d () const noexcept {
        Dim3 lo3d = smallend.dim3(0);
        return {{lo3d.x, lo3d.y, lo3d.z}};
    }

    //! Return high-end coordinates padded to 3D (missing dims default to 0).
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> hiVect3d () const noexcept {
        Dim3 hi3d = bigend.dim3(0);
        return {{hi3d.x, hi3d.y, hi3d.z}};
    }

    //! Return a pointer to the low-end coordinates (useful for Fortran calls).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const int* loVect () const& noexcept { return smallend.getVect(); }
    AMREX_GPU_HOST_DEVICE
    const int* loVect () && = delete;
    //! Return a pointer to the high-end coordinates (useful for Fortran calls).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const int* hiVect () const& noexcept { return bigend.getVect(); }
    AMREX_GPU_HOST_DEVICE
    const int* hiVect () && = delete;

    //! Return the coordinate normal to the given face.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int operator[] (Orientation face) const noexcept {
        const int dir = face.coordDir();
        return face.isLow() ? smallend[dir] : bigend[dir];
    }

    //! Checks if it is an empty BoxND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool isEmpty () const noexcept { return !ok(); }

    //! Checks if it is a proper BoxND (including a valid type).
    //! Return true if high bounds are ≥ low bounds and the index type is valid.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool ok () const noexcept { return bigend.allGE(smallend) && btype.ok(); }

    //! Return true if argument is contained within BoxND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const IntVectND<dim>& p) const noexcept {
        return p.allGE(smallend) && p.allLE(bigend);
    }

    //! Return true if argument is contained within BoxND.
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const Dim3& p) const noexcept {
        IntVectND<dim> piv{p};
        return contains(piv);
    }

    //! Return true if argument is contained within BoxND.
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (int i, int j, int k) const noexcept {
        Dim3 p3d{i, j, k};
        return contains(p3d);
    }

    /** \brief Return true if argument is contained within BoxND.
    * It is an error if the Boxes have different types.
    */
    /// \return True when \p b lies within this box.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const BoxND& b) const noexcept
    {
        BL_ASSERT(sameType(b));
        return b.smallend.allGE(smallend) && b.bigend.allLE(bigend);
    }

    //! Return true if argument is strictly contained within BoxND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool strictly_contains (const IntVectND<dim>& p) const noexcept {
        return p.allGT(smallend) && p.allLT(bigend);
    }

    /**
    * \brief Return true if argument is strictly contained within BoxND.
    * It is an error if the Boxes have different types.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool strictly_contains (const BoxND& b) const noexcept
    {
        BL_ASSERT(sameType(b));
        return b.smallend.allGT(smallend) && b.bigend.allLT(bigend);
    }

    //! Return true if argument is strictly contained within BoxND.
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool strictly_contains (const Dim3& p) const noexcept {
        IntVectND<dim> piv{p};
        return strictly_contains(piv);
    }

    //! Return true if argument is strictly contained within BoxND.
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool strictly_contains (int i, int j, int k) const noexcept {
        Dim3 p3d{i, j, k};
        return strictly_contains(p3d);
    }

    /**
    * \brief Return true if Boxes have non-null intersections.
    * It is an error if the Boxes have different types.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool intersects (const BoxND& b) const noexcept { BoxND isect(*this); isect &= b; return isect.ok(); }

    /**
    * \brief Return true is Boxes same size, ie translates of each other,.
    * It is an error if they have different types.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool sameSize (const BoxND& b) const noexcept {
        BL_ASSERT(sameType(b));
        return length() == b.length();
    }

    //! Return true if Boxes have same type.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool sameType (const BoxND &b) const noexcept { return btype == b.btype; }

    //! Return true if Boxes are identical (including type).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator== (const BoxND& b) const noexcept { return smallend == b.smallend && bigend == b.bigend && b.btype == btype; }

    //! Return true if Boxes differ (including type).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator!= (const BoxND& b) const noexcept { return !operator==(b); }

    //! Lexicographic ordering (type, then lower, then upper corners).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator< (const BoxND& rhs) const noexcept
    {
        return btype < rhs.btype ||
             ((btype == rhs.btype) &&
                 (  (smallend < rhs.smallend) ||
                   ((smallend == rhs.smallend) && (bigend < rhs.bigend)) ));
    }
    //! Return true if this box is <= rhs in the lexicographic ordering.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator <= (const BoxND& rhs) const noexcept {
        return !(rhs < *this);
    }
    //! Return true if this box is > rhs in the lexicographic ordering.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator> (const BoxND& rhs) const noexcept {
        return rhs < *this;
    }
    //! Return true if this box is >= rhs in the lexicographic ordering.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool operator>= (const BoxND& rhs) const noexcept {
        return !(*this < rhs);
    }

    //! Return true if BoxND is cell-centered in all indexing directions.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool cellCentered () const noexcept { return !btype.any(); }

    /// \cond DOXYGEN_IGNORE
    //! Assert that there are no int/Long overflows when calling length or numPts.
    void checkOverflow () const noexcept {
        if (ok()) {
            for (int i = 0; i < dim; ++i) {
                auto lo = static_cast<Long>(smallend[i]);
                auto hi = static_cast<Long>(bigend[i]);
                Long len = hi - lo + 1;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(len>=0 && len<std::numeric_limits<int>::max(),
                    "Overflow when computing length of box");
            }
            auto num_pts = static_cast<Long>(length(0));
            for (int i = 1; i < dim; ++i) {
                auto len = static_cast<Long>(length(i));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(num_pts == 0 || len == 0 ||
                    num_pts <= std::numeric_limits<Long>::max() / len,
                    "Overflow when computing numPts of box");
                num_pts *= len;
            }
        }
    }
    /// \endcond

    /**
    * \brief Return the number of points contained in the BoxND.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Long numPts () const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
        AMREX_IF_ON_HOST((checkOverflow();))
#endif
        if (ok()) {
            auto num_pts = static_cast<Long>(length(0));
            for (int i = 1; i < dim; ++i) {
                num_pts *= static_cast<Long>(length(i));
            }
            return num_pts;
        } else {
            return Long(0);
        }
    }

    /**
    * \brief Return the number of points contained in the BoxND.
    * This is intended for use only in diagnostic messages.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    double d_numPts () const noexcept {
        if (ok()) {
            auto num_pts = static_cast<double>(length(0));
            for (int i = 1; i < dim; ++i) {
                num_pts *= static_cast<double>(length(i));
            }
            return num_pts;
        } else {
            return 0.0;
        }
    }

    /**
    * \brief Return the volume, in indexing space, of region enclosed by
    * this BoxND.   This is identical to numPts() for CELL centered
    * BoxND; otherwise, numPts() > volume().
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Long volume () const noexcept {
        if (ok()) {
            auto num_pts = static_cast<Long>(length(0)-btype[0]);
            for (int i = 1; i < dim; ++i) {
                num_pts *= static_cast<Long>(length(i)-btype[i]);
            }
            return num_pts;
        } else {
            return Long(0);
        }
    }

    /**
    * \brief Return length of longest side.  dir is modified to give
    * direction with longest side: 0...dim-1. Ignores type.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int longside (int& dir) const noexcept {
        int maxlen = length(0);
        dir = 0;
        for (int i = 1; i < dim; i++)
        {
            if (length(i) > maxlen)
            {
                maxlen = length(i);
                dir = i;
            }
        }
        return maxlen;
    }

    //! Return length of longest side.  Ignores type.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int longside () const noexcept {
        int ignore = 0;
        return longside(ignore);
    }

    /**
    * \brief Return length of shortest side.  dir is modified to give
    *  direction with shortest side: 0...dim-1.  Ignores type.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int shortside (int& dir) const noexcept {
        int minlen = length(0);
        dir = 0;
        for (int i = 1; i < dim; i++)
        {
            if (length(i) < minlen)
            {
                minlen = length(i);
                dir = i;
            }
        }
        return minlen;
    }

    //! Return length of shortest side.  Ignores type.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int shortside () const noexcept {
        int ignore = 0;
        return shortside(ignore);
    }

    /**
    * \brief Return offset of point from smallend; i.e.
    * index(smallend) -> 0, bigend would return numPts()-1.
    * Is used in accessing FArrayBox.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Long index (const IntVectND<dim>& v) const noexcept;

    //! Given the offset, compute IntVectND<dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVectND<dim> atOffset (Long offset) const noexcept;

    //! Given the offset, return its corresponding index vector in
    //! GpuArray<int,3>. This function is available for 1D, 2D and 3D Boxes
    //! only,
    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> atOffset3d (Long offset) const noexcept;

    //! Redefine the small end of the BoxND.
    AMREX_GPU_HOST_DEVICE
    BoxND& setSmall (const IntVectND<dim>& sm) noexcept { smallend = sm; return *this; }

    //! Redefine the small end of the BoxND.
    AMREX_GPU_HOST_DEVICE
    BoxND& setSmall (int dir, int sm_index) noexcept { smallend.setVal(dir,sm_index); return *this; }

    //! Redefine the big end of the BoxND.
    AMREX_GPU_HOST_DEVICE
    BoxND& setBig (const IntVectND<dim>& bg) noexcept { bigend = bg; return *this; }

    //! Redefine the big end of the BoxND.
    AMREX_GPU_HOST_DEVICE
    BoxND& setBig (int dir, int bg_index) noexcept { bigend.setVal(dir,bg_index); return *this; }

    /**
    * \brief Set the entire range in a given direction, starting at
    * sm_index with length n_cells.  NOTE: This will yield an
    * illegal BoxND if n_cells <= 0.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& setRange (int dir,
                   int sm_index,
                   int n_cells = 1) noexcept;

    //! Set indexing type
    AMREX_GPU_HOST_DEVICE
    BoxND& setType (const IndexTypeND<dim>& t) noexcept { btype = t; return *this; }

    //! Shift this BoxND nzones indexing positions in coordinate direction dir.
    AMREX_GPU_HOST_DEVICE
    BoxND& shift (int dir, int nzones) noexcept { smallend.shift(dir,nzones); bigend.shift(dir,nzones); return *this; }

    //! Equivalent to b.shift(0,iv[0]).shift(1,iv[1]) ....
    AMREX_GPU_HOST_DEVICE
    BoxND& shift (const IntVectND<dim>& iv) noexcept { smallend.shift(iv); bigend.shift(iv); return *this; }

    /**
    * \brief This member shifts the BoxND by "half" indices, thereby
    * converting the BoxND from type CELL to NODE and visa-versa.
    * b.shiftHalf(0,1)  shifts b to the right by 1/2 cells.
    * b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
    * NOTE: If num_halfs is EVEN the shift is num_halfs/2 full
    * zones and hence will not change the type.
    * This is: b.shifthalf(4) == b.shift(2).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& shiftHalf (int dir, int num_halfs) noexcept;

    //! Equivalent to b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...
    AMREX_GPU_HOST_DEVICE
    BoxND& shiftHalf (const IntVectND<dim>& iv) noexcept;

    /**
    * \brief Convert the BoxND from the current type into the
    * argument type.  This may change the BoxND coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& convert (IndexTypeND<dim> typ) noexcept;

    /**
    * \brief Convert the BoxND from the current type into the
    * argument type.  This may change the BoxND coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& convert (const IntVectND<dim>& typ) noexcept;

    //! Convert to NODE type in all directions.
    AMREX_GPU_HOST_DEVICE
    BoxND& surroundingNodes () noexcept;

    //! Convert to NODE type in given direction.
    AMREX_GPU_HOST_DEVICE
    BoxND& surroundingNodes (int dir) noexcept;

    //! Convert to NODE type in given direction.
    AMREX_GPU_HOST_DEVICE
    BoxND& surroundingNodes (Direction d) noexcept { return surroundingNodes(static_cast<int>(d)); }

    //! Convert to CELL type in all directions.
    AMREX_GPU_HOST_DEVICE
    BoxND& enclosedCells () noexcept;

    //! Convert to CELL type in given direction.
    AMREX_GPU_HOST_DEVICE
    BoxND& enclosedCells (int dir) noexcept;

    //! Convert to CELL type in given Direction `d`.
    AMREX_GPU_HOST_DEVICE
    BoxND& enclosedCells (Direction d) noexcept { return enclosedCells(static_cast<int>(d)); }

    /**
    * \brief Return BoxND that is intersection of this BoxND
    * and argument.  The Boxes MUST be of same type.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND operator& (const BoxND& rhs) const noexcept { BoxND lhs(*this); lhs &= rhs; return lhs; }

    //! Intersect this BoxND with its argument. The Boxes MUST be of the same type.
    AMREX_GPU_HOST_DEVICE
    BoxND& operator&= (const BoxND& rhs) noexcept
    {
        BL_ASSERT(sameType(rhs));
        smallend.max(rhs.smallend);
        bigend.min(rhs.bigend);
        return *this;
    }

    /**
    * \brief Modify BoxND to that of the minimum BoxND containing both
    * the original BoxND and the argument.
    * Both Boxes must have identical type.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& minBox (const BoxND& b) noexcept {
        // BoxArray may call this with not ok boxes.  BL_ASSERT(b.ok() && ok());
        BL_ASSERT(sameType(b));
        smallend.min(b.smallend);
        bigend.max(b.bigend);
        return *this;
    }

    //! Shift BoxND (relative) by given IntVectND<dim>.
    AMREX_GPU_HOST_DEVICE
    BoxND& operator+= (const IntVectND<dim>& v) noexcept { smallend += v; bigend   += v; return *this; }

    //! Shift BoxND (relative) by given IntVectND<dim>.
    AMREX_GPU_HOST_DEVICE
    BoxND operator+ (const IntVectND<dim>& v) const noexcept { BoxND r(*this); r += v; return r; }

    //! Shift BoxND (relative) by given IntVectND<dim>.
    AMREX_GPU_HOST_DEVICE
    BoxND& operator-= (const IntVectND<dim>& v) noexcept { smallend -= v; bigend -= v; return *this; }

    //! Shift BoxND (relative) by given IntVectND<dim>.
    AMREX_GPU_HOST_DEVICE
    BoxND operator- (const IntVectND<dim>& v) const noexcept { BoxND r(*this); r -= v; return r; }

    /**
    * \brief Chop the BoxND at \p chop_pnt in the \p dir direction
    * return one BoxND and modify the object BoxND.
    * The union of the two is the original BoxND.
    * The modified BoxND is the low end, the returned BoxND
    * is the high end.  If type(dir) = CELL, the Boxes are disjoint
    * with the chop_pnt included in the high end (new BoxND).
    * It is an ERROR if chop_pnt is the low end of the orig BoxND.
    * If type(dir) = NODE, the chop_pnt is included in both Boxes
    * but is the only point in common.  It is also an error if the
    * chop_pnt is an end node of the BoxND.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND chop (int dir, int chop_pnt) noexcept;

    //! Grow in all directions by \p i cells (negative shrinks).
    AMREX_GPU_HOST_DEVICE
    BoxND& grow (int i) noexcept { smallend.diagShift(-i); bigend.diagShift(i); return *this; }

    //! Grow in each direction by the supplied vector (negative shrinks).
    AMREX_GPU_HOST_DEVICE
    BoxND& grow (const IntVectND<dim>& v) noexcept { smallend -= v; bigend += v; return *this;}

    /**
    * \brief Grow on both ends in direction \p idir by \p n_cell cells (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& grow (int idir, int n_cell) noexcept { smallend.shift(idir, -n_cell); bigend.shift(idir, n_cell); return *this; }

    /**
    * \brief Grow on both ends in direction d by \p n_cell cells (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& grow (Direction d, int n_cell) noexcept { return grow(static_cast<int>(d), n_cell); }

    /**
    * \brief Grow on the low end by \p n_cell cells in direction \p idir (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& growLo (int idir, int n_cell = 1) noexcept { smallend.shift(idir, -n_cell); return *this; }

    /**
    * \brief Grow on the low end by \p n_cell cells in direction \p d (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& growLo (Direction d, int n_cell = 1) noexcept { return growLo(static_cast<int>(d), n_cell); }

    /**
    * \brief Grow on the high end by \p n_cell cells in direction \p idir (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& growHi (int idir, int n_cell = 1) noexcept { bigend.shift(idir,n_cell); return *this; }

    /**
    * \brief Grow on the high end by \p n_cell cells in direction \p d (negative shrinks).
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& growHi (Direction d, int n_cell = 1) noexcept { return growHi(static_cast<int>(d), n_cell); }

    //! Grow in the direction of the given face.
    AMREX_GPU_HOST_DEVICE
    BoxND& grow (Orientation face, int n_cell = 1) noexcept {
        int idir = face.coordDir();
        if (face.isLow()) {
            smallend.shift(idir, -n_cell);
        } else {
            bigend.shift(idir,n_cell);
        }
        return *this;
    }

    /**
    * \brief Refine BoxND by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& refine (int ref_ratio) noexcept {
        return this->refine(IntVectND<dim>(ref_ratio));
    }

    /*
    * \brief Refine BoxND by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& refine (const IntVectND<dim>& ref_ratio) noexcept;

    /**
    * \brief Coarsen BoxND by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened BoxND must contain
    * the original BoxND.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& coarsen (int ref_ratio) noexcept {
        return this->coarsen(IntVectND<dim>(ref_ratio));
    }

    /**
    * \brief Coarsen BoxND by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened BoxND must contain
    * the original BoxND.
    */
    AMREX_GPU_HOST_DEVICE
    BoxND& coarsen (const IntVectND<dim>& ref_ratio) noexcept;

    /**
    * \brief Step through the rectangle.  It is a runtime error to give
    * a point not inside rectangle.  Iteration may not be efficient.
    */
    AMREX_GPU_HOST_DEVICE
    void next (IntVectND<dim> &) const noexcept;

    /**
    * \brief This static member function returns a constant reference to
    * an object of type BoxND representing the unit BoxND in
    * dim-dimensional space.
    */
    AMREX_GPU_HOST_DEVICE
    static BoxND TheUnitBox () noexcept {
        return BoxND(IntVectND<dim>::TheZeroVector(),IntVectND<dim>::TheZeroVector());
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool isSquare() const noexcept;

    /**
     * \brief Return whether this Box is coarsenable.
     *
     * A coasrenable Box means `box.coarsen(ratio).refine(ratio) == box`. A
     * Box is considered uncoarsenable if the coarsened Box would have a
     * length smaller than the given minimum width.
     *
     * \param refrat coarsening ratio
     * \param min_width minimum width of coarsened Box.
     * \return whether this Box is coarsenable.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool coarsenable (const IntVectND<dim>& refrat, const IntVectND<dim>& min_width) const noexcept
    {
        if (!size().allGE(refrat*min_width)) {
            return false;
        } else {
            BoxND testBox = *this;
            testBox.coarsen(refrat);
            testBox.refine (refrat);
            return (*this == testBox);
        }
    }

    /**
     * \brief Return whether this Box is coarsenable.
     *
     * A coasrenable Box means `box.coarsen(ratio).refine(ratio) == box`. A
     * Box is considered uncoarsenable if the coarsened Box would have a
     * length smaller than the given minimum width.
     *
     * \param refrat coarsening ratio
     * \param min_width minimum width of coarsened Box.
     * \return whether this Box is coarsenable.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool coarsenable (int refrat, int min_width=1) const noexcept {
        return coarsenable(IntVectND<dim>(refrat), IntVectND<dim>(min_width));
    }

    /**
     * \brief Return whether this Box is coarsenable.
     *
     * A coasrenable Box means `box.coarsen(ratio).refine(ratio) == box`. A
     * Box is considered uncoarsenable if the coarsened Box would have a
     * length smaller than the given minimum width.
     *
     * \param refrat coarsening ratio
     * \param min_width minimum width of coarsened Box.
     * \return whether this Box is coarsenable.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool coarsenable (const IntVectND<dim>& refrat, int min_width=1) const noexcept
    {
        return coarsenable(refrat, IntVectND<dim>(min_width));
    }

    //! Ensure each dimension has at least one cell or node by expanding empty directions.
    AMREX_GPU_HOST_DEVICE
    void normalize () noexcept
    {
        for (int idim=0; idim < dim; ++idim) {
            if (this->length(idim) == 0) {
                this->growHi(idim,1);
            }
        }
    }

    //! Collapse the box to a single slab at coordinate \p slab_index along \p direction.
    AMREX_GPU_HOST_DEVICE
    BoxND& makeSlab (int direction, int slab_index) noexcept
    {
        smallend[direction] = slab_index;
        bigend[direction] = slab_index;
        return *this;
    }

    /**
    * \brief Return a BoxIteratorND that can be used to loop over the IntVects
    * contained by the BoxND.
    * \code{.cpp}
    *   Box b;
    *   ...
    *   for (IntVect iv : b.iterator()) {
    *       // do operations involving iv
    *   }
    * \endcode
    */
    [[nodiscard]] BoxIteratorND<dim> iterator () const noexcept {
        return BoxIteratorND<dim>{*this};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr std::size_t ndims () noexcept {
        return static_cast<std::size_t>(dim);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int indims () noexcept {
        return dim;
    }

    /**
    * \brief Return a new BoxND of dimension new_dim and
    * assigns the first new_dim dimension of this BoxND to it.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    BoxND<new_dim> shrink () const noexcept {
        static_assert(new_dim <= dim);
        auto lo = smallend.template shrink<new_dim>();
        auto hi = bigend.template shrink<new_dim>();
        auto typ = btype.template shrink<new_dim>();
        return BoxND<new_dim>(lo, hi, typ);
    }

    /**
    * \brief Return a new BoxND of size new_dim and
    * assigns all values of this BoxND to it and
    * (small=0, big=0, typ=CELL) to the remaining elements.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    BoxND<new_dim> expand () const noexcept {
        static_assert(new_dim >= dim);
        auto lo = smallend.template expand<new_dim>(0);
        auto hi = bigend.template expand<new_dim>(0);
        auto typ = btype.template expand<new_dim>(IndexType::CellIndex::CELL);
        return BoxND<new_dim>(lo, hi, typ);
    }

    /**
    * \brief Return a new BoxND of size new_dim
    * by either shrinking or expanding this BoxND
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    BoxND<new_dim> resize () const noexcept {
        if constexpr (new_dim > dim) {
            return expand<new_dim>();
        } else {
            return shrink<new_dim>();
        }
    }

private:
    IntVectND<dim>         smallend;
    IntVectND<dim>         bigend;
    IndexTypeND<dim>       btype;
};

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::refine (const IntVectND<dim>& ref_ratio) noexcept
{
    if (ref_ratio != 1) {
        IntVectND<dim> shft(1);
        shft -= btype.ixType();
        smallend *= ref_ratio;
        bigend += shft;
        bigend *= ref_ratio;
        bigend -= shft;
    }
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::coarsen (const IntVectND<dim>& ref_ratio) noexcept
{
    if (ref_ratio != 1)
    {
        smallend.coarsen(ref_ratio);

        if (btype.any())
        {
            IntVectND<dim> off(0);
            for (int dir = 0; dir < dim; dir++)
            {
                if (btype[dir]) {
                    if (bigend[dir]%ref_ratio[dir]) {
                        off.setVal(dir,1);
                    }
                }
            }
            bigend.coarsen(ref_ratio);
            bigend += off;
        }
        else
        {
            bigend.coarsen(ref_ratio);
        }
    }

    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::convert (const IntVectND<dim>& typ) noexcept
{
    BL_ASSERT(typ.allGE(0) && typ.allLE(1));
    IntVectND<dim> shft(typ - btype.ixType());
    bigend += shft;
    btype = IndexTypeND<dim>(typ);
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::convert (IndexTypeND<dim> t) noexcept
{
   for (int dir = 0; dir < dim; dir++)
   {
      const auto typ = t[dir];
      const auto bitval = btype[dir];
      const int off = typ - bitval;
      bigend.shift(dir,off);
      btype.setType(dir, static_cast<IndexType::CellIndex>(typ));
   }
   return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::surroundingNodes (int dir) noexcept
{
    if (!(btype[dir]))
    {
        bigend.shift(dir,1);
        //
        // Set dir'th bit to 1 = IndexType::NODE.
        //
        btype.set(dir);
    }
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::surroundingNodes () noexcept
{
    for (int i = 0; i < dim; ++i) {
        if ((btype[i] == 0)) {
            bigend.shift(i,1);
        }
    }
    btype.setall();
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::enclosedCells (int dir) noexcept
{
    if (btype[dir])
    {
        bigend.shift(dir,-1);
        //
        // Set dir'th bit to 0 = IndexType::CELL.
        //
        btype.unset(dir);
    }
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::enclosedCells () noexcept
{
    for (int i = 0 ; i < dim; ++i) {
        if (btype[i]) {
            bigend.shift(i,-1);
        }
    }
    btype.clear();
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Long
BoxND<dim>::index (const IntVectND<dim>& v) const noexcept
{
    IntVectND<dim> vz = v - smallend;
    Long result = vz[0];
    Long mult = length(0);
    for (int i = 1 ; i < dim; ++i) {
        result += mult * vz[i];
        mult *= length(i);
    }
    return result;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
BoxND<dim>::atOffset (Long offset) const noexcept
{
    IntVectND<dim> result = smallend;

    if constexpr (dim > 1) {
        GpuArray<Long, dim-1> mult{};
        mult[0] = length(0);
        for (int i = 1 ; i < dim-1; ++i) {
            mult[i] = mult[i-1] * length(i);
        }
        for (int i = dim-1 ; i > 0; --i) {
            Long idx = offset / mult[i-1];
            offset -= idx * mult[i-1];
            result[i] += static_cast<int>(idx);
        }
    }

    result[0] += static_cast<int>(offset);

    return result;
}

template<int dim>
template <int N, std::enable_if_t<( 1<=N && N<=3 ), int>>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
GpuArray<int,3>
BoxND<dim>::atOffset3d (Long offset) const noexcept
{
    Dim3 iv3d = atOffset(offset).dim3(0);
    return {{iv3d.x, iv3d.y, iv3d.z}};
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim>&
BoxND<dim>::setRange (int dir,
               int sm_index,
               int n_cells) noexcept
{
    smallend.setVal(dir,sm_index);
    bigend.setVal(dir,sm_index+n_cells-1);
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
BoxND<dim>::next (IntVectND<dim>& p) const noexcept // NOLINT(readability-convert-member-functions-to-static)
{
    BL_ASSERT(contains(p));

    ++p[0];

    for (int i = 0 ; i < dim-1; ++i) {
        if (p[i] > bigend[i]) {
            p[i] = smallend[i];
            ++p[i+1];
        } else {
            break;
        }
    }
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
BoxND<dim>::isSquare () const noexcept // NOLINT(readability-convert-member-functions-to-static)
{
    if constexpr (dim == 1) {
        return false; // can't build a square in 1-D
    } else {
        bool is_square = true;
        const IntVectND<dim>& sz = this->size();
        for (int i = 0 ; i < dim-1; ++i) {
            is_square = is_square && (sz[i] == sz[i+1]);
        }
        return is_square;
    }
}

//
// Modified BoxND is low end, returned BoxND is high end.
// If CELL: chop_pnt included in high end.
// If NODE: chop_pnt included in both Boxes.
//

template<int dim>
AMREX_GPU_HOST_DEVICE
inline
BoxND<dim>
BoxND<dim>::chop (int dir, int chop_pnt) noexcept
{
    //
    // Define new high end BoxND including chop_pnt.
    //
    IntVectND<dim> sm(smallend);
    IntVectND<dim> bg(bigend);
    sm.setVal(dir,chop_pnt);
    if (btype[dir])
    {
        //
        // NODE centered BoxND.
        //
        BL_ASSERT(chop_pnt > smallend[dir] && chop_pnt < bigend[dir]);
        //
        // Shrink original BoxND to just contain chop_pnt.
        //
        bigend.setVal(dir,chop_pnt);
    }
    else
    {
        //
        // CELL centered BoxND.
        //
        BL_ASSERT(chop_pnt > smallend[dir] && chop_pnt <= bigend[dir]);
        //
        // Shrink original BoxND to one below chop_pnt.
        //
        bigend.setVal(dir,chop_pnt-1);
    }
    return BoxND<dim>(sm,bg,btype);
}

template<int dim>
AMREX_GPU_HOST_DEVICE
inline
BoxND<dim>&
BoxND<dim>::shiftHalf (int dir, int num_halfs) noexcept
{
    const int nbit = (num_halfs<0 ? -num_halfs : num_halfs)%2;
    int nshift = num_halfs/2;
    //
    // Toggle btyp bit if num_halfs is odd.
    //
    const unsigned int bit_dir = btype[dir];
    if (nbit) {
        btype.flip(dir);
    }
    if (num_halfs < 0) {
        nshift -= (bit_dir ? nbit : 0);
    } else {
        nshift += (bit_dir ? 0 : nbit);
    }
    smallend.shift(dir,nshift);
    bigend.shift(dir,nshift);
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
inline
BoxND<dim>&
BoxND<dim>::shiftHalf (const IntVectND<dim>& iv) noexcept
{
    for (int i = 0; i < dim; i++) {
        shiftHalf(i,iv[i]);
    }
    return *this;
}

/// \cond DOXYGEN_IGNORE
class BoxCommHelper
{
public:

    explicit BoxCommHelper (const Box& bx, int* p_ = nullptr);

    [[nodiscard]] int* data () const& noexcept { return p; }
    int* data () && = delete;

    [[nodiscard]] Box make_box () const noexcept {
        return Box(IntVect(p), IntVect(p+AMREX_SPACEDIM), IntVect(p+2*AMREX_SPACEDIM));
    }

    [[nodiscard]] static int size () noexcept { return 3*AMREX_SPACEDIM; }

private:
    int* p;
    std::vector<int> v;
};
/// \endcond

class BoxConverter { // NOLINT
public:
    virtual Box doit (const Box& fine) const = 0; // NOLINT
    virtual BoxConverter* clone () const = 0; // NOLINT
    virtual ~BoxConverter () = default;
};

/**
 * \brief Gather boxes from all MPI ranks into \p bxs.
 *
 * The vector is resized to hold the union of every process' boxes.  If \p bxs
 * does not already have enough capacity, \p n_extra_reserve can be used to
 * reserve additional storage prior to the collective communication.
 *
 * \param bxs On input, the local boxes; on output, the concatenation of all ranks.
 * \param n_extra_reserve Additional slots to reserve before gathering.
 */
void AllGatherBoxes (Vector<Box>& bxs, int n_extra_reserve=0);

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b grown uniformly by \p i cells in every direction.
 *
 * Negative \p i values shrink the box by that many cells per direction.
 *
 * \param b Box to grow.
 * \param i Number of cells to add (negative shrinks).
 * \return A copy of \p b after the grow operation.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> grow (const BoxND<dim>& b, int i) noexcept
{
    BoxND<dim> result = b;
    result.grow(i);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to grow.
 * \param v Direction-wise cell counts to add (negative shrinks).
 * \return A copy of \p b after the grow operation.
*/
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> grow (const BoxND<dim>& b, const IntVectND<dim>& v) noexcept
{
    BoxND<dim> result = b;
    result.grow(v);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to grow.
 * \param idir Coordinate direction to grow.
 * \param n_cell Number of cells to add to the given direction (negative shrinks).
 * \return A copy of \p b after the grow operation.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> grow (const BoxND<dim>& b, int idir, int n_cell) noexcept
{
    BoxND<dim> result = b;
    result.grow(idir, n_cell);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> grow (const BoxND<dim>& b, Direction d, int n_cell) noexcept
{
    return grow(b, static_cast<int>(d), n_cell);
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b with the low end grown by \p n_cell cells in direction \p idir.
 *
 * Negative values shrink the low end.
 *
 * \param b Box to modify.
 * \param idir Coordinate direction to grow.
 * \param n_cell Number of cells to add (negative shrinks).
 * \return A copy of \p b after the asymmetric grow.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> growLo (const BoxND<dim>& b, int idir, int n_cell) noexcept
{
    BoxND<dim> result = b;
    result.growLo(idir, n_cell);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> growLo (const BoxND<dim>& b, Direction d, int n_cell) noexcept
{
    return growLo(b, static_cast<int>(d), n_cell);
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b with the high end grown by \p n_cell cells in direction \p idir.
 *
 * Negative values shrink the high end.
 *
 * \param b Box to modify.
 * \param idir Coordinate direction to grow.
 * \param n_cell Number of cells to add (negative shrinks).
 * \return A copy of \p b after the asymmetric grow.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> growHi (const BoxND<dim>& b, int idir, int n_cell) noexcept
{
    BoxND<dim> result = b;
    result.growHi(idir, n_cell);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> growHi (const BoxND<dim>& b, Direction d, int n_cell) noexcept
{
    return growHi(b, static_cast<int>(d), n_cell);
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b coarsened by the isotropic ratio \p ref_ratio.
 *
 * The coarsened box is constructed so that refining it by the same ratio
 * recovers the original box (taking nodal directions into account).
 *
 * \param b Box to coarsen.
 * \param ref_ratio Positive refinement ratio applied to every direction.
 * \return A copy of \p b after coarsening.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> coarsen (const BoxND<dim>& b, int ref_ratio) noexcept
{
    BoxND<dim> result = b;
    result.coarsen(IntVectND<dim>(ref_ratio));
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to coarsen.
 * \param ref_ratio Direction-wise positive refinement ratios.
 * \return A copy of \p b after coarsening.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> coarsen (const BoxND<dim>& b, const IntVectND<dim>& ref_ratio) noexcept
{
    BoxND<dim> result = b;
    result.coarsen(ref_ratio);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b refined by the isotropic ratio \p ref_ratio.
 *
 * Each direction is scaled by \p ref_ratio while respecting whether the box is
 * cell- or node-centered so that coarsening by the same ratio recovers the
 * original extent.
 *
 * \param b Box to refine.
 * \param ref_ratio Positive refinement ratio applied to every direction.
 * \return A copy of \p b after refinement.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> refine (const BoxND<dim>& b, int ref_ratio) noexcept
{
    BoxND<dim> result = b;
    result.refine(IntVectND<dim>(ref_ratio));
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to refine.
 * \param ref_ratio Direction-wise positive refinement ratios.
 * \return A copy of \p b after refinement.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> refine (const BoxND<dim>& b, const IntVectND<dim>& ref_ratio) noexcept
{
    BoxND<dim> result = b;
    result.refine(ref_ratio);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b shifted by \p nzones cells in direction \p dir.
 *
 * \param b Box to shift.
 * \param dir Coordinate direction to shift.
 * \param nzones Signed number of cells to shift by.
 * \return A copy of \p b after shifting.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> shift (const BoxND<dim>& b, int dir, int nzones) noexcept
{
    BoxND<dim> result = b;
    result.shift(dir, nzones);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to shift.
 * \param nzones Direction-wise signed cell offsets.
 * \return A copy of \p b after shifting.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> shift (const BoxND<dim>& b, const IntVectND<dim>& nzones) noexcept
{
    BoxND<dim> result = b;
    result.shift(nzones);
    return result;
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b expanded to nodal in direction \p dir.
 *
 * \param b Box to expand.
 * \param dir Coordinate direction promoted to node-centered indexing.
 * \return A copy of \p b with the requested direction nodal.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> surroundingNodes (const BoxND<dim>& b, int dir) noexcept
{
    BoxND<dim> bx(b);
    bx.surroundingNodes(dir);
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \overload
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> surroundingNodes (const BoxND<dim>& b, Direction d) noexcept
{
    return surroundingNodes(b, static_cast<int>(d));
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b converted to nodal indexing in every direction.
 *
 * \param b Box to expand.
 * \return A copy of \p b after promoting all directions to nodal.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> surroundingNodes (const BoxND<dim>& b) noexcept
{
    BoxND<dim> bx(b);
    bx.surroundingNodes();
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b converted to the nodal flags \p typ.
 *
 * \param b Box to convert.
 * \param typ Direction-wise nodal flags (0 cell, 1 node).
 * \return A copy of \p b after conversion.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> convert (const BoxND<dim>& b, const IntVectND<dim>& typ) noexcept
{
    BoxND<dim> bx(b);
    bx.convert(typ);
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \overload
 * \param b Box to convert.
 * \param typ Target index type descriptor.
 * \return A copy of \p b after conversion.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> convert (const BoxND<dim>& b, const IndexTypeND<dim>& typ) noexcept
{
    BoxND<dim> bx(b);
    bx.convert(typ);
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b converted to cell-centered in direction \p dir.
 *
 * \param b Box to convert.
 * \param dir Coordinate direction demoted to cell-centered indexing.
 * \return A copy of \p b with the requested direction cell-centered.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> enclosedCells (const BoxND<dim>& b, int dir) noexcept
{
    BoxND<dim> bx(b);
    bx.enclosedCells(dir);
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \overload
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> enclosedCells (const BoxND<dim>& b, Direction d) noexcept
{
    return enclosedCells(b, static_cast<int>(d));
}

/**
 * \ingroup amrex_index_space
 * \brief Return a copy of \p b converted to cell-centered in every direction.
 *
 * \param b Box to convert.
 * \return A copy of \p b after demoting all directions to cell-centered.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> enclosedCells (const BoxND<dim>& b) noexcept
{
    BoxND<dim> bx(b);
    bx.enclosedCells();
    return bx;
}

/**
 * \ingroup amrex_index_space
 * \brief Return the nodal boundary region adjacent to the low side of \p b.
 *
 * \param b Reference box.
 * \param dir Coordinate direction of the face.
 * \param len Number of nodal layers to include (positive).
 * \return A nodal box covering the low boundary.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> bdryLo (const BoxND<dim>& b, int dir, int len=1) noexcept
{
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    int sm = low[dir];
    low.setVal(dir,sm-len+1);
    hi.setVal(dir,sm);
    //
    // set dir'th bit to 1 = IndexType::NODE.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.set(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the nodal boundary region adjacent to the high side of \p b.
 *
 * \param b Reference box.
 * \param dir Coordinate direction of the face.
 * \param len Number of nodal layers to include (positive).
 * \return A nodal box covering the high boundary.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> bdryHi (const BoxND<dim>& b, int dir, int len=1) noexcept
{
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    auto const bitval = b.type()[dir];
    int bg = hi[dir]  + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg+len-1);
    //
    // Set dir'th bit to 1 = IndexType::NODE.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.set(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the nodal boundary region adjacent to \p face of \p b.
 *
 * \param b Reference box.
 * \param face Orientation of the boundary.
 * \param len Number of nodal layers to include (positive).
 * \return A nodal box covering the requested boundary.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> bdryNode (const BoxND<dim>& b, Orientation face, int len=1) noexcept
{
    int dir = face.coordDir();
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    if (face.isLow())
    {
        int sm = low[dir];
        low.setVal(dir,sm-len+1);
        hi.setVal(dir,sm);
    }
    else
    {
        int bitval = b.type()[dir];
        int bg = hi[dir]  + 1 - bitval%2;
        low.setVal(dir,bg);
        hi.setVal(dir,bg+len-1);
    }
    //
    // Set dir'th bit to 1 = IndexType::NODE.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.set(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the cell-centered region adjacent to the low side of \p b.
 *
 * \param b Reference box.
 * \param dir Coordinate direction of the face.
 * \param len Number of cell layers to include (positive).
 * \return A cell-centered box disjoint from \p b on the low side.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> adjCellLo (const BoxND<dim>& b, int dir, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    int sm = low[dir];
    low.setVal(dir,sm - len);
    hi.setVal(dir,sm - 1);
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.unset(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the cell-centered region adjacent to the high side of \p b.
 *
 * \param b Reference box.
 * \param dir Coordinate direction of the face.
 * \param len Number of cell layers to include (positive).
 * \return A cell-centered box disjoint from \p b on the high side.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> adjCellHi (const BoxND<dim>& b, int dir, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    int bitval = b.type()[dir];
    int bg = hi[dir] + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg + len - 1);
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.unset(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the cell-centered region adjacent to \p face of \p b.
 *
 * \param b Reference box.
 * \param face Orientation of the boundary.
 * \param len Number of cell layers to include (positive).
 * \return A cell-centered box disjoint from \p b on the requested side.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> adjCell (const BoxND<dim>& b, Orientation face, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVectND<dim> low(b.smallEnd());
    IntVectND<dim> hi(b.bigEnd());
    int dir = face.coordDir();
    if (face.isLow())
    {
        int sm = low[dir];
        low.setVal(dir,sm - len);
        hi.setVal(dir,sm - 1);
    }
    else
    {
        int bitval = b.type()[dir];
        int bg = hi[dir] + 1 - bitval%2;
        low.setVal(dir,bg);
        hi.setVal(dir,bg + len - 1);
    }
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexTypeND<dim> typ(b.ixType());
    typ.unset(dir);
    return BoxND<dim>(low,hi,typ);
}

/**
 * \ingroup amrex_index_space
 * \brief Return the minimal BoxND containing both \p b1 and \p b2.
 *
 * The boxes must share the same index type.
 *
 * \param b1 First box.
 * \param b2 Second box.
 * \return The minimal box that contains both inputs.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> minBox (const BoxND<dim>& b1, const BoxND<dim>& b2) noexcept
{
    BoxND<dim> result = b1;
    result.minBox(b2);
    return result;
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    std::ostream& box_write (std::ostream& os, const int * smallend, const int * bigend,
                             const int * type, int dim);
    std::istream& box_read (std::istream& is, int * smallend, int * bigend, int * type, int dim);

    template<std::size_t...Ns, class T, class U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    auto BoxSplit_imp (std::index_sequence<Ns...>,
                       const T& lo, const T& hi, const U& typ) noexcept {
        return makeTuple(BoxND(get<Ns>(lo), get<Ns>(hi), get<Ns>(typ))...);
    }
}
/// \endcond

/**
 * \brief Write an ASCII representation of \p bx to \p os.
 *
 * \param os Output stream.
 * \param bx Box to print.
 * \return The output stream.
 */
template<int dim>
std::ostream& operator<< (std::ostream& os, const BoxND<dim>& bx)
{
    IntVectND<dim> type = bx.type();
    return detail::box_write(os, bx.smallEnd().begin(), bx.bigEnd().begin(), type.begin(), dim);
}

/**
 * \brief Read a BoxND from \p is.
 *
 * \param is Input stream.
 * \param bx Box filled by the read.
 * \return The input stream.
 */
template<int dim>
std::istream& operator>> (std::istream& is, BoxND<dim>& bx) {
    IntVectND<dim> small;
    IntVectND<dim> big;
    IntVectND<dim> type;
    detail::box_read(is, small.begin(), big.begin(), type.begin(), dim);
    bx = BoxND<dim>{small, big, type};
    return is;
}

/**
 * \brief Concatenate boxes along their dimensions.
 *
 * The returned box has dimension equal to the sum of the input dimensions.
 *
 * \param bx First box.
 * \param boxes Remaining boxes.
 * \return The concatenated box.
 */
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr BoxND<detail::get_sum<d, dims...>()>
BoxCat (const BoxND<d>& bx, const BoxND<dims>&...boxes) noexcept {
    auto lo = IntVectCat(bx.smallEnd(), boxes.smallEnd()...);
    auto hi = IntVectCat(bx.bigEnd(), boxes.bigEnd()...);
    auto typ = IndexTypeCat(bx.ixType(), boxes.ixType()...);
    return BoxND<detail::get_sum<d, dims...>()>{lo, hi, typ};
}

/**
 * \brief Split \p bx into boxes whose dimensions match the template arguments.
 *
 * \param bx Source box to split.
 * \return A tuple containing the split boxes.
 */
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr GpuTuple<BoxND<d>, BoxND<dims>...>
BoxSplit (const BoxND<detail::get_sum<d, dims...>()>& bx) noexcept {
    auto lo = IntVectSplit<d, dims...>(bx.smallEnd());
    auto hi = IntVectSplit<d, dims...>(bx.bigEnd());
    auto typ = IndexTypeSplit<d, dims...>(bx.ixType());
    return detail::BoxSplit_imp(std::make_index_sequence<1 + sizeof...(dims)>(), lo, hi, typ);
}

/**
 * \brief Shrink \p bx to its first \p new_dim directions.
 *
 * \param bx Box to shrink.
 * \return A lower-dimensional view of \p bx.
 */
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr BoxND<new_dim>
BoxShrink (const BoxND<old_dim>& bx) noexcept {
    return bx.template shrink<new_dim>();
}

/**
 * \brief Expand \p bx to \p new_dim directions, padding with zeros.
 *
 * Additional dimensions are initialized to cell-centered [0,0].
 *
 * \param bx Box to expand.
 * \return A higher-dimensional box containing \p bx.
 */
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr BoxND<new_dim>
BoxExpand (const BoxND<old_dim>& bx) noexcept {
    return bx.template expand<new_dim>();
}

/**
 * \brief Resize \p bx to \p new_dim directions via shrink or expand.
 *
 * \param bx Box to resize.
 * \return The resized box.
 */
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr BoxND<new_dim>
BoxResize (const BoxND<old_dim>& bx) noexcept {
    return bx.template resize<new_dim>();
}

/**
 * \brief Return the inclusive lower corner of \p box as an IntVectND.
 *
 * \param box Box to query.
 * \return The lower corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> lbound_iv (BoxND<dim> const& box) noexcept
{
    return box.smallEnd();
}

/**
 * \brief Return the inclusive upper corner of \p box as an IntVectND.
 *
 * \param box Box to query.
 * \return The upper corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> ubound_iv (BoxND<dim> const& box) noexcept
{
    return box.bigEnd();
}

/**
 * \brief Return the iterator begin IntVectND for \p box.
 *
 * \param box Box to query.
 * \return The first valid index.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> begin_iv (BoxND<dim> const& box) noexcept
{
    return box.smallEnd();
}

/**
 * \brief Return the iterator end IntVectND for \p box.
 *
 * The end corresponds to one past the last valid index.
 *
 * \param box Box to query.
 * \return One-past-the-end index.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> end_iv (BoxND<dim> const& box) noexcept
{
    return box.bigEnd() + 1;
}

/**
 * \brief Return the side lengths of \p box as an IntVectND.
 *
 * \param box Box to query.
 * \return Direction-wise lengths.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> length_iv (BoxND<dim> const& box) noexcept
{
    return box.bigEnd() - box.smallEnd() + 1;
}

/**
 * \brief Return the component-wise maximum of the lower corners of two boxes.
 *
 * \param b1 First box.
 * \param b2 Second box.
 * \return The maximum lower corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> max_lbound_iv (BoxND<dim> const& b1, BoxND<dim> const& b2) noexcept
{
    return max(b1.smallEnd(), b2.smallEnd());
}

/**
 * \brief Return the component-wise maximum of the lower corner of \p b1 and \p lo.
 *
 * \param b1 Box to query.
 * \param lo Lower corner candidate.
 * \return The maximum lower corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> max_lbound_iv (BoxND<dim> const& b1, IntVectND<dim> const& lo) noexcept
{
    return max(b1.smallEnd(), lo);
}

/**
 * \brief Return the component-wise minimum of the upper corners of two boxes.
 *
 * \param b1 First box.
 * \param b2 Second box.
 * \return The minimum upper corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> min_ubound_iv (BoxND<dim> const& b1, BoxND<dim> const& b2) noexcept
{
    return min(b1.bigEnd(), b2.bigEnd());
}

/**
 * \brief Return the component-wise minimum of the upper corner of \p b1 and \p hi.
 *
 * \param b1 Box to query.
 * \param hi Upper corner candidate.
 * \return The minimum upper corner.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> min_ubound_iv (BoxND<dim> const& b1, IntVectND<dim> const& hi) noexcept
{
    return min(b1.bigEnd(), hi);
}

/**
 * \brief Return the inclusive lower corner of \p box padded to Dim3.
 *
 * \param box Box to query.
 * \return The lower corner as Dim3.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 lbound (BoxND<dim> const& box) noexcept
{
    return box.smallEnd().dim3();
}

/**
 * \brief Return the inclusive upper corner of \p box padded to Dim3.
 *
 * \param box Box to query.
 * \return The upper corner as Dim3.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 ubound (BoxND<dim> const& box) noexcept
{
    return box.bigEnd().dim3();
}

/**
 * \brief Return the iterator begin coordinate of \p box as Dim3.
 *
 * \param box Box to query.
 * \return The first valid index.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 begin (BoxND<dim> const& box) noexcept
{
    return box.smallEnd().dim3();
}

/**
 * \brief Return the iterator end coordinate of \p box as Dim3.
 *
 * The end coordinate is one past the last valid index.
 *
 * \param box Box to query.
 * \return One-past-the-end index.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 end (BoxND<dim> const& box) noexcept
{
    return (box.bigEnd() + 1).dim3(1);
}

/**
 * \brief Return the side lengths of \p box as Dim3.
 *
 * \param box Box to query.
 * \return Direction-wise lengths.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 length (BoxND<dim> const& box) noexcept
{
    return (box.bigEnd() - box.smallEnd() + 1).dim3(1);
}

// Max of lower bound
/**
 * \brief Return the Dim3 maximum of the lower corners of two boxes.
 *
 * \param b1 First box.
 * \param b2 Second box.
 * \return Component-wise maximum lower corner.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 max_lbound (BoxND<dim> const& b1, BoxND<dim> const& b2) noexcept
{
    return max(b1.smallEnd(), b2.smallEnd()).dim3();
}

/**
 * \brief Return the Dim3 maximum of the lower corner of \p b1 and \p lo.
 *
 * \param b1 Box to query.
 * \param lo Candidate lower corner.
 * \return Component-wise maximum lower corner.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 max_lbound (BoxND<dim> const& b1, Dim3 const& lo) noexcept
{
    return max(b1.smallEnd(), IntVectND<dim>(lo)).dim3();
}

// Min of upper bound
/**
 * \brief Return the Dim3 minimum of the upper corners of two boxes.
 *
 * \param b1 First box.
 * \param b2 Second box.
 * \return Component-wise minimum upper corner.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 min_ubound (BoxND<dim> const& b1, BoxND<dim> const& b2) noexcept
{
    return min(b1.bigEnd(), b2.bigEnd()).dim3();
}

/**
 * \brief Return the Dim3 minimum of the upper corner of \p b1 and \p hi.
 *
 * \param b1 Box to query.
 * \param hi Candidate upper corner.
 * \return Component-wise minimum upper corner.
 */
template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 min_ubound (BoxND<dim> const& b1, Dim3 const& hi) noexcept
{
    return min(b1.bigEnd(), IntVectND<dim>(hi)).dim3();
}

/**
 * \ingroup amrex_index_space
 * \brief Return the index-space extent of \p b1.
 *
 * \param b1 Box to query.
 * \return The same box, treated as an index bound helper.
 */
template<int dim>
[[nodiscard]]
AMREX_FORCE_INLINE
BoxND<dim> getIndexBounds (BoxND<dim> const& b1) noexcept
{
    return b1;
}

/**
 * \ingroup amrex_index_space
 * \brief Return the minimal index-space box covering \p b1 and \p b2.
 *
 * The result adopts the index type of \p b2.
 *
 * \param b1 First box.
 * \param b2 Second box (also supplies the target index type).
 * \return A box covering both inputs.
 */
template<int dim>
[[nodiscard]]
AMREX_FORCE_INLINE
BoxND<dim> getIndexBounds (BoxND<dim> const& b1, BoxND<dim> const& b2) noexcept
{
    BoxND<dim> b = b1;
    b.setType(b2.ixType());
    b.minBox(b2);
    return b;
}

/**
 * \brief Return the minimal index-space box covering all inputs.
 *
 * \param b1 First box-like object.
 * \param b2 Second box-like object.
 * \param b3 Additional box-like objects.
 * \return A box covering every input.
 */
template <class T, class ... Ts>
[[nodiscard]]
AMREX_FORCE_INLINE
auto getIndexBounds (T const& b1, T const& b2, Ts const& ... b3) noexcept
{
    return getIndexBounds(getIndexBounds(b1,b2),b3...);
}


/**
 * \brief Return the cell coordinate corresponding to the flattened index \p icell.
 *
 * \param boxes Pointer to an array of boxes.
 * \param nboxes Number of boxes in the array.
 * \param icell Zero-based flattened cell index across all boxes.
 * \return The IntVectND coordinate of the requested cell.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> getCell (BoxND<dim> const* boxes, int nboxes, Long icell) noexcept
{
    int ibox;
    Long ncells_subtotal = 0;
    Long offset = 0;
    for (ibox = 0; ibox < nboxes; ++ibox) {
        const Long n = boxes[ibox].numPts();
        ncells_subtotal += n;
        if (icell < ncells_subtotal) {
            offset = icell - (ncells_subtotal - n);
            break;
        }
    }
    return boxes[ibox].atOffset(offset);
}

/**
 * \brief Return a copy of \p b collapsed to a slab perpendicular to \p direction.
 *
 * \param b Box to slice.
 * \param direction Coordinate direction of the slab.
 * \param slab_index Coordinate value at which to slice.
 * \return A slab box extracted from \p b.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> makeSlab (BoxND<dim> const& b, int direction, int slab_index) noexcept
{
    BoxND<dim> r = b;
    r.makeSlab(direction,slab_index);
    return r;
}

/**
 * \ingroup amrex_index_space
 * \brief Create a single-cell box at coordinates (\p i,\p j,\p k).
 *
 * \param i Cell index along x.
 * \param j Cell index along y.
 * \param k Cell index along z.
 * \param typ Target index type (cell by default).
 * \return A single-cell box.
 */
template<int dim=AMREX_SPACEDIM, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> makeSingleCellBox (int i, int j, int k, IndexTypeND<dim> typ = IndexTypeND<dim>::TheCellType())
{
    Dim3 p3d{i, j, k};
    IntVectND<dim> vect{p3d};
    return BoxND<dim>{vect, vect, typ};
}

/**
 * \ingroup amrex_index_space
 * \brief Create a single-cell box at \p vect.
 *
 * \param vect Cell indices of the single cell.
 * \param typ Target index type (cell by default).
 * \return A single-cell box.
 */
template<int dim>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
BoxND<dim> makeSingleCellBox (IntVectND<dim> const& vect, IndexTypeND<dim> typ = IndexTypeND<dim>::TheCellType())
{
    return BoxND<dim>{vect, vect, typ};
}

/**
 * \brief Utility that maps flattened cell indices back to IntVectND coordinates.
 *
 * The constructor captures a sequence of \ref Math::FastDivmodU64 helpers that
 * implement fast division by the running products of box lengths, so repeated
 * calls to \ref intVect or the Dim3 functor avoid costly integer divisions while
 * walking a \ref BoxND.
 *
 * \tparam dim Box dimensionality.
 */
template<int dim>
struct BoxIndexerND
{
    std::uint64_t npts;
    Math::FastDivmodU64 fdm[dim-1];
    IntVectND<dim> lo;

    BoxIndexerND (BoxND<dim> const& box)
        : npts(box.numPts()),
          lo  (box.smallEnd())
        {
            std::uint64_t mult = 1;
            for (int i=0; i<dim-1; ++i) {
                mult *= box.length(i);
                fdm[i] = Math::FastDivmodU64(mult);
            }
        }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Convert flattened cell index \p icell to its IntVectND coordinate.
     *
     * \param icell Zero-based flattened cell index.
     * \return Coordinate within the tracked box.
     */
    IntVectND<dim> intVect (std::uint64_t icell) const
    {
        IntVectND<dim> retval = lo;

        for (int i=dim-1; i>0; --i) {
            std::uint64_t quotient, remainder;
            fdm[i-1](quotient, remainder, icell);
            retval[i] += quotient;
            icell = remainder;
        }

        retval[0] += icell;

        return retval;
    }

    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Functor returning the coordinate padded to Dim3.
     *
     * \param icell Zero-based flattened cell index.
     * \return Coordinate as Dim3.
     */
    Dim3 operator() (std::uint64_t icell) const
    {
        return intVect(icell).dim3();
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Return the number of cells covered by the indexed box.
     *
     * \return Total number of points.
     */
    std::uint64_t numPts () const { return npts; }
};

/**
 * \brief Specialization for the 1D case where no fast-division table is needed.
 */
template<>
struct BoxIndexerND<1>
{
    std::uint64_t npts;

    int lo;

    /**
     * \brief Construct the 1D indexer; no fast divisors are needed because only one axis exists.
     *
     * \param box Box that supplies the span.
     */
    BoxIndexerND (BoxND<1> const& box)
        : npts(box.numPts()),
          lo(box.smallEnd(0))
        {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Convert flattened cell index \p icell to the 1D coordinate.
     *
     * \param icell Zero-based flattened cell index.
     * \return Coordinate within the tracked box.
     */
    IntVectND<1> intVect (std::uint64_t icell) const
    {
        return IntVectND<1>{int(icell)+lo};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Return the coordinate padded to Dim3.
     *
     * \param icell Zero-based flattened cell index.
     * \return Coordinate as Dim3.
     */
    Dim3 operator() (std::uint64_t icell) const
    {
        return {int(icell)+lo, 0, 0};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    /**
     * \brief Return the number of cells covered by the indexed box.
     *
     * \return Total number of points.
     */
    std::uint64_t numPts () const { return npts; }
};

using BoxIndexer = BoxIndexerND<AMREX_SPACEDIM>;

}

#endif /*AMREX_BOX_H*/
