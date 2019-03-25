
#ifndef AMREX_BOX_H_
#define AMREX_BOX_H_

#include <iosfwd>
#include <array>

#include <AMReX_ArrayLim.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Orientation.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>
#include <AMReX_CudaMemory.H>

namespace amrex
{
class BoxCommHelper;

/**
* \brief A Rectangular Domain on an Integer Lattice
*
* A Box is an abstraction for defining discrete regions of
* SPACEDIM indexing space.  Boxes have an IndexType, which defines
* IndexType::CELL or IndexType::NODE based points for each direction
* and a low and high INTVECT which defines the lower and upper corners
* of the Box.  Boxes can exist in positive and negative indexing space.
*
* Box is a dimension dependent class, so SPACEDIM must be
* defined as either 1, 2, or 3 when compiling.
*/
class Box // : public Gpu::Managed
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<Box>::type();
    friend class BoxCommHelper;

public:
    /*
    * \brief The default constructor.  For safety, the constructed Box is
    * invalid and may be tested for validity with ok().
    * DO NOT CHANGE THIS BEHAVIOR!
    */
    AMREX_GPU_HOST_DEVICE
    constexpr Box () noexcept
        : smallend(1),
          bigend(0),
          btype()
        {}

    //! Construct cell-centered type Box.
    AMREX_GPU_HOST_DEVICE
    constexpr Box (const IntVect& small, const IntVect& big) noexcept
        : smallend(small),
          bigend(big)
        {}

    //! Construct box with specified lengths.
    AMREX_GPU_HOST_DEVICE
    Box (const IntVect& small, const int* vec_len) noexcept
        : smallend(small),
          bigend(AMREX_D_DECL(small[0]+vec_len[0]-1,
                              small[1]+vec_len[1]-1,
                              small[2]+vec_len[2]-1))
        {}

    /**
    * \brief Construct Box with given type.  small and big are expected
    * to be consistent with given type.
    */
    AMREX_GPU_HOST_DEVICE
    Box (const IntVect& small, const IntVect& big, const IntVect& typ) noexcept
        : smallend(small),
          bigend(big),
          btype(typ)
        {
            BL_ASSERT(typ.allGE(IntVect::TheZeroVector()) && typ.allLE(IntVect::TheUnitVector()));
        }

    //! Construct dimension specific Boxes.
    AMREX_GPU_HOST_DEVICE
    Box (const IntVect& small, const IntVect& big, IndexType t) noexcept
        : smallend(small),
          bigend(big),
          btype(t)
        {}

    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    //! Get the smallend of the box.
    AMREX_GPU_HOST_DEVICE
    const IntVect& smallEnd () const& noexcept { return smallend; }
    const IntVect& smallEnd () && = delete;

    //! Returns the coordinate of the low end in the given direction.
    AMREX_GPU_HOST_DEVICE
    int smallEnd (int dir) const& noexcept { return smallend[dir]; }

    //! Get the bigend.
    AMREX_GPU_HOST_DEVICE
    const IntVect& bigEnd () const& noexcept { return bigend; }
    const IntVect& bigEnd () && = delete;

    //! Returns the coordinate of the high end in the given direction.
    AMREX_GPU_HOST_DEVICE
    int bigEnd (int dir) const noexcept { return bigend[dir]; }

    //! Returns the indexing type.
    AMREX_GPU_HOST_DEVICE
    IndexType ixType () const noexcept { return btype; }

    //! Returns the indexing type.
    AMREX_GPU_HOST_DEVICE
    IntVect type () const noexcept { return btype.ixType(); }

    //! Returns the indexing type in the specified direction.
    AMREX_GPU_HOST_DEVICE
    IndexType::CellIndex type (int dir) const noexcept { return btype.ixType(dir); }

    //! Return the length of the Box.
    AMREX_GPU_HOST_DEVICE
    IntVect size () const noexcept
    {
        return IntVect(AMREX_D_DECL(bigend[0]-smallend[0] + 1,
                                    bigend[1]-smallend[1] + 1,
                                    bigend[2]-smallend[2] + 1));
    }

    //! Return the length of the Box.
    AMREX_GPU_HOST_DEVICE
    IntVect length () const noexcept
    {
        return IntVect(AMREX_D_DECL(bigend[0]-smallend[0] + 1,
                                    bigend[1]-smallend[1] + 1,
                                    bigend[2]-smallend[2] + 1));
    }

    //! Return the length of the Box in given direction.
    AMREX_GPU_HOST_DEVICE
    int length (int dir) const noexcept { return bigend[dir] - smallend[dir] + 1; }

    AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> length3d () const noexcept {
#if (AMREX_SPACEDIM == 1)
        return {bigend[0]-smallend[0]+1, 1, 1};
#elif (AMREX_SPACEDIM == 2)
        return {bigend[0]-smallend[0]+1, bigend[1]-smallend[1]+1, 1};
#elif (AMREX_SPACEDIM == 3)
        return {bigend[0]-smallend[0]+1, bigend[1]-smallend[1]+1, bigend[2]-smallend[2]+1};
#endif
    }

    AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> loVect3d () const noexcept {
#if (AMREX_SPACEDIM == 1)
        return {smallend[0], 0, 0};
#elif (AMREX_SPACEDIM == 2)
        return {smallend[0], smallend[1], 0};
#elif (AMREX_SPACEDIM == 3)
        return {smallend[0], smallend[1], smallend[2]};
#endif
    }

    AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> hiVect3d () const noexcept {
#if (AMREX_SPACEDIM == 1)
        return {bigend[0], 0, 0};
#elif (AMREX_SPACEDIM == 2)
        return {bigend[0], bigend[1], 0};
#elif (AMREX_SPACEDIM == 3)
        return {bigend[0], bigend[1], bigend[2]};
#endif
    }

    //! Returns a constant pointer the array of low end coordinates. Useful for calls to FORTRAN.
    AMREX_GPU_HOST_DEVICE
    const int* loVect () const& noexcept { return smallend.getVect(); }
    AMREX_GPU_HOST_DEVICE
    const int* loVect () && = delete;
    //! Returns a constant pointer the array of high end coordinates. Useful for calls to FORTRAN.
    AMREX_GPU_HOST_DEVICE
    const int* hiVect () const& noexcept { return bigend.getVect(); }
    AMREX_GPU_HOST_DEVICE
    const int* hiVect () && = delete;

    //! Returns the coordinate normal to given face.
    AMREX_GPU_HOST_DEVICE
    int operator[] (Orientation face) const noexcept {
        const int dir = face.coordDir();
        return face.isLow() ? smallend[dir] : bigend[dir];
    }

    //! Checks if it is an empty box.
    AMREX_GPU_HOST_DEVICE
    bool isEmpty () const noexcept { return !ok(); }

    //! Checks if it is a proper Box (including a valid type).
    AMREX_GPU_HOST_DEVICE
    bool ok () const noexcept { return bigend.allGE(smallend) && btype.ok(); }

    //! Returns true if argument is contained within Box.
    AMREX_GPU_HOST_DEVICE
    bool contains (const IntVect& p) const noexcept { return p.allGE(smallend) && p.allLE(bigend); }

    /** \brief Returns true if argument is contained within Box.
    * It is an error if the Boxes have different types.
    */
    AMREX_GPU_HOST_DEVICE
    bool contains (const Box& b) const noexcept
    {
	BL_ASSERT(sameType(b));
	return b.smallend.allGE(smallend) && b.bigend.allLE(bigend);
    }

    //! Returns true if argument is strictly contained within Box.
    AMREX_GPU_HOST_DEVICE
    bool strictly_contains (const IntVect& p) const noexcept { return p.allGT(smallend) && p.allLT(bigend); }

    /**
    * \brief Returns true if argument is strictly contained within Box.
    * It is an error if the Boxes have different types.
    */
    AMREX_GPU_HOST_DEVICE
    bool strictly_contains (const Box& b) const noexcept
    {
	BL_ASSERT(sameType(b));
	return b.smallend.allGT(smallend) && b.bigend.allLT(bigend);
    }
    /**
    * \brief Returns true if Boxes have non-null intersections.
    * It is an error if the Boxes have different types.
    */
    AMREX_GPU_HOST_DEVICE
    bool intersects (const Box& b) const noexcept { Box isect(*this); isect &= b; return isect.ok(); }

    /**
    * \brief Returns true is Boxes same size, ie translates of each other,.
    * It is an error if they have different types.
    */
    AMREX_GPU_HOST_DEVICE
    bool sameSize (const Box& b) const noexcept {
        BL_ASSERT(sameType(b));
        return AMREX_D_TERM(length(0) == b.length(0),
                         && length(1) == b.length(1),
                         && length(2) == b.length(2));
    }

    //! Returns true if Boxes have same type.
    AMREX_GPU_HOST_DEVICE
    bool sameType (const Box &b) const noexcept { return btype == b.btype; }

    //! Returns true if Boxes are identical (including type).
    AMREX_GPU_HOST_DEVICE
    bool operator== (const Box& b) const noexcept { return smallend == b.smallend && bigend == b.bigend && b.btype == btype; }

    //! Returns true if Boxes differ (including type).
    AMREX_GPU_HOST_DEVICE
    bool operator!= (const Box& b) const noexcept { return !operator==(b); }

    AMREX_GPU_HOST_DEVICE
    bool operator< (const Box& rhs) const noexcept
    {
        return btype < rhs.btype ||
             ((btype == rhs.btype) &&
	         (  (smallend < rhs.smallend) ||
                   ((smallend == rhs.smallend) && (bigend < rhs.bigend)) ));
    }
    AMREX_GPU_HOST_DEVICE
    bool operator <= (const Box& rhs) const noexcept {
        return !(rhs < *this);
    }
    AMREX_GPU_HOST_DEVICE
    bool operator> (const Box& rhs) const noexcept {
        return rhs < *this;
    }
    AMREX_GPU_HOST_DEVICE
    bool operator>= (const Box& rhs) const noexcept {
        return !(*this < rhs);
    }

    //! Returns true if Box is cell-centered in all indexing directions.
    AMREX_GPU_HOST_DEVICE
    bool cellCentered () const noexcept { return !btype.any(); }

    /**
    * \brief Returns the number of points contained in the Box, else
    * abort()s if the number cannot be represented in a long.
    */
    AMREX_GPU_HOST_DEVICE
    long numPts () const noexcept {
        auto r = AMREX_D_TERM( static_cast<unsigned long>(length(0)),
                              *static_cast<unsigned long>(length(1)),
                              *static_cast<unsigned long>(length(2)));
        return static_cast<long>(r);
    }

    /**
    * \brief Returns the number of points contained in the Box.
    * This is intended for use only in diagnostic messages.
    */
    AMREX_GPU_HOST_DEVICE
    double d_numPts () const noexcept {
        BL_ASSERT(ok());
        return AMREX_D_TERM(double(length(0)), *double(length(1)), *double(length(2)));
    }

    /**
    * \brief Return the volume, in indexing space, of region enclosed by
    * this Box.   This is identical to numPts() for CELL centered
    * Box; otherwise, numPts() > volume().  Abort()s if the number
    * cannot be represented in a long.
    */
    AMREX_GPU_HOST_DEVICE
    long volume () const noexcept {
        return AMREX_D_TERM( static_cast<long>(length(0)-btype[0]),
                            *static_cast<long>(length(1)-btype[1]),
                            *static_cast<long>(length(2)-btype[2]));
    }

    /**
    * \brief Returns length of longest side.  dir is modified to give
    * direction with longest side: 0...SPACEDIM-1. Ignores type.
    */
    AMREX_GPU_HOST_DEVICE
    int longside (int& dir) const noexcept {
        int maxlen = length(0);
        dir = 0;
        for (int i = 1; i < AMREX_SPACEDIM; i++)
        {
            if (length(i) > maxlen)
            {
                maxlen = length(i);
                dir = i;
            }
        }
        return maxlen;
    }

    //! Returns length of longest side.  Ignores type.
    AMREX_GPU_HOST_DEVICE
    int longside () const noexcept {
        int ignore = 0;
        return longside(ignore);
    }

    /**
    * \brief Returns length of shortest side.  dir is modified to give
    *  direction with shortest side: 0...SPACEDIM-1.  Ignores type.
    */
    AMREX_GPU_HOST_DEVICE
    int shortside (int& dir) const noexcept {
        int minlen = length(0);
        dir = 0;
        for (int i = 1; i < AMREX_SPACEDIM; i++)
        {
            if (length(i) < minlen)
            {
                minlen = length(i);
                dir = i;
            }
        }
        return minlen;
    }

    //! Returns length of shortest side.  Ignores type.
    AMREX_GPU_HOST_DEVICE
    int shortside () const noexcept {
        int ignore = 0;
        return shortside(ignore);
    }

    /**
    * \brief Returns offset of point from smallend; i.e.
    * index(smallend) -> 0, bigend would return numPts()-1.
    * Is used in accessing FArrayBox.
    */
    AMREX_GPU_HOST_DEVICE
    long index (const IntVect& v) const noexcept;

    //! Given the offset, compute IntVect
    AMREX_GPU_HOST_DEVICE
    IntVect atOffset (long offset) const noexcept;

    AMREX_GPU_HOST_DEVICE
    GpuArray<int,3> atOffset3d (long offset) const noexcept;

    //! Redefine the small end of the Box.
    AMREX_GPU_HOST_DEVICE
    Box& setSmall (const IntVect& sm) noexcept { smallend = sm; return *this; }

    //! Redefine the small end of the Box.
    AMREX_GPU_HOST_DEVICE
    Box& setSmall (int dir, int sm_index) noexcept { smallend.setVal(dir,sm_index); return *this; }

    //! Redefine the big end of the Box.
    AMREX_GPU_HOST_DEVICE
    Box& setBig (const IntVect& bg) noexcept { bigend = bg; return *this; }

    //! Redefine the big end of the Box.
    AMREX_GPU_HOST_DEVICE
    Box& setBig (int dir, int bg_index) noexcept { bigend.setVal(dir,bg_index); return *this; }

    /**
    * \brief Set the entire range in a given direction, starting at
    * sm_index with length n_cells.  NOTE: This will yield an
    * illegal Box if n_cells <= 0.
    */
    AMREX_GPU_HOST_DEVICE
    Box& setRange (int dir,
                   int sm_index,
                   int n_cells = 1) noexcept;

    //! Set indexing type
    AMREX_GPU_HOST_DEVICE
    Box& setType (const IndexType& t) noexcept { btype = t; return *this; }

    //! Shift this Box nzones indexing positions in coordinate direction dir.
    AMREX_GPU_HOST_DEVICE
    Box& shift (int dir, int nzones) noexcept { smallend.shift(dir,nzones); bigend.shift(dir,nzones); return *this; }

    //! Equivalent to b.shift(0,iv[0]).shift(1,iv[1]) ....
    AMREX_GPU_HOST_DEVICE
    Box& shift (const IntVect& iv) noexcept { smallend.shift(iv); bigend.shift(iv); return *this; }

    /**
    * \brief This member shifts the Box by "half" indicies, thereby
    * converting the Box from type CELL to NODE and visa-versa.
    * b.shiftHalf(0,1)  shifts b to the right by 1/2 cells.
    * b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
    * NOTE: If num_halfs is EVEN the shift is num_halfs/2 full
    * zones and hence will not change the type.
    * This is: b.shifthalf(4) == b.shift(2).
    */
    AMREX_GPU_HOST_DEVICE
    Box& shiftHalf (int dir, int num_halfs) noexcept;

    //! Equivalent to b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...
    AMREX_GPU_HOST_DEVICE
    Box& shiftHalf (const IntVect& iv) noexcept;

    /**
    * \brief Convert the Box from the current type into the
    * argument type.  This may change the Box coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    AMREX_GPU_HOST_DEVICE
    Box& convert (IndexType typ) noexcept;

    /**
    * \brief Convert the Box from the current type into the
    * argument type.  This may change the Box coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    AMREX_GPU_HOST_DEVICE
    Box& convert (const IntVect& typ) noexcept;

    //! Convert to NODE type in all directions.
    AMREX_GPU_HOST_DEVICE
    Box& surroundingNodes () noexcept;

    //! Convert to NODE type in given direction.
    AMREX_GPU_HOST_DEVICE
    Box& surroundingNodes (int dir) noexcept;

    //! Convert to CELL type in all directions.
    AMREX_GPU_HOST_DEVICE
    Box& enclosedCells () noexcept;

    //! Convert to CELL type in given direction.
    AMREX_GPU_HOST_DEVICE
    Box& enclosedCells (int dir) noexcept;

    /**
    * \brief Return Box that is intersection of this Box
    * and argument.  The Boxes MUST be of same type.
    */
    AMREX_GPU_HOST_DEVICE
    Box operator& (const Box& rhs) const noexcept { Box lhs(*this); lhs &= rhs; return lhs; }

    //! Intersect this Box with its argument. The Boxes MUST be of the same type.
    AMREX_GPU_HOST_DEVICE
    Box& operator&= (const Box& rhs) noexcept
    {
        BL_ASSERT(sameType(rhs));
        smallend.max(rhs.smallend);
        bigend.min(rhs.bigend);
        return *this;
    }

    //! for serialization
    static size_t linearSize() noexcept
    {
        size_t retval = 2*IntVect::linearSize();
        return retval;
    }

    //! for serialization
    void linearOut(void* a_buffer ) const noexcept;
    void linearIn(void* a_buffer ) noexcept;

    /**
    * \brief Modify Box to that of the minimum Box containing both
    * the original Box and the argument.
    * Both Boxes must have identical type.
    */
    AMREX_GPU_HOST_DEVICE
    Box& minBox (const Box& b) noexcept {
        // BoxArray may call this with not ok boxes.  BL_ASSERT(b.ok() && ok());
        BL_ASSERT(sameType(b));
        smallend.min(b.smallend);
        bigend.max(b.bigend);
        return *this;
    }

    //! Shift Box (relative) by given IntVect.
    AMREX_GPU_HOST_DEVICE
    Box& operator+= (const IntVect& v) noexcept { smallend += v; bigend   += v; return *this; }

    //! Shift Box (relative) by given IntVect.
    AMREX_GPU_HOST_DEVICE
    Box operator+ (const IntVect& v) const noexcept { Box r(*this); r += v; return r; }

    //! Shift Box (relative) by given IntVect.
    AMREX_GPU_HOST_DEVICE
    Box& operator-= (const IntVect& v) noexcept { smallend -= v; bigend -= v; return *this; }

    //! Shift Box (relative) by given IntVect.
    AMREX_GPU_HOST_DEVICE
    Box operator- (const IntVect& v) const noexcept { Box r(*this); r -= v; return r; }

    /**
    * \brief Chop the Box at the chop_pnt in the dir direction
    * returns one Box, modifies the object Box.
    * The union of the two is the original Box.
    * The modified Box is the low end, the returned Box
    * is the high end.  If type(dir) = CELL, the Boxes are disjoint
    * with the chop_pnt included in the high end (new Box).
    * It is an ERROR if chop_pnt is the low end of the orig Box.
    * If type(dir) = NODE, the chop_pnt is included in both Boxes
    * but is the only point in common.  It is also an error if the
    * chop_pnt is an end node of the Box.
    */
    AMREX_GPU_HOST_DEVICE
    Box chop (int dir, int chop_pnt) noexcept;

    /*
    * \brief Grow Box in all directions by given amount.
    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
    AMREX_GPU_HOST_DEVICE
    Box& grow (int i) noexcept { smallend.diagShift(-i); bigend.diagShift(i); return *this; }

    //! Grow Box in each direction by specified amount.
    AMREX_GPU_HOST_DEVICE
    Box& grow (const IntVect& v) noexcept { smallend -= v; bigend += v; return *this;}

    /**
    * \brief Grow the Box on the low and high end by n_cell cells
    * in direction idir.
    */
    AMREX_GPU_HOST_DEVICE
    Box& grow (int idir, int n_cell) noexcept { smallend.shift(idir, -n_cell); bigend.shift(idir, n_cell); return *this; }

    /**
    * \brief Grow the Box on the low end by n_cell cells in direction idir.
    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
    AMREX_GPU_HOST_DEVICE
    Box& growLo (int idir, int n_cell = 1) noexcept { smallend.shift(idir, -n_cell); return *this; }

    /**
    * \brief Grow the Box on the high end by n_cell cells in
    * direction idir.  NOTE: n_cell negative shrinks the Box by that
    * number of cells.
    */
    AMREX_GPU_HOST_DEVICE
    Box& growHi (int idir, int n_cell = 1) noexcept { bigend.shift(idir,n_cell); return *this; }

    //! Grow in the direction of the given face.
    AMREX_GPU_HOST_DEVICE
    Box& grow (Orientation face, int n_cell = 1) noexcept {
        int idir = face.coordDir();
        if (face.isLow())
            smallend.shift(idir, -n_cell);
        else
            bigend.shift(idir,n_cell);
        return *this;
    }

    /**
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    AMREX_GPU_HOST_DEVICE
    Box& refine (int ref_ratio) noexcept {
        return this->refine(IntVect(ref_ratio));
    }

    /*
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    AMREX_GPU_HOST_DEVICE
    Box& refine (const IntVect& ref_ratio) noexcept;

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    AMREX_GPU_HOST_DEVICE
    Box& coarsen (int ref_ratio) noexcept {
        return this->coarsen(IntVect(ref_ratio));
    }

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    AMREX_GPU_HOST_DEVICE
    Box& coarsen (const IntVect& ref_ratio) noexcept;

    /**
    * \brief Step through the rectangle.  It is a runtime error to give
    * a point not inside rectangle.  Iteration may not be efficient.
    */
    AMREX_GPU_HOST_DEVICE
    void next (IntVect &) const noexcept;

    /**
    * \brief This static member function returns a constant reference to
    * an object of type Box representing the unit box in
    * AMREX_SPACEDIM-dimensional space.
    */
    AMREX_GPU_HOST_DEVICE
    static Box TheUnitBox () noexcept {
        return Box(IntVect::TheZeroVector(),IntVect::TheZeroVector());
    }

    AMREX_GPU_HOST_DEVICE
    bool isSquare() const noexcept;

    AMREX_GPU_HOST_DEVICE
    bool coarsenable(int refrat, int min_width=1) const noexcept {
        return coarsenable(IntVect(refrat),min_width);
    }

    AMREX_GPU_HOST_DEVICE
    bool coarsenable(const IntVect& refrat, int min_width=1) const noexcept
    {
        if (!size().allGE(refrat*min_width)) {
            return false;
        } else {
            Box testBox = *this;
            testBox.coarsen(refrat);
            testBox.refine (refrat);
            return (*this == testBox);
        }
    }

    AMREX_GPU_HOST_DEVICE
    void normalize () noexcept
    {
        for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            if (this->length(idim) == 0) {
                this->growHi(idim,1);
            }
        }
    }

    AMREX_GPU_HOST_DEVICE friend AMREX_FORCE_INLINE Dim3 lbound (Box const& box) noexcept;
    AMREX_GPU_HOST_DEVICE friend AMREX_FORCE_INLINE Dim3 ubound (Box const& box) noexcept;
    AMREX_GPU_HOST_DEVICE friend AMREX_FORCE_INLINE Dim3 begin  (Box const& box) noexcept;
    AMREX_GPU_HOST_DEVICE friend AMREX_FORCE_INLINE Dim3 end    (Box const& box) noexcept;
    AMREX_GPU_HOST_DEVICE friend AMREX_FORCE_INLINE Dim3 length (Box const& box) noexcept;

private:
    IntVect         smallend;
    IntVect         bigend;
    IndexType       btype;
};

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::refine (const IntVect& ref_ratio) noexcept
{
    if (ref_ratio != 1) {
        IntVect shft(1);
        shft -= btype.ixType();
        smallend *= ref_ratio;
        bigend += shft;
        bigend *= ref_ratio;
        bigend -= shft;
    }
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::coarsen (const IntVect& ref_ratio) noexcept
{
    if (ref_ratio != 1)
    {
        smallend.coarsen(ref_ratio);

        if (btype.any())
        {
            IntVect off(0);
            for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
            {
                if (btype[dir])
                    if (bigend[dir]%ref_ratio[dir])
                        off.setVal(dir,1);
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

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::convert (const IntVect& typ) noexcept
{
    BL_ASSERT(typ.allGE(IntVect::TheZeroVector()) && typ.allLE(IntVect::TheUnitVector()));
    IntVect shft(typ - btype.ixType());
    bigend += shft;
    btype = IndexType(typ);
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::convert (IndexType t) noexcept
{
   for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
   {
      const unsigned int typ = t[dir];
      const unsigned int bitval = btype[dir];
      const int off = typ - bitval;
      bigend.shift(dir,off);
      btype.setType(dir, (IndexType::CellIndex) typ);
   }
   return *this;
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::surroundingNodes (int dir) noexcept
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

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::surroundingNodes () noexcept
{
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
        if ((btype[i] == 0))
            bigend.shift(i,1);
    btype.setall();
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::enclosedCells (int dir) noexcept
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

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::enclosedCells () noexcept
{
    for (int i = 0 ; i < AMREX_SPACEDIM; ++i)
        if (btype[i])
            bigend.shift(i,-1);
    btype.clear();
    return *this;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
long
Box::index (const IntVect& v) const noexcept
{
    long result = v[0]-smallend[0];
#if   AMREX_SPACEDIM==2
    result += length(0)*(v[1]-smallend[1]);
#elif AMREX_SPACEDIM==3
    result += length(0)*(v[1]-smallend[1]
                         +(v[2]-smallend[2])*length(1));
#endif
    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVect
Box::atOffset (long offset) const noexcept
{
    auto len = this->length3d();
    auto lo  = this->loVect3d();
    long k = offset / (len[0]*len[1]);
    long j = (offset - k*(len[0]*len[1])) / len[0];
    long i = (offset - k*(len[0]*len[1])) - j*len[0];
    return {AMREX_D_DECL(static_cast<int>(i+lo[0]),
                         static_cast<int>(j+lo[1]),
                         static_cast<int>(k+lo[2])) };
}

AMREX_GPU_HOST_DEVICE
inline
GpuArray<int,3>
Box::atOffset3d (long offset) const noexcept
{
    auto len = this->length3d();
    auto lo  = this->loVect3d();
    long k = offset / (len[0]*len[1]);
    long j = (offset - k*(len[0]*len[1])) / len[0];
    long i = (offset - k*(len[0]*len[1])) - j*len[0];
    return {static_cast<int>(i+lo[0]),
            static_cast<int>(j+lo[1]),
            static_cast<int>(k+lo[2]) };
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::setRange (int dir,
               int sm_index,
               int n_cells) noexcept
{
    smallend.setVal(dir,sm_index);
    bigend.setVal(dir,sm_index+n_cells-1);
    return *this;
}

inline
void
Box::linearOut(void* a_buffer ) const noexcept
{
    unsigned char* buf = (unsigned char*) a_buffer;
    const IntVect& ivlo = this->smallEnd();
    const IntVect& ivhi = this->bigEnd();
    ivlo.linearOut(buf);
    buf += ivlo.linearSize();
    ivhi.linearOut(buf);
    buf += ivhi.linearSize();
}

inline
void
Box::linearIn(void* a_buffer ) noexcept
{
    unsigned char* buf = (unsigned char*) a_buffer;
    IntVect ivlo;
    IntVect ivhi;
    ivlo.linearIn(buf);
    buf += ivlo.linearSize();
    ivhi.linearIn(buf);
    buf += ivhi.linearSize();

    *this = Box(ivlo, ivhi);
}

AMREX_GPU_HOST_DEVICE
inline
void
Box::next (IntVect& p) const noexcept
{
    BL_ASSERT(contains(p));

    ++p[0];

#if (AMREX_SPACEDIM >= 2)
    if (p[0] > bigend[0])
    {
	p[0] = smallend[0];
	++p[1];
#if (AMREX_SPACEDIM == 3)
	if (p[1] > bigend[1])
	{
	    p[1] = smallend[1];
	    ++p[2];
	}
#endif
    }
#endif
}

AMREX_GPU_HOST_DEVICE
inline
bool
Box::isSquare () const noexcept
{
    const IntVect& sz = this->size();
#if AMREX_SPACEDIM==1
    return false; // can't build a square in 1-D
#elif AMREX_SPACEDIM==2
    return (sz[0] == sz[1]);
#elif AMREX_SPACEDIM==3
    return (sz[0] == sz[1] && (sz[1] == sz[2]));
#endif
}

//
// Modified Box is low end, returned Box is high end.
// If CELL: chop_pnt included in high end.
// If NODE: chop_pnt included in both Boxes.
//

AMREX_GPU_HOST_DEVICE
inline
Box
Box::chop (int dir, int chop_pnt) noexcept
{
    //
    // Define new high end Box including chop_pnt.
    //
    IntVect sm(smallend);
    IntVect bg(bigend);
    sm.setVal(dir,chop_pnt);
    if (btype[dir])
    {
        //
        // NODE centered Box.
        //
        BL_ASSERT(chop_pnt > smallend[dir] && chop_pnt < bigend[dir]);
        //
        // Shrink original Box to just contain chop_pnt.
        //
        bigend.setVal(dir,chop_pnt);
    }
    else
    {
        //
        // CELL centered Box.
        //
        BL_ASSERT(chop_pnt > smallend[dir] && chop_pnt <= bigend[dir]);
        //
        // Shrink origional Box to one below chop_pnt.
        //
        bigend.setVal(dir,chop_pnt-1);
    }
    return Box(sm,bg,btype);
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::shiftHalf (int dir, int nzones) noexcept
{
    const int nbit = (nzones<0 ? -nzones : nzones)%2;
    int nshift = nzones/2;
    //
    // Toggle btyp bit if nzones is odd.
    //
    const unsigned int bit_dir = btype[dir];
    if (nbit)
        btype.flip(dir);
    if (nzones < 0)
        nshift -= (bit_dir ? nbit : 0);
    else
        nshift += (bit_dir ? 0 : nbit);
    smallend.shift(dir,nshift);
    bigend.shift(dir,nshift);
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
Box&
Box::shiftHalf (const IntVect& nz) noexcept
{
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        shiftHalf(i,nz[i]);
   return *this;
}

class BoxCommHelper
{
public:

    explicit BoxCommHelper (const Box& bx, int* p_ = 0);

    int* data () const& noexcept { return p; }
    int* data () && = delete;

    Box make_box () const noexcept {
	return Box(IntVect(p), IntVect(p+AMREX_SPACEDIM), IntVect(p+2*AMREX_SPACEDIM));
    }

    static int size () noexcept { return 3*AMREX_SPACEDIM; }

private:
    int* p;
    std::vector<int> v;
};

class BoxConverter {
public:
    virtual Box doit (const Box& fine) const = 0;
    virtual BoxConverter* clone () const = 0;
    virtual ~BoxConverter () = 0;
};

inline BoxConverter::~BoxConverter () { }

void AllGatherBoxes (Vector<Box>& bxs);

    /**
    * \brief Grow Box in all directions by given amount.

    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
AMREX_GPU_HOST_DEVICE
inline
Box grow (const Box& b, int i) noexcept
{
    Box result = b;
    result.grow(i);
    return result;
}

    //! Grow Box in each direction by specified amount.
AMREX_GPU_HOST_DEVICE
inline
Box grow (const Box& b, const IntVect& v) noexcept
{
    Box result = b;
    result.grow(v);
    return result;
}

    //! Grow Box in direction idir be n_cell cells
AMREX_GPU_HOST_DEVICE
inline
Box grow (const Box& b, int idir, int n_cell) noexcept
{
    Box result = b;
    result.grow(idir, n_cell);
    return result;
}

AMREX_GPU_HOST_DEVICE
inline
Box growLo (const Box& b, int idir, int n_cell) noexcept
{
    Box result = b;
    result.growLo(idir, n_cell);
    return result;
}

AMREX_GPU_HOST_DEVICE
inline
Box growHi (const Box& b, int idir, int n_cell) noexcept
{
    Box result = b;
    result.growHi(idir, n_cell);
    return result;
}

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
AMREX_GPU_HOST_DEVICE
inline
Box coarsen (const Box& b, int ref_ratio) noexcept
{
    Box result = b;
    result.coarsen(IntVect(ref_ratio));
    return result;
}

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
AMREX_GPU_HOST_DEVICE
inline
Box coarsen (const Box& b, const IntVect& ref_ratio) noexcept
{
    Box result = b;
    result.coarsen(ref_ratio);
    return result;
}

    /**
    * Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
AMREX_GPU_HOST_DEVICE
inline
Box refine (const Box& b, int ref_ratio) noexcept
{
    Box result = b;
    result.refine(IntVect(ref_ratio));
    return result;
}

    /**
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
AMREX_GPU_HOST_DEVICE
inline
Box refine (const Box& b, const IntVect& ref_ratio) noexcept
{
    Box result = b;
    result.refine(ref_ratio);
    return result;
}

    //! Return a Box with indices shifted by nzones in dir direction.
AMREX_GPU_HOST_DEVICE
inline
Box shift (const Box& b, int dir, int nzones) noexcept
{
    Box result = b;
    result.shift(dir, nzones);
    return result;
}

AMREX_GPU_HOST_DEVICE
inline
Box shift (const Box& b, const IntVect& nzones) noexcept
{
    Box result = b;
    result.shift(nzones);
    return result;
}

    /**
    * \brief Returns a Box with NODE based coordinates in direction dir
    * that encloses Box b.  NOTE: equivalent to b.convert(dir,NODE)
    * NOTE: error if b.type(dir) == NODE.
    */
AMREX_GPU_HOST_DEVICE
inline
Box surroundingNodes (const Box& b, int dir) noexcept
{
    Box bx(b);
    bx.surroundingNodes(dir);
    return bx;
}

    /**
    * \brief Returns a Box with NODE based coordinates in all
    * directions that encloses Box b.
    */
AMREX_GPU_HOST_DEVICE
inline
Box surroundingNodes (const Box& b) noexcept
{
    Box bx(b);
    bx.surroundingNodes();
    return bx;
}

    //! Returns a Box with different type
AMREX_GPU_HOST_DEVICE
inline
Box convert (const Box& b, const IntVect& typ) noexcept
{
    Box bx(b);
    bx.convert(typ);
    return bx;
}

AMREX_GPU_HOST_DEVICE
inline
Box convert (const Box& b, const IndexType& typ) noexcept
{
    Box bx(b);
    bx.convert(typ);
    return bx;
}

    /**
    * \brief Returns a Box with CELL based coordinates in
    * direction dir that is enclosed by b.
    * NOTE: equivalent to b.convert(dir,CELL)
    * NOTE: error if b.type(dir) == CELL.
    */
AMREX_GPU_HOST_DEVICE
inline
Box enclosedCells (const Box& b, int dir) noexcept
{
    Box bx(b);
    bx.enclosedCells(dir);
    return bx;
}

    /**
    * \brief Returns a Box with CELL based coordinates in all
    * directions that is enclosed by b.
    */
AMREX_GPU_HOST_DEVICE
inline
Box enclosedCells (const Box& b) noexcept
{
    Box bx(b);
    bx.enclosedCells();
    return bx;
}

    /**
    * \brief Returns the edge-centered Box (in direction dir) defining
    * the low side of Box b.
    */
AMREX_GPU_HOST_DEVICE
inline
Box bdryLo (const Box& b, int dir, int len=1) noexcept
{
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    int sm = low[dir];
    low.setVal(dir,sm-len+1);
    hi.setVal(dir,sm);
    //
    // set dir'th bit to 1 = IndexType::NODE.
    //
    IndexType typ(b.ixType());
    typ.set(dir);
    return Box(low,hi,typ);
}

    /**
    * \brief Returns the edge-centered Box (in direction dir) defining
    * the high side of Box b.
    */
AMREX_GPU_HOST_DEVICE
inline
Box bdryHi (const Box& b, int dir, int len=1) noexcept
{
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    unsigned int bitval = b.type()[dir];
    int bg              = hi[dir]  + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg+len-1);
    //
    // Set dir'th bit to 1 = IndexType::NODE.
    //
    IndexType typ(b.ixType());
    typ.set(dir);
    return Box(low,hi,typ);
}

    /**
    * \brief Similar to bdryLo and bdryHi except that it operates on the
    * given face of  box b.
    */
AMREX_GPU_HOST_DEVICE
inline
Box bdryNode (const Box& b, Orientation face, int len=1) noexcept
{
    int dir = face.coordDir();
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    if (face.isLow())
    {
        int sm = low[dir];
	low.setVal(dir,sm-len+1);
        hi.setVal(dir,sm);
    }
    else
    {
        unsigned int bitval = b.type()[dir];
        int bg              = hi[dir]  + 1 - bitval%2;
        low.setVal(dir,bg);
        hi.setVal(dir,bg+len-1);
    }
    //
    // Set dir'th bit to 1 = IndexType::NODE.
    //
    IndexType typ(b.ixType());
    typ.set(dir);
    return Box(low,hi,typ);
}

    /**
    * \brief Returns the cell centered Box of length len adjacent
    * to b on the low end along the coordinate direction dir.
    * The return Box is identical to b in the other directions.
    * The return Box and b have an empty intersection.
    * NOTE:  len >= 1
    * NOTE:  Box retval = b.adjCellLo(b,dir,len)
    * is equivalent to the following set of operations:
    * Box retval(b);
    * retval.convert(dir,Box::CELL);
    * retval.setrange(dir,retval.smallEnd(dir)-len,len);
    */
AMREX_GPU_HOST_DEVICE
inline
Box adjCellLo (const Box& b, int dir, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    int sm = low[dir];
    low.setVal(dir,sm - len);
    hi.setVal(dir,sm - 1);
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexType typ(b.ixType());
    typ.unset(dir);
    return Box(low,hi,typ);
}

    //! Similar to adjCellLo but builds an adjacent Box on the high end.
AMREX_GPU_HOST_DEVICE
inline
Box adjCellHi (const Box& b, int dir, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    unsigned int bitval = b.type()[dir];
    int bg              = hi[dir]  + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg + len - 1);
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexType typ(b.ixType());
    typ.unset(dir);
    return Box(low,hi,typ);
}

    //! Similar to adjCellLo and adjCellHi; operates on given face.
AMREX_GPU_HOST_DEVICE
inline
Box adjCell (const Box& b, Orientation face, int len=1) noexcept
{
    BL_ASSERT(len > 0);
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    int dir = face.coordDir();
    if (face.isLow())
    {
        int sm = low[dir];
        low.setVal(dir,sm - len);
        hi.setVal(dir,sm - 1);
    }
    else
    {
        unsigned int bitval = b.type()[dir];
        int bg              = hi[dir]  + 1 - bitval%2;
        low.setVal(dir,bg);
        hi.setVal(dir,bg + len - 1);
    }
    //
    // Set dir'th bit to 0 = IndexType::CELL.
    //
    IndexType typ(b.ixType());
    typ.unset(dir);
    return Box(low,hi,typ);
}

    /**
    * \brief Modify Box to that of the minimum Box containing both
    * the original Box and the argument.
    * Both Boxes must have identical type.
    */
AMREX_GPU_HOST_DEVICE
inline
Box minBox (const Box& b1, const Box& b2) noexcept
{
    Box result = b1;
    result.minBox(b2);
    return result;
}

    //! Write an ASCII representation to the ostream.
std::ostream& operator<< (std::ostream& os, const Box& bx);

    //! Read from istream.
std::istream& operator>> (std::istream& os, Box& bx);

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 lbound (Box const& box) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {box.smallend[0], 0, 0};
#elif (AMREX_SPACEDIM == 2)
    return {box.smallend[0], box.smallend[1], 0};
#elif (AMREX_SPACEDIM == 3)
    return {box.smallend[0], box.smallend[1], box.smallend[2]};
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 ubound (Box const& box) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {box.bigend[0], 0, 0};
#elif (AMREX_SPACEDIM == 2)
    return {box.bigend[0], box.bigend[1], 0};
#elif (AMREX_SPACEDIM == 3)
    return {box.bigend[0], box.bigend[1], box.bigend[2]};
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 begin (Box const& box) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {box.smallend[0], 0, 0};
#elif (AMREX_SPACEDIM == 2)
    return {box.smallend[0], box.smallend[1], 0};
#elif (AMREX_SPACEDIM == 3)
    return {box.smallend[0], box.smallend[1], box.smallend[2]};
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 end (Box const& box) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {box.bigend[0]+1, 1, 1};
#elif (AMREX_SPACEDIM == 2)
    return {box.bigend[0]+1, box.bigend[1]+1, 1};
#elif (AMREX_SPACEDIM == 3)
    return {box.bigend[0]+1, box.bigend[1]+1, box.bigend[2]+1};
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Dim3 length (Box const& box) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {box.bigend[0]-box.smallend[0]+1, 1, 1};
#elif (AMREX_SPACEDIM == 2)
    return {box.bigend[0]-box.smallend[0]+1,
            box.bigend[1]-box.smallend[1]+1, 1};
#elif (AMREX_SPACEDIM == 3)
    return {box.bigend[0]-box.smallend[0]+1,
            box.bigend[1]-box.smallend[1]+1,
            box.bigend[2]-box.smallend[2]+1};
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
Dim3 length (Box index_bounds, Box const& work_domain) noexcept
{
    index_bounds.setType(work_domain.ixType());
    index_bounds &= work_domain;
    return length(index_bounds);
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
Box getThreadBox (Box index_bounds, Box const& work_domain) noexcept
{
    index_bounds.setType(work_domain.ixType());
    index_bounds &= work_domain;
    return index_bounds;
}

    /**
    * Returns a Box that covers all the argument Boxes in index
    * space. The types are ignored. Thus, the arguments can have
    * different index types, and the returned Box's index type has no
    * meaning.
    */
inline
Box getIndexBounds (Box const& b1, Box const& b2) noexcept
{
    Box b = b1;
    b.setType(b2.ixType());
    b.minBox(b2);
    return b;
}

inline
Box getIndexBounds (Box const& b1, Box const& b2, Box const& b3) noexcept
{
    return getIndexBounds(getIndexBounds(b1,b2),b3);
}

inline
Box getIndexBounds (Box const& b1) noexcept
{
    return b1;
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect getCell (Box const* boxes, int nboxes, long icell) noexcept
{
    int ibox;
    long ncells_subtotal = 0;
    long offset = 0;
    for (ibox = 0; ibox < nboxes; ++ibox) {
        const long n = boxes[ibox].numPts();
        ncells_subtotal += n;
        if (icell < ncells_subtotal) {
            offset = icell - (ncells_subtotal - n);
            break;
        }
    }
    return boxes[ibox].atOffset(offset);
}

}

#endif /*BL_BOX_H*/
