
#ifndef BL_BOX_H
#define BL_BOX_H

#include <iosfwd>

#include <AMReX_ccse-mpi.H>
#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Orientation.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

namespace amrex
{
    class Box;
    class BoxCommHelper;

    /**
    * \brief Grow Box in all directions by given amount.
  
    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
    Box grow (const Box& b, int i);

    //! Grow Box in each direction by specified amount.
    Box grow (const Box& b, const IntVect& v);

    //! Grow Box in direction idir be n_cell cells
    Box grow (const Box& b, int idir, int n_cell);

    Box growLo (const Box& b, int idir, int n_cell);
    Box growHi (const Box& b, int idir, int n_cell);

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    Box coarsen (const Box& b, int refinement_ratio);

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    Box coarsen (const Box& b, const IntVect& refinement_ratio);

    /**
    * Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    Box refine (const Box& b, int refinement_ratio);

    //! Return a Box with indices shifted by nzones in dir direction.
    Box shift (const Box& b, int dir, int nzones);

    /**
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    Box refine (const Box& b, const IntVect& refinement_ratio);

    /**
    * \brief Returns a Box with NODE based coordinates in direction dir
    * that encloses Box b.  NOTE: equivalent to b.convert(dir,NODE)
    * NOTE: error if b.type(dir) == NODE.
    */
    Box surroundingNodes (const Box& b, int dir);

    /**
    * \brief Returns a Box with NODE based coordinates in all
    * directions that encloses Box b.
    */
    Box surroundingNodes (const Box& b);

    //! Returns a Box with different type
    Box convert (const Box& b, const IntVect& typ);
    Box convert (const Box& b, const IndexType& typ);

    /**
    * \brief Returns a Box with CELL based coordinates in 
    * direction dir that is enclosed by b.
    * NOTE: equivalent to b.convert(dir,CELL)
    * NOTE: error if b.type(dir) == CELL.
    */
    Box enclosedCells (const Box& b, int dir);

    /**
    * \brief Returns a Box with CELL based coordinates in all
    * directions that is enclosed by b.
    */
    Box enclosedCells (const Box& b);

    /**
    * \brief Returns the edge-centered Box (in direction dir) defining
    * the low side of Box b.
    */
    Box bdryLo (const Box& b, int dir, int len=1);

    /**
    * \brief Returns the edge-centered Box (in direction dir) defining
    * the high side of Box b.
    */
    Box bdryHi (const Box& b, int dir, int len=1);

    /**
    * \brief Similar to bdryLo and bdryHi except that it operates on the
    * given face of  box b.
    */
    Box bdryNode (const Box& b, Orientation face, int len=1);

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
    Box adjCellLo (const Box& b, int dir, int len=1);

    //! Similar to adjCellLo but builds an adjacent Box on the high end.
    Box adjCellHi (const Box& b, int dir, int len=1);

    //! Similar to adjCellLo and adjCellHi; operates on given face.
    Box adjCell (const Box& b, Orientation face, int len=1);

    /**
    * \brief Modify Box to that of the minimum Box containing both
    * the original Box and the argument.
    * Both Boxes must have identical type.
    */
    Box minBox (const Box& b1, const Box& b2);

    Array<int> SerializeBox(const Box &b);
    int SerializeBoxSize();
    Box UnSerializeBox(const Array<int> &serarray);


    //! Write an ASCII representation to the ostream.
    std::ostream& operator<< (std::ostream& os, const Box& bx);

    //! Read from istream.
    std::istream& operator>> (std::istream& os, Box& bx);


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
class Box
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<Box>::type();
    friend class BoxCommHelper;

public:
    /*
    * \brief The default constructor.  For safety, the constructed Box is  
    * invalid and may be tested for validity with ok(). 
    * DO NOT CHANGE THIS BEHAVIOR!
    */
    Box ();

    //! Construct cell-centered type Box.
    Box (const IntVect& small,
         const IntVect& big);

    //! Construct box with specified lengths.
    Box (const IntVect& small,
         const int*     vec_len);

    /**
    * \brief Construct Box with given type.  small and big are expected
    * to be consistent with given type.
    */
    Box (const IntVect& small,
         const IntVect& big,
         const IntVect& typ);
    //
    // Construct dimension specific Boxes.
    //
    Box (const IntVect& small,
         const IntVect& big,
         IndexType      t);

    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    //! Get the smallend of the box.
    const IntVect& smallEnd () const& { return smallend; }
    const IntVect& smallEnd () && = delete;

    //! Returns the coordinate of the low end in the given direction.
    int smallEnd (int dir) const& { return smallend[dir]; }

    //! Get the bigend.
    const IntVect& bigEnd () const& { return bigend; }
    const IntVect& bigEnd () && = delete;

    //! Returns the coordinate of the high end in the given direction.
    int bigEnd (int dir) const { return bigend[dir]; }

    //! Returns the indexing type.
    IndexType ixType () const { return btype; }

    //! Returns the indexing type.
    IntVect type () const { return btype.ixType(); }

    //! Returns the indexing type in the specified direction.
    IndexType::CellIndex type (int dir) const { return btype.ixType(dir); }

    //! Return the length of the Box.
    IntVect size () const
    {
        return IntVect(AMREX_D_DECL(bigend[0]-smallend[0] + 1,
                              bigend[1]-smallend[1] + 1,
                              bigend[2]-smallend[2] + 1));
    }

    //! Return the length of the Box in given direction.
    int length (int dir) const { return bigend[dir] - smallend[dir] + 1; }

    //! Returns a constant pointer the array of low end coordinates. Useful for calls to FORTRAN.
    const int* loVect () const& { return smallend.getVect(); }
    const int* loVect () && = delete;
    //! Returns a constant pointer the array of high end coordinates. Useful for calls to FORTRAN.
    const int* hiVect () const& { return bigend.getVect(); }
    const int* hiVect () && = delete;

    //! Returns the coordinate normal to given face.
    int operator[] (Orientation face) const;

    //! Checks if it is an empty box.
    bool isEmpty () const { return !ok(); }

    //! Checks if it is a proper Box (including a valid type).
    bool ok () const { return bigend.allGE(smallend) && btype.ok(); }

    //! Returns true if argument is contained within Box.
    bool contains (const IntVect& p) const { return p.allGE(smallend) && p.allLE(bigend); }

    /** \brief Returns true if argument is contained within Box.
    * It is an error if the Boxes have different types.
    */
    bool contains (const Box& b) const
    { 
	BL_ASSERT(sameType(b)); 
	return b.smallend.allGE(smallend) && b.bigend.allLE(bigend);
    } 

    //! Returns true if argument is strictly contained within Box.
    bool strictly_contains (const IntVect& p) const { return p.allGT(smallend) && p.allLT(bigend); }

    /**
    * \brief Returns true if argument is strictly contained within Box.
    * It is an error if the Boxes have different types.
    */
    bool strictly_contains (const Box& b) const 
    { 
	BL_ASSERT(sameType(b));
	return b.smallend.allGT(smallend) && b.bigend.allLT(bigend);
    }
    /**
    * \brief Returns true if Boxes have non-null intersections.
    * It is an error if the Boxes have different types.
    */
    bool intersects (const Box& b) const { Box isect(*this); isect &= b; return isect.ok(); }

    /**
    * \brief Returns true is Boxes same size, ie translates of each other,.
    * It is an error if they have different types.
    */
    bool sameSize (const Box& b) const;

    //! Returns true if Boxes have same type.
    bool sameType (const Box &b) const { return btype == b.btype; }

    //! Returns true if Boxes are identical (including type).
    bool operator== (const Box& b) const { return smallend == b.smallend && bigend == b.bigend && b.btype == btype; }

    //! Returns true if Boxes differ (including type).
    bool operator!= (const Box& b) const { return !operator==(b); }

    //! Returns true if Box is cell-centered in all indexing directions.
    bool cellCentered () const { return !btype.any(); }

    /**
    * \brief Returns the number of points contained in the Box, else
    * abort()s if the number cannot be represented in a long.
    */
    long numPts () const;

    /**
    * \brief Returns the number of points contained in the Box.
    * This is inteded for use only in diagnostic messages.
    */
    double d_numPts () const;

    /**
    * \brief Return the volume, in indexing space, of region enclosed by
    * this Box.   This is identical to numPts() for CELL centered
    * Box; otherwise, numPts() > volume().  Abort()s if the number
    * cannot be represented in a long.
    */
    long volume () const;

    /**
    * \brief Returns length of longest side.  dir is modified to give
    * direction with longest side: 0...SPACEDIM-1. Ignores type.
    */
    int longside (int& dir) const;

    //! Returns length of longest side.  Ignores type.
    int longside () const;

    /**
    * \brief Returns length of shortest side.  dir is modified to give
    *  direction with shortest side: 0...SPACEDIM-1.  Ignores type.
    */
    int shortside (int& dir) const;

    //! Returns length of shortest side.  Ignores type.
    int shortside () const;

    /**
    * \brief Returns offset of point from smallend; i.e.
    * index(smallend) -> 0, bigend would return numPts()-1.
    * Is used in accessing FArrayBox.
    */
    long index (const IntVect& v) const;

    //! Redefine the small end of the Box.
    Box& setSmall (const IntVect& sm) { smallend = sm; return *this; }

    //! Redefine the small end of the Box.
    Box& setSmall (int dir, int sm_index)  { smallend.setVal(dir,sm_index); return *this; }

    //! Redefine the big end of the Box.
    Box& setBig (const IntVect& bg) { bigend = bg; return *this; }

    //! Redefine the big end of the Box.
    Box& setBig (int dir, int bg_index) { bigend.setVal(dir,bg_index); return *this; }

    /**
    * \brief Set the entire range in a given direction, starting at
    * sm_index with length n_cells.  NOTE: This will yield an
    * illegal Box if n_cells <= 0.
    */
    Box& setRange (int dir,
                   int sm_index,
                   int n_cells = 1);

    //! Set indexing type
    Box& setType (const IndexType& t) { btype = t; return *this; }

    //! Shift this Box nzones indexing positions in coordinate direction dir.
    Box& shift (int dir, int nzones) { smallend.shift(dir,nzones); bigend.shift(dir,nzones); return *this; }

    //! Equivalent to b.shift(0,iv[0]).shift(1,iv[1]) ....
    Box& shift (const IntVect& iv) { smallend.shift(iv); bigend.shift(iv); return *this; }

    /**
    * \brief This member shifts the Box by "half" indicies, thereby
    * converting the Box from type CELL to NODE and visa-versa.
    * b.shiftHalf(0,1)  shifts b to the right by 1/2 cells.
    * b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
    * NOTE: If num_halfs is EVEN the shift is num_halfs/2 full
    * zones and hence will not change the type.
    * This is: b.shifthalf(4) == b.shift(2).
    */
    Box& shiftHalf (int dir, int num_halfs);

    //! Equivalent to b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...
    Box& shiftHalf (const IntVect& iv);

    /**
    * \brief Convert the Box from the current type into the
    * argument type.  This may change the Box coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    Box& convert (IndexType typ);

    /**
    * \brief Convert the Box from the current type into the
    * argument type.  This may change the Box coordinates:
    * type CELL -> NODE : increase coordinate by one on high end
    * type NODE -> CELL : reduce coordinate by one on high end
    * other type mappings make no change.
    */
    Box& convert (const IntVect& typ);

    //! Convert to NODE type in all directions.
    Box& surroundingNodes ();

    //! Convert to NODE type in given direction.
    Box& surroundingNodes (int dir);

    //! Convert to CELL type in all directions.
    Box& enclosedCells ();

    //! Convert to CELL type in given direction.
    Box& enclosedCells (int dir);

    /**
    * \brief Return Box that is intersection of this Box
    * and argument.  The Boxes MUST be of same type.
    */
    Box operator& (const Box& rhs) const { Box lhs(*this); lhs &= rhs; return lhs; }

    //! Intersect this Box with its argument. The Boxes MUST be of the same type.
    Box& operator&= (const Box& rhs)
    {
        BL_ASSERT(sameType(rhs));
        smallend.max(rhs.smallend);
        bigend.min(rhs.bigend);
        return *this;
    }

    /// for serialization
    static size_t 
    linearSize()
    {
      size_t retval = 2*IntVect::linearSize();
      return retval;
    }

    /// for serialization
    void 
    linearOut(void* a_buffer ) const
    {
      unsigned char* buf = (unsigned char*) a_buffer;
      const IntVect& ivlo = this->smallEnd();
      const IntVect& ivhi = this->bigEnd();
      ivlo.linearOut(buf);
      buf += ivlo.linearSize();
      ivhi.linearOut(buf);
      buf += ivhi.linearSize();
    }

    /// for serialization
    void 
    linearIn(void* a_buffer )
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

    /**
    * \brief Modify Box to that of the minimum Box containing both
    * the original Box and the argument.
    * Both Boxes must have identical type.
    */
    Box& minBox (const Box&);

    //! Shift Box (relative) by given IntVect.
    Box& operator+= (const IntVect& v) { smallend += v; bigend   += v; return *this; }

    //! Shift Box (relative) by given IntVect.
    Box operator+ (const IntVect& v) const { Box r(*this); r += v; return r; }

    //! Shift Box (relative) by given IntVect.
    Box& operator-= (const IntVect& v) { smallend -= v; bigend -= v; return *this; }

    //! Shift Box (relative) by given IntVect.
    Box operator- (const IntVect& v) const { Box r(*this); r -= v; return r; }

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
    Box chop (int dir, int chop_pnt);

    /*
    * \brief Grow Box in all directions by given amount.
    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
    Box& grow (int i) { smallend.diagShift(-i); bigend.diagShift(i); return *this; }

    //! Grow Box in each direction by specified amount.
    Box& grow (const IntVect& v) { smallend -= v; bigend += v; return *this;}

    /**
    * \brief Grow the Box on the low and high end by n_cell cells
    * in direction idir.
    */
    Box& grow (int idir, int n_cell) { smallend.shift(idir, -n_cell); bigend.shift(idir, n_cell); return *this; }

    /**
    * \brief Grow the Box on the low end by n_cell cells in direction idir.
    * NOTE: n_cell negative shrinks the Box by that number of cells.
    */
    Box& growLo (int idir, int n_cell = 1) { smallend.shift(idir, -n_cell); return *this; }

    /**
    * \brief Grow the Box on the high end by n_cell cells in
    * direction idir.  NOTE: n_cell negative shrinks the Box by that
    * number of cells.
    */
    Box& growHi (int idir, int n_cell = 1) { bigend.shift(idir,n_cell); return *this; }

    //! Grow in the direction of the given face.
    Box& grow (Orientation face, int n_cell = 1);

    /**
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    Box& refine (int refinement_ratio);

    /*
    * \brief Refine Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
    * hi <- (hi+1)*ratio - 1.
    * NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
    * hi <- hi*ratio.
    */
    Box& refine (const IntVect& refinement_ratio);

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    Box& coarsen (int refinement_ratio);

    /**
    * \brief Coarsen Box by given (positive) refinement ratio.
    * NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
    * hi <- hi/ratio.
    * NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
    * hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
    * That is, refinement of coarsened Box must contain
    * the original Box.
    */
    Box& coarsen (const IntVect& refinement_ratio);

    /**
    * \brief Step through the rectangle.  It is a runtime error to give
    * a point not inside rectangle.  Iteration may not be efficient.
    */
    void next (IntVect &) const;

    /**
    * \brief This static member function returns a constant reference to 
    * an object of type Box representing the unit box in
    * BL_SPACEDIM-dimensional space.
    */
    static const Box& TheUnitBox ();

    bool isSquare() const;

    bool coarsenable(int refrat) const {
        return coarsenable(IntVect(AMREX_D_DECL(refrat,refrat,refrat)));
    }

    bool coarsenable(const IntVect& refrat) const
    {
      Box testBox = *this;
      testBox.coarsen(refrat);
      testBox.refine (refrat);
      return (*this == testBox);
    }

    void normalize ()
    {
        for (int idim=0; idim < BL_SPACEDIM; ++idim) {
            if (this->length(idim) == 0) {
                this->growHi(idim,1);
            }
        }
    }

private:
    IntVect         smallend;
    IntVect         bigend;
    IndexType       btype;
};

inline
long
Box::index (const IntVect& v) const
{
    long result = v[0]-smallend[0];
#if   BL_SPACEDIM==2
    result += length(0)*(v[1]-smallend[1]);
#elif BL_SPACEDIM==3
    result += length(0)*(v[1]-smallend[1]
                         +(v[2]-smallend[2])*length(1));
#endif
    return result;
}

class BoxCommHelper
{
public:

    explicit BoxCommHelper (const Box& bx, int* p_ = 0);

    int* data () const& { return p; }
    int* data () && = delete;

    Box make_box () const {
	return Box(IntVect(p), IntVect(p+BL_SPACEDIM), IntVect(p+2*BL_SPACEDIM));
    }

    static int size () { return 3*BL_SPACEDIM; }

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

void AllGatherBoxes (Array<Box>& bxs);

}

#endif /*BL_BOX_H*/
