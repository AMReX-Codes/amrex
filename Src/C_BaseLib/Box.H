#ifndef BL_BOX_H
#define BL_BOX_H
//
// $Id: Box.H,v 1.25 2001-10-10 20:12:44 car Exp $
//
#ifdef NODEV
#undef NODEV
#endif

#include <iosfwd>

#include <ccse-mpi.H>
#include <IntVect.H>
#include <IndexType.H>
#include <Orientation.H>
#include <SPACE.H>
//
// Forward declaration.
//
class Box;

namespace BoxLib
{
    /*@ManDoc: Grow Box in all directions by given amount.
               NOTE: n\_cell negative shrinks the Box by that number of cells.
    */
    Box grow (const Box& b, int i);
    //
    //@ManDoc: Grow Box in each direction by specified amount.
    //
    Box grow (const Box& b, const IntVect& v);

    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    Box coarsen (const Box& b, int refinement_ratio);

    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    Box coarsen (const Box& b, const IntVect& refinement_ratio);

    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    Box refine (const Box& b, int refinement_ratio);

    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    Box refine (const Box& b, const IntVect& refinement_ratio);

    /*@ManDoc: Return a Box with NODE based coordinates in direction dir
               that encloses Box b.  NOTE: equivalent to b.convert(dir,NODE)
               NOTE: error if b.type(dir) == NODE.
    */
    Box surroundingNodes (const Box& b, int dir);

    /*@ManDoc: Return a Box with NODE based coordinates in all
               directions that encloses Box b.
    */
    Box surroundingNodes (const Box& b);

    /*@ManDoc: Return a Box with CELL based coordinates in 
               direction dir that is enclosed by b.
               NOTE: equivalent to b.convert(dir,CELL)
               NOTE: error if b.type(dir) == CELL.
    */
    Box enclosedCells (const Box& b, int dir);

    /*@ManDoc: Return a Box with CELL based coordinates in all
               directions that is enclosed by b.
    */
    Box enclosedCells (const Box& b);

    /*@ManDoc: Returns the edge-centered Box (in direction dir) defining
               the low side of Box b.
    */
    Box bdryLo (const Box& b, int dir, int len=1);

    /*@ManDoc: Returns the edge-centered Box (in direction dir) defining
               the high side of Box b.
    */
    Box bdryHi (const Box& b, int dir, int len=1);

    /*@ManDoc: Similar to bdryLo and bdryHi except that it operates on the
               given face of  box b.
    */
    Box bdryNode (const Box& b, const Orientation& face, int len=1);

   /*@ManDoc: Return the cell centered Box of length len adjacent
               to b on the low end along the coordinate direction dir.
               The return Box is identical to b in the other directions.
               The return Box and b have an empty intersection.
               NOTE:  len >= 1
               NOTE:  Box retval = b.adjCellLo(b,dir,len)
               is equivalent to the following set of operations:
               Box retval(b);
               retval.convert(dir,Box::CELL);
               retval.setrange(dir,retval.smallEnd(dir)-len,len);
    */
    Box adjCellLo (const Box& b, int dir, int len=1);

    /*@ManDoc: Similar to adjCellLo but builds an adjacent Box on the high end.
    */
    Box adjCellHi (const Box& b, int dir, int len=1);

   /*@ManDoc: Similar to adjCellLo and adjCellHi; operates on given face.
    */
    Box adjCell (const Box& b, const Orientation& face, int len=1);

    /*@ManDoc: Modify Box to that of the minimum Box containing both
               the original Box and the argument.
               Both Boxes must have identical type.
    */
    Box minBox (const Box& b1, const Box& b2);
}

//
//@ManDoc: Write an ASCII representation to the ostream.
//
std::ostream& operator<< (std::ostream& os, const Box& bx);
//
//@ManDoc: Read from istream.
//
std::istream& operator>> (std::istream& os, Box& bx);

//
//@Man:
//@Memo: A Rectangular Domain on an Integer Lattice
/*@Doc: 

  A Box is an abstraction for defining discrete regions of
  SPACEDIM indexing space.  Boxes have an IndexType, which defines
  IndexType::CELL or IndexType::NODE based points for each direction
  and a low and high INTVECT which defines the lower and upper corners
  of the Box.  Boxes can exist in positive and negative indexing space.

  Box is a dimension dependent class, so SPACEDIM must be 
  defined as either 1, 2, or 3 when compiling.

  For historical reasons, BOX is synonomous with Box.  That usage is
  deprecated, however.
*/
class Box
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<Box>::type();
public:

    /*@ManDoc: The default constructor.  For safety, the constructed Box is  
               invalid  and may be tested for validity with ok().
    */
    Box ();
    //
    //@ManDoc: Construct cell-centered type Box.
    //
    Box (const IntVect& small,
         const IntVect& big);
    //
    //@ManDoc: Construct box with specified lengths.
    //
    Box (const IntVect& small,
         const int*     vec_len);
 
    /*@ManDoc: Construct Box with given type.  small and big are expected
               to be consistent with given type.
    */
    Box (const IntVect& small,
         const IntVect& big,
         const IntVect& typ);
    //
    //@ManDoc: Construct dimension specific Boxes.
    //
    Box (const IntVect&   small,
         const IntVect&   big,
         const IndexType& t);
    //
    //@ManDoc: Get the smallend of the box.
    //
    const IntVect& smallEnd () const;
    //
    //@ManDoc: Return the coordinate of the low end in the given direction.
    //
    int smallEnd (int dir) const;
    //
    //@ManDoc: Get the bigend.
    //
    const IntVect& bigEnd () const;
    //
    //@ManDoc: Return the coordinate of the high end in the given direction.
    //
    int bigEnd (int dir) const;
    //
    //@ManDoc: Return the indexing type.
    //
    IndexType ixType () const;
    //
    //@ManDoc: Return the indexing type.
    //
    IntVect type () const;
    //
    //@ManDoc: Return the indexing type in the specified direction.
    //
    IndexType::CellIndex type (int dir) const;
    //
    //@ManDoc: Return the length of the Box.
    //
    const IntVect& length () const;
    //
    //@ManDoc: Return the length of the Box in given direction.
    //
    int length (int dir) const;
 
    /*@ManDoc: Returns a constant pointer the array of low end coordinates.
               Useful for calls to FORTRAN.
    */
    const int* loVect () const;
 
    /*@ManDoc: Returns a constant pointer the array of high end coordinates.
               Useful for calls to FORTRAN.
    */
    const int* hiVect () const;
 
    /*@ManDoc: Returns a constant pointer to the array of coordinates in the
               Box.  Useful for fortran, but otherwise too dangerous for use.
    */
    const int* getVect () const;
    //
    //@ManDoc: Returns the coordinate normal to given face.
    //
    int operator[] (const Orientation& face) const;
    //
    //@ManDoc: Checks if it is a proper box (including a valid type).
    //
    bool isEmpty () const;
    //
    //@ManDoc: Checks if it is a proper Box (including a valid type).
    //
    bool ok () const;
    //
    //@ManDoc: Returns true if argument is contained within Box.
    //
    bool contains (const IntVect& p) const;
 
    /*@ManDoc: Returns true if argument is contained within Box.
               It is an error if the Boxes have different types.
    */
    bool contains (const Box& b) const;
 
    /*@ManDoc: Returns true if Boxes have non-null intersections.
               It is an error if the Boxes have different types.
    */
    bool intersects (const Box& b) const;
 
    /*@ManDoc: Returns true is Boxes same size, ie translates of each other,.
               It is an error if they have different types.
    */
    bool sameSize (const Box& b) const;
    //
    //@ManDoc: Returns true if Boxes have same type.
    //
    bool sameType (const Box &b) const;
    //
    //@ManDoc: Returns true if Boxes are identical (including type).
    //
    bool operator== (const Box& b) const;
    //
    //@ManDoc: Returns true if Boxes differ (including type).
    //
    bool operator!= (const Box& b) const;
    //
    //@ManDoc: Returns true if Box is cell-centered in all indexing directions.
    //
    bool cellCentered () const;
    //
    //@ManDoc: Is the number calculated by numPts() representable in a long?
    //
    bool numPtsOK () const;

    /*@ManDoc: Returns the number of points contained in the Box, else
               abort()s if the number cannot be represented in a long.
     */
    long numPts () const;
    //
    //@ManDoc: Is the number calculated by volume() representable in a long?
    //
    bool volumeOK () const;
 
    /*@ManDoc: Return the volume, in indexing space, of region enclosed by
               this Box.   This is identical to numPts() for CELL centered
               Box; i.e. consider a edge-based box in 2D which has two edges
               in it.  This box has one cell center, located between the
               two edges.  numPts() would return 2 as would volume().
               Otherwise, numPts() >= volume().  Abort()s if the number
               cannot be represented in a long.
    */
    long volume () const;
 
    /*@ManDoc: Returns length of longest side.  dir is modified to give
               direction with longest side: 0...SPACEDIM-1. Ignores type.
    */
    int longside (int& dir) const;
    //
    //@ManDoc: Returns length of longest side.  Ignores type.
    //
    int longside () const;
 
    /*@ManDoc: Returns length of shortest side.  dir is modified to give
               direction with shortest side: 0...SPACEDIM-1.  Ignores type.
    */
    int shortside (int& dir) const;
    //
    //@ManDoc: Returns length of shortest side.  Ignores type.
    //
    int shortside () const;
 
    /*@ManDoc: Returns offset of point from smallend; i.e.
               index(smallend) -> 0, bigend would return volume()-1.
               Is used in accessing FArrayBox.
    */
    long index (const IntVect& v) const;
    //
    //@ManDoc: Redefine the small end of the Box.
    //
    Box& setSmall (const IntVect& sm);
    //
    //@ManDoc: Redefine the small end of the Box.
    //
    Box& setSmall (int dir,
                   int sm_index);
    //
    //@ManDoc: Redefine the big end of the Box.
    //
    Box& setBig (const IntVect& bg);
    //
    //@ManDoc: Redefine the big end of the Box.
    //
    Box& setBig (int dir,
                 int bg_index);
 
    /*@ManDoc: Set the entire range in a given direction, starting at
               sm\_index with length n\_cells.  NOTE: This will yield an
               illegal Box if n\_cells <= 0.
    */
    Box& setRange (int dir,
                   int sm_index,
                   int n_cells = 1);
 
    /*@ManDoc: Shift this Box nzones indexing positions in coordinate
               direction dir.
    */
    Box& shift (int dir,
                int nzones);
    //
    //@ManDoc: Equivalent to b.shift(0,iv[0]).shift(1,iv[1]) ....
    //
    Box& shift (const IntVect& iv);
 
    /*@ManDoc: This member shifts the Box by "half" indicies, thereby
               converting the Box from type CELL to NODE and visa-versa.
               b.shiftHalf(0,1)  shifts b to the right by 1/2 cells.
               b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
               NOTE: If num\_halfs is EVEN the shift is num\_halfs/2 full
               zones and hence will not change the type.
               This is: b.shifthalf(4) == b.shift(2).
    */
    Box& shiftHalf (int dir,
                    int num_halfs);
    //
    //@ManDoc: Equivalent to b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...
    //
    Box& shiftHalf (const IntVect& iv);
 
    /*@ManDoc: Convert the Box from the current type into the
               argument type.  This may change the Box coordinates:
               type CELL -> NODE : increase coordinate by one on high end
               type NODE -> CELL : reduce coordinate by one on high end
               other type mappings make no change.
    */
    Box& convert (IndexType typ);

    /*@ManDoc: Convert the Box from the current type into the
               argument type.  This may change the Box coordinates:
               type CELL -> NODE : increase coordinate by one on high end
               type NODE -> CELL : reduce coordinate by one on high end
               other type mappings make no change.
    */
    Box& convert (const IntVect& typ);
    //
    //@ManDoc: Convert to NODE type in all directions.
    //
    Box& surroundingNodes ();
    //
    //@ManDoc: Convert to NODE type in given direction.
    //
    Box& surroundingNodes (int dir);
    //
    //@ManDoc: Convert to CELL type in all directions.
    //
    Box& enclosedCells ();
    //
    //@ManDoc: Convert to CELL type in given direction.
    //
    Box& enclosedCells (int dir);
 
    /*@ManDoc: Return Box that is intersection of this Box
               and argument.  The Boxes MUST be of same type.
    */
    Box operator& (const Box&) const;

    /*@ManDoc: Intersect this Box with its argument.
               The Boxes MUST be of the same type.
    */
    Box& operator&= (const Box&);
 
    /*@ManDoc: Modify Box to that of the minimum Box containing both
               the original Box and the argument.
               Both Boxes must have identical type.
    */
    Box& minBox (const Box &);
    //
    //@ManDoc: Shift Box (relative) by given IntVect.
    //
    Box& operator+= (const IntVect& v);
    //
    //@ManDoc: Shift Box (relative) by given IntVect.
    //
    Box  operator+  (const IntVect& v) const;
    //
    //@ManDoc: Shift Box (relative) by given IntVect.
    //
    Box& operator-= (const IntVect& v);
    //
    //@ManDoc: Shift Box (relative) by given IntVect.
    //
    Box  operator-  (const IntVect& v) const;

    /*@ManDoc: Chop the Box at the chop\_pnt in the dir direction
               returns one Box, modifies the object Box.
               The union of the two is the original Box.
               The modified Box is the low end, the returned Box
               is the high end.  If type(dir) = CELL, the Boxes are disjoint
               with the chop\_pnt included in the high end (new Box).
               It is an ERROR if chop\_pnt is the low end of the orig Box.
               If type(dir) = NODE, the chop\_pnt is included in both Boxes
               but is the only point in common.  It is also an error if the
               chop\_pnt is an end node of the Box.
    */
    Box chop (int dir,
              int chop_pnt);

    /*@ManDoc: Grow Box in all directions by given amount.
               NOTE: n\_cell negative shrinks the Box by that number of cells.
    */
    Box& grow (int i);
    //
    //@ManDoc: Grow Box in each direction by specified amount.
    //
    Box& grow (const IntVect& v);

    /*@ManDoc: Grow the Box on the low and high end by n\_cell cells
               in direction idir.
    */
    Box& grow (int idir,
               int n_cell);

    /*@ManDoc: Grow the Box on the low end by n\_cell cells in direction idir.
              NOTE: n\_cell negative shrinks the Box by that number of cells.

    */
    Box& growLo (int idir,
                 int n_cell=1);

    /*@ManDoc: Grow the Box on the high end by n\_cell cells in
               direction idir.  NOTE: n\_cell negative shrinks the Box by that
               number of cells.
    */
    Box& growHi (int idir,
                 int n_cell=1);
    //
    //@ManDoc: Grow in the direction of the given face.
    //
    Box& grow (const Orientation& face,
               int                n_cell=1);

    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    Box& refine (int refinement_ratio);

    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    Box& refine (const IntVect& refinement_ratio);

    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    Box& coarsen (int refinement_ratio);

    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    Box& coarsen (const IntVect& refinement_ratio);

    /*@ManDoc: Step through the rectangle.  It is a runtime error to give
               a point not inside rectangle.  Iteration may not be efficient.
    */
    void next (IntVect &) const;

    /*@ManDoc: Scan argument IntVect over object second arg is
               increment vector.  Runtime error if IntVect is not
               contained in object Box.  Iteration may not be efficient.
    */
    void next (IntVect&   p,
               const int* shv) const;

    /*@ManDoc: This static member function returns a constant reference to 
               an object of type Box representing the unit box in
               BL\_SPACEDIM-dimensional space.
      */
    static const Box& TheUnitBox ();

protected:
    //
    // A helper function for numPtsOK() and numPts().
    //
    bool numPtsOK (long& N) const;
    //
    // A helper function for volumeOK() and volume().
    //
    bool volumeOK (long& N) const;

    IntVect         smallend;
    IntVect         bigend;
    mutable IntVect len;
    IndexType       btype;
};

inline
const IntVect&
Box::smallEnd () const
{
    return smallend;
}

inline
int
Box::smallEnd (int dir) const
{
    return smallend[dir];
}

inline
const IntVect&
Box::bigEnd () const
{
    return bigend;
}

inline
int
Box::bigEnd (int dir) const
{
    return bigend[dir];
}

inline
IndexType
Box::ixType () const
{
    return btype;
}

inline
IntVect
Box::type () const
{
    return btype.ixType();
}

inline
IndexType::CellIndex
Box::type (int dir) const
{
    return btype.ixType(dir);
}

inline
int
Box::length (int dir) const
{
    return bigend[dir] - smallend[dir] + 1;
}

inline
const int*
Box::loVect () const
{
    return smallend.getVect();
}

inline
const int*
Box::hiVect () const
{
    return bigend.getVect();
}

inline
const int*
Box::getVect () const
{
    return smallend.getVect();
}

inline
int
Box::operator[] (const Orientation& face) const
{
    const int dir = face.coordDir();
    return face.isLow() ? smallend[dir] : bigend[dir];
}

inline
bool
Box::numPtsOK () const
{
    long ignore;
    return numPtsOK(ignore);
}

inline
bool
Box::isEmpty () const
{
    return numPts() == 0;
}

inline
bool
Box::sameType (const Box &b) const
{
    return btype == b.btype;
}

inline
bool
Box::contains (const IntVect& p) const
{
    return p >= smallend && p <= bigend;
}

inline
bool
Box::contains (const Box& b) const
{
    BL_ASSERT(sameType(b));
    return b.smallend >= smallend && b.bigend <= bigend;
}

inline
bool
Box::sameSize (const Box& b) const
{
    BL_ASSERT(sameType(b));
    return D_TERM(length(0) == b.length(0),
                  && length(1)==b.length(1),
                  && length(2)==b.length(2));
}

inline
bool
Box::operator== (const Box& b) const
{
    return smallend == b.smallend && bigend == b.bigend && b.btype == btype;
}

inline
bool
Box::operator!= (const Box& b) const
{
    return !operator==(b);
}

inline
bool
Box::cellCentered () const
{
    return !btype.any();
}

inline
bool
Box::ok () const
{
    return (bigend >= smallend) && btype.ok();
}

inline
const IntVect&
Box::length () const
{
    D_EXPR(len[0] = bigend[0]-smallend[0] + 1,
           len[1] = bigend[1]-smallend[1] + 1,
           len[2] = bigend[2]-smallend[2] + 1);

    return len;
}

#endif /*BL_BOX_H*/
