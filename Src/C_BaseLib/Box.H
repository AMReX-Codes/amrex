//BL_COPYRIGHT_NOTICE

#ifndef BL_BOX_H
#define BL_BOX_H

//
// $Id: Box.H,v 1.3 1997-09-18 20:12:44 lijewski Exp $
//

#ifdef NODEV
#undef NODEV
#endif

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::ws;
#else
#include <iostream.h>
#endif

#include <Boolean.H>
#include <IntVect.H>
#include <Misc.H>
#include <IndexType.H>
#include <Orientation.H>
#include <SPACE.H>

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
public:

    /*@ManDoc: The default constructor.  For safety, the constructed Box is  
               invalid  and may be tested for validity with ok() or
               isValid().
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
    //@ManDoc: The copy constructor.
    //
    Box (const Box& b);
    //
    //@ManDoc: The assignment operator.
    //
    Box& operator= (const Box& b);
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
    bool isValid () const;
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
    //@ManDoc: Equivalent to b.shift(0,iv[0]).shiftHalf(1,iv[1]) ....
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

    /*@ManDoc: Convert the Box from the current type into the
               argument type.  This may change the Box coordinates.
               type CELL -> NODE : increase coordinate by one on high end
               type NODE -> CELL : reduce coordinate by one on high end
               other type mappings make no change.
    */
    Box& convert (int                  dir,
                  IndexType::CellIndex typ);
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

    /*@ManDoc: Return a Box with NODE based coordinates in direction dir
               that encloses Box b.  NOTE: equivalent to b.convert(dir,NODE)
               NOTE: error if b.type(dir) == NODE.
    */
    friend Box surroundingNodes (const Box& b,
                                 int        dir);

    /*@ManDoc: Return a Box with NODE based coordinates in all
               directions that encloses Box b.
    */
    friend Box surroundingNodes (const Box& b);
 
    /*@ManDoc: Return a Box with CELL based coordinates in 
               direction dir that is enclosed by b.
               NOTE: equivalent to b.convert(dir,CELL)
               NOTE: error if b.type(dir) == CELL.
    */
    friend Box enclosedCells (const Box& b,
                              int        dir);
 
    /*@ManDoc: Return a Box with CELL based coordinates in all
               directions that is enclosed by b.
    */
    friend Box enclosedCells (const Box& b);
 
    /*@ManDoc: Returns the NODE based Box (in direction dir) defining
               the low side of Box b.
    */
    friend Box bdryLo (const Box& b,
                       int        dir,
                       int        len=1);
 
    /*@ManDoc: Returns the NODE based Box (in direction dir) defining
               the high side of Box b.
    */
    friend Box bdryHi (const Box& b,
                       int        dir,
                       int        len=1);
 
    /*@ManDoc: Similar to bdryLo and bdryHi except that it operates on the
               given face of  box b.
    */
    friend Box bdryNode (const Box&         b,
                         const Orientation& face,
                         int                len=1);
 
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
    friend Box adjCellLo (const Box& b,
                          int        dir,
                          int        len=1);
 
    /*@ManDoc: Similar to adjCellLo but builds an adjacent Box on
               the high end.
    */
    friend Box adjCellHi (const Box& b,
                          int        dir,
                          int        len=1);
 
    /*@ManDoc: Similar to adjCellLo and adjCellHi except that it operates
               on given face.
    */
    friend Box adjCell (const Box&         b,
                        const Orientation& face,
                        int                len=1);

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

    /*@ManDoc: Modify Box to that of the minimum Box containing both
               the original Box and the argument.
               Both Boxes must have identical type.
    */
    friend Box minBox (const Box& b1,
                       const Box& b2);
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

    /*@ManDoc: Grow Box in all directions by given amount.
               NOTE: n\_cell negative shrinks the Box by that number of cells.
    */
    friend inline Box grow (const Box& b,
                            int        i);
    //
    //@ManDoc: Grow Box in each direction by specified amount.
    //
    Box& grow (const IntVect& v);
    //
    //@ManDoc: Grow Box in each direction by specified amount.
    //
    friend inline Box grow (const Box&     b,
                            const IntVect& v);

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
    friend Box refine (const Box& b,
                       int        refinement_ratio);
    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    Box& refine (const IntVect& refinement_ratio);

    /*@ManDoc: Refine Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo*ratio and
               hi <- (hi+1)*ratio - 1.
               NOTE: if type(dir) = NODE centered: lo <- lo*ratio and
               hi <- hi*ratio.
    */
    friend Box refine (const Box&     b,
                       const IntVect& refinement_ratio);

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
    friend Box coarsen (const Box& b,
                        int        refinement_ratio);
    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    Box& coarsen (const IntVect& refinement_ratio);

    /*@ManDoc: Coarsen Box by given (positive) refinement ratio.
               NOTE: if type(dir) = CELL centered: lo <- lo/ratio and
               hi <- hi/ratio.
               NOTE: if type(dir) = NODE centered: lo <- lo/ratio and
               hi <- hi/ratio + ((hi%ratio)==0 ? 0 : 1).
               That is, refinement of coarsened Box must contain
               the original Box.
    */
    friend Box coarsen (const Box&     b,
                        const IntVect& refinement_ratio);

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
    //
    //@ManDoc: Write an ASCII representation to the ostream.
    //
    friend ostream& operator<< (ostream&   os,
                                const Box& bx);
    //
    //@ManDoc: Read from istream.
    //
    friend istream& operator>> (istream& os,
                                Box&     bx);

    /*@ManDoc: Gives more detail than printOn.
               Useful for exiting due to an error.
    */
    void dumpOn (ostream& strm = cerr) const;

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

    IntVect   smallend;
    IntVect   bigend;
    IntVect   len;
    IndexType btype;

private:
    //
    // Sets the 'len' element of the Box.  Aborts on integer overflow.
    //
    void computeBoxLen ();
};

//
//@ManDoc: BOX is a typedef for Box provided for backward compatibility.
//
typedef Box BOX;

//
// Inlines.
//

inline
Box::Box (const Box& b)
    : smallend(b.smallend),
      bigend(b.bigend),
      btype(b.btype)
{
    D_EXPR(len[0] = b.len[0],
           len[1] = b.len[1],
           len[2] = b.len[2]);
}

inline
Box&
Box::operator= (const Box& b)
{
    smallend = b.smallend;
    bigend = b.bigend;
    btype = b.btype;
    D_EXPR(len[0] = b.len[0],
           len[1] = b.len[1],
           len[2] = b.len[2]);
    return *this;
}

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
const IntVect&
Box::length () const
{
    return len;
}

inline
int
Box::length (int dir) const
{
    return len[dir];
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
    int dir = face.coordDir();
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
Box::ok () const
{
    bool typ_ok = btype.ok();
    bool len_ok = D_TERM((len[0] == (bigend[0]-smallend[0]+1)),
                          && (len[1] == (bigend[1]-smallend[1]+1)),
                          && (len[2] == (bigend[2]-smallend[2]+1)));
    return (bigend >= smallend) && len_ok && typ_ok;
}

inline
bool
Box::isValid () const
{
    return ok();
}

inline
bool
Box::contains (const IntVect& p) const
{
    return p >= smallend && p <= bigend;
}

inline
bool
Box::sameType (const Box &b) const
{
    return btype == b.btype;
}

inline
bool
Box::contains (const Box& b) const
{
    assert(sameType(b));
    return b.smallend >= smallend && b.bigend <= bigend;
}

inline
bool
Box::sameSize (const Box& b) const
{
    assert(sameType(b));
    return D_TERM(len[0] == b.len[0],
                  && len[1]==b.len[1],
                  && len[2]==b.len[2]);
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
Box::volumeOK () const
{
    long ignore;
    return volumeOK(ignore);
}

inline
long
Box::index (const IntVect& v) const
{
    long result = v.vect[0]-smallend.vect[0];
#if   BL_SPACEDIM==2
    result += len[0]*(v.vect[1]-smallend.vect[1]);
#elif BL_SPACEDIM==3
    result += len[0]*(v.vect[1]-smallend.vect[1]
                      +(v.vect[2]-smallend.vect[2])*len[1]);
#endif
    return result;
}

inline
void
Box::computeBoxLen ()
{
    D_EXPR(len[0] = bigend[0]-smallend[0] + 1,
           len[1] = bigend[1]-smallend[1] + 1,
           len[2] = bigend[2]-smallend[2] + 1);
}

inline
Box&
Box::setSmall (const IntVect& sm)
{
    smallend = sm;
    computeBoxLen();
    return *this;
}

inline
Box&
Box::setSmall (int dir,
               int sm_index)
{
    smallend.setVal(dir,sm_index);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::setBig (const IntVect& bg)
{
    bigend = bg;
    computeBoxLen();
    return *this;
}

inline
Box&
Box::setBig (int dir,
             int bg_index)
{
    bigend.setVal(dir,bg_index);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::setRange (int dir,
               int sm_index,
               int n_cells)
{
    smallend.setVal(dir,sm_index);
    bigend.setVal(dir,sm_index+n_cells-1);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::shift (int dir,
            int nzones)
{
    smallend.shift(dir,nzones);
    bigend.shift(dir,nzones);
    return *this;
}

inline
Box&
Box::shift (const IntVect& iv)
{
    smallend.shift(iv);
    bigend.shift(iv);
    return *this;
}

inline
Box&
Box::convert (const IntVect& typ)
{
    assert(typ >= IntVect::TheZeroVector() && typ <= IntVect::TheUnitVector());
    IntVect shft(typ - btype.ixType());
    bigend += shft;
    btype = typ;
    computeBoxLen();
    return *this;
}

inline
Box&
Box::surroundingNodes (int dir)
{
    if (!(btype[dir]))
    {
        bigend.shift(dir,1);
        //
        // Set dir'th bit to 1 = IndexType::NODE.
        //
        btype.set(dir);
        computeBoxLen();
    }
    return *this;
}

inline
Box&
Box::enclosedCells (int dir)
{
    if (btype[dir])
    {
        bigend.shift(dir,-1);
        //
        // Set dir'th bit to 0 = IndexType::CELL.
        //
        btype.unset(dir);
        computeBoxLen();
    }
    return *this;
}

inline
Box&
Box::operator+= (const IntVect& v)
{
    smallend += v;
    bigend += v;
    return *this;
}

inline
Box
Box::operator+  (const IntVect& v) const
{
    IntVect small(smallend);
    small += v;
    IntVect big(bigend);
    big += v;
    return Box(small,big,btype);
}

inline
Box&
Box::operator-= (const IntVect& v)
{
    smallend -= v;
    bigend -= v;
    return *this;
}

inline
Box
Box::operator-  (const IntVect& v) const
{
    IntVect small = smallend;
    small -= v;
    IntVect big = bigend;
    big -= v;
    return Box(small,big,btype);
}

inline
Box&
Box::grow (int i)
{
    smallend.diagShift(-i);
    bigend.diagShift(i);
    computeBoxLen();
    return *this;
}

inline
Box
grow (const Box& b,
      int        i)
{
    IntVect small = diagShift(b.smallend,-i);
    IntVect big   = diagShift(b.bigend,i);
    return Box(small,big,b.btype);
}

inline
Box&
Box::grow (const IntVect& v)
{
    smallend -= v;
    bigend   += v;
    computeBoxLen();
    return *this;
}

inline
Box
grow (const Box&     b,
      const IntVect& v)
{
    IntVect small = b.smallend - v;
    IntVect big   = b.bigend   + v;
    return Box(small,big,b.btype);
}

inline
Box&
Box::grow (int idir,
           int n_cell)
{
    smallend.shift(idir, -n_cell);
    bigend.shift(idir, n_cell);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::growLo (int idir,
             int n_cell)
{
    smallend.shift(idir, -n_cell);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::growHi (int idir,
             int n_cell)
{
    bigend.shift(idir,n_cell);
    computeBoxLen();
    return *this;
}

inline
Box&
Box::grow (const Orientation& face,
           int                n_cell)
{
    int idir = face.coordDir();
    if (face.isLow())
        smallend.shift(idir, -n_cell);
    else
        bigend.shift(idir,n_cell);
    computeBoxLen();
    return *this;
}

#endif /*BL_BOX_H*/
