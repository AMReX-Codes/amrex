//
// $Id: Box.cpp,v 1.15 2001-07-19 17:51:48 lijewski Exp $
//

#include <climits>

#include <BLassert.H>
#include <BoxLib.H>
#include <Box.H>

const Box&
Box::TheUnitBox ()
{
    static const Box Unit(IntVect::TheZeroVector(), IntVect::TheUnitVector());
    return Unit;
}

Box::Box (const Box& b)
    :
    smallend(b.smallend),
    bigend(b.bigend),
    btype(b.btype)
{ }

Box::Box ()
    :
    smallend(IntVect::TheUnitVector()),
    bigend(IntVect::TheZeroVector()),
    btype()
{ }

Box::Box (const IntVect& small,
          const int*     vec_len)
    :
    smallend(small),
    bigend(D_DECL(small[0]+vec_len[0]-1,
                  small[1]+vec_len[1]-1,
                  small[2]+vec_len[2]-1))
{ }

Box::Box (const IntVect&   small,
          const IntVect&   big,
          const IndexType& t)
    :
    smallend(small),
    bigend(big),
    btype(t)
{ }

Box::Box (const IntVect& small,
          const IntVect& big)
    :
    smallend(small),
    bigend(big)
{ }

Box::Box (const IntVect& small,
          const IntVect& big,
          const IntVect& typ)
    :
    smallend(small),
    bigend(big),
    btype(typ)
{
    BL_ASSERT(typ >= IntVect::TheZeroVector() && typ <= IntVect::TheUnitVector());
}

Box&
Box::operator= (const Box& b)
{
    smallend = b.smallend;
    bigend   = b.bigend;
    btype    = b.btype;
    return *this;
}

bool
Box::ok () const
{
    return (bigend >= smallend) && btype.ok();
}

const IntVect&
Box::smallEnd () const
{
    return smallend;
}

int
Box::smallEnd (int dir) const
{
    return smallend[dir];
}

const IntVect&
Box::bigEnd () const
{
    return bigend;
}

int
Box::bigEnd (int dir) const
{
    return bigend[dir];
}

IndexType
Box::ixType () const
{
    return btype;
}

IntVect
Box::type () const
{
    return btype.ixType();
}

IndexType::CellIndex
Box::type (int dir) const
{
    return btype.ixType(dir);
}

int
Box::length (int dir) const
{
    return bigend[dir] - smallend[dir] + 1;
}

const int*
Box::loVect () const
{
    return smallend.getVect();
}

const int*
Box::hiVect () const
{
    return bigend.getVect();
}

const int*
Box::getVect () const
{
    return smallend.getVect();
}

int
Box::operator[] (const Orientation& face) const
{
    const int dir = face.coordDir();
    return face.isLow() ? smallend[dir] : bigend[dir];
}

bool
Box::numPtsOK () const
{
    long ignore;
    return numPtsOK(ignore);
}

bool
Box::isEmpty () const
{
    return numPts() == 0;
}

bool
Box::contains (const IntVect& p) const
{
    return p >= smallend && p <= bigend;
}

bool
Box::sameType (const Box &b) const
{
    return btype == b.btype;
}

bool
Box::contains (const Box& b) const
{
    BL_ASSERT(sameType(b));
    return b.smallend >= smallend && b.bigend <= bigend;
}

bool
Box::sameSize (const Box& b) const
{
    BL_ASSERT(sameType(b));
    return D_TERM(length(0) == b.length(0),
                  && length(1)==b.length(1),
                  && length(2)==b.length(2));
}

bool
Box::operator== (const Box& b) const
{
    return smallend == b.smallend && bigend == b.bigend && b.btype == btype;
}

bool
Box::operator!= (const Box& b) const
{
    return !operator==(b);
}

bool
Box::cellCentered () const
{
    return !btype.any();
}

bool
Box::volumeOK () const
{
    long ignore;
    return volumeOK(ignore);
}

long
Box::index (const IntVect& v) const
{
    long result = v.vect[0]-smallend.vect[0];
#if   BL_SPACEDIM==2
    result += length(0)*(v.vect[1]-smallend.vect[1]);
#elif BL_SPACEDIM==3
    result += length(0)*(v.vect[1]-smallend.vect[1]
                      +(v.vect[2]-smallend.vect[2])*length(1));
#endif
    return result;
}

const IntVect&
Box::length () const
{
    D_EXPR(len[0] = bigend[0]-smallend[0] + 1,
           len[1] = bigend[1]-smallend[1] + 1,
           len[2] = bigend[2]-smallend[2] + 1);

    return len;
}

Box&
Box::setSmall (const IntVect& sm)
{
    smallend = sm;
    return *this;
}

Box&
Box::setSmall (int dir,
               int sm_index)
{
    smallend.setVal(dir,sm_index);
    return *this;
}

Box&
Box::setBig (const IntVect& bg)
{
    bigend = bg;
    return *this;
}

Box&
Box::setBig (int dir,
             int bg_index)
{
    bigend.setVal(dir,bg_index);
    return *this;
}

Box&
Box::setRange (int dir,
               int sm_index,
               int n_cells)
{
    smallend.setVal(dir,sm_index);
    bigend.setVal(dir,sm_index+n_cells-1);
    return *this;
}

Box&
Box::shift (int dir,
            int nzones)
{
    smallend.shift(dir,nzones);
    bigend.shift(dir,nzones);
    return *this;
}

Box&
Box::shift (const IntVect& iv)
{
    smallend.shift(iv);
    bigend.shift(iv);
    return *this;
}

Box&
Box::convert (const IntVect& typ)
{
    BL_ASSERT(typ >= IntVect::TheZeroVector() && typ <= IntVect::TheUnitVector());
    IntVect shft(typ - btype.ixType());
    bigend += shft;
    btype = IndexType(typ);
    return *this;
}

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
    }
    return *this;
}

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
    }
    return *this;
}

Box
BoxLib::surroundingNodes (const Box& b,
                          int        dir)
{
    Box bx(b);
    return bx.surroundingNodes(dir);
}

Box
BoxLib::surroundingNodes (const Box& b)
{
    Box bx(b);
    return bx.surroundingNodes();
}

Box
BoxLib::enclosedCells (const Box& b,
                       int        dir)
{
    Box bx(b);
    return bx.enclosedCells(dir);
}

Box
BoxLib::enclosedCells (const Box& b)
{
    Box bx(b);
    return bx.enclosedCells();
}

Box&
Box::operator+= (const IntVect& v)
{
    smallend += v;
    bigend   += v;
    return *this;
}

Box
Box::operator+  (const IntVect& v) const
{
    IntVect small(smallend);
    IntVect big(bigend);
    small += v;
    big   += v;
    return Box(small,big,btype);
}

Box&
Box::operator-= (const IntVect& v)
{
    smallend -= v;
    bigend   -= v;
    return *this;
}

Box
Box::operator-  (const IntVect& v) const
{
    IntVect small = smallend;
    IntVect big = bigend;
    small -= v;
    big   -= v;
    return Box(small,big,btype);
}

Box&
Box::grow (int i)
{
    smallend.diagShift(-i);
    bigend.diagShift(i);
    return *this;
}

Box
BoxLib::grow (const Box& b,
              int        i)
{
    IntVect small = diagShift(b.smallEnd(),-i);
    IntVect big   = diagShift(b.bigEnd(),i);
    return Box(small,big,b.ixType());
}

Box&
Box::grow (const IntVect& v)
{
    smallend -= v;
    bigend   += v;
    return *this;
}

Box
BoxLib::grow (const Box&     b,
              const IntVect& v)
{
    IntVect small = b.smallEnd() - v;
    IntVect big   = b.bigEnd()   + v;
    return Box(small,big,b.ixType());
}

Box&
Box::grow (int idir,
           int n_cell)
{
    smallend.shift(idir, -n_cell);
    bigend.shift(idir, n_cell);
    return *this;
}

Box&
Box::growLo (int idir,
             int n_cell)
{
    smallend.shift(idir, -n_cell);
    return *this;
}

Box&
Box::growHi (int idir,
             int n_cell)
{
    bigend.shift(idir,n_cell);
    return *this;
}

Box&
Box::grow (const Orientation& face,
           int                n_cell)
{
    int idir = face.coordDir();
    if (face.isLow())
        smallend.shift(idir, -n_cell);
    else
        bigend.shift(idir,n_cell);
    return *this;
}

Box
Box::operator& (const Box& rhs) const
{
    Box lhs(*this);
    return lhs &= rhs;
}

bool
Box::numPtsOK (long& N) const
{
    BL_ASSERT(ok());

    N = length(0);

    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        if (length(i) == 0)
        {
            N = 0;
            return true;
        }
        else if (N <= LONG_MAX/length(i))
        {
            N *= length(i);
        }
        else
        {
            //
            // The return value `N' will be undefined.
            //
            return false;
        }
    }

    return true;
}

long
Box::numPts () const
{
    long result;
    if (!numPtsOK(result))
        BoxLib::Error("Arithmetic overflow in Box::numPts()");
    return result;
}

bool
Box::volumeOK (long& N) const
{
    BL_ASSERT(ok());

    N = length(0)-btype[0];

    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        long diff = (length(i)-btype[i]);

        if (diff == 0)
        {
            N = 0;
            return true;
        }
        else if (N <= LONG_MAX/diff)
        {
            N *= diff;
        }
        else
        {
            //
            // The return value of N will be undefined.
            //
            return false;
        }
    }

    return true;
}

long
Box::volume () const
{
    long result;
    if (!volumeOK(result))
        BoxLib::Error("Arithmetic overflow in Box::volume()");
    return result;
}

Box&
Box::shiftHalf (int dir,
                int nzones)
{
    int nbit = (nzones<0 ? -nzones : nzones)%2;
    int nshift = nzones/2;
    unsigned int bit_dir = btype[dir];
    //
    // Toggle btyp bit if nzones is odd.
    //
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

//
// Boolean functions.
//

bool
Box::intersects (const Box& b) const
{
    BL_ASSERT(sameType(b));
    IntVect low(smallend);
    IntVect hi(bigend);
    low.max(b.smallend);
    hi.min(b.bigend);
    return low <= hi;
}

//
// Intersection functions.
//

Box&
Box::operator&= (const Box& b)
{
    BL_ASSERT(sameType(b));
    smallend.max(b.smallend);
    bigend.min(b.bigend);
    return *this;
}

Box&
Box::surroundingNodes ()
{
    for (int i = 0; i < BL_SPACEDIM; ++i)
        if ((btype[i] == 0))
            bigend.shift(i,1);
    btype.setall();
    return *this;
}

Box&
Box::enclosedCells ()
{
    for (int i = 0 ; i < BL_SPACEDIM; ++i)
        if (btype[i])
            bigend.shift(i,-1);
    btype.clear();
    return *this;
}

//
// Next:step through the rectangle with unit increment.
//

void
Box::next (IntVect& p) const
{
    BL_ASSERT(contains(p));

    p.shift(0,1);
#if BL_SPACEDIM==2
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
    }
#elif BL_SPACEDIM==3
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
        }
    }
#endif
}

//
// Scan point over region of object Box with a vector incrment
// point incrments by 0 direction portion of vector.  When end of
// Box is reached, an increment is made with the 1 direction portion
// of the vector, and the 0 direction scan is resumed.
// effectively, we are scanning a Box, whose length vector is the argument
// vector over the object Box.
// when scan is over, the argument point is over edge of object Box
// this is the signal that we can go no further.
//

void
Box::next (IntVect&   p,
           const int* shv) const
{
    BL_ASSERT(contains(p));

#if   BL_SPACEDIM==1
    p.shift(0,shv[0]);
#elif BL_SPACEDIM==2
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
    }
#elif BL_SPACEDIM==3
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if(!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
        }
    }
#endif
}

Box&
Box::refine (int refinement_ratio)
{
    IntVect shft(IntVect::TheUnitVector());
    shft -= btype.ixType();
    smallend.scale(refinement_ratio);
    bigend += shft; // Bigend does more than just multiply.
    bigend.scale(refinement_ratio);
    bigend -= shft;
    return *this;
}

Box&
Box::refine (const IntVect& refinement_ratio)
{
    IntVect shft(IntVect::TheUnitVector());
    shft -= btype.ixType();
    smallend *= refinement_ratio;
    bigend += shft;
    bigend   *= refinement_ratio;
    bigend -= shft;
    return *this;
}

//
// Define a macro which will compute an object's length vector from
// the smallend and bigend.  Do several versions according to dimension
// requires you to be in a member functions.

int
Box::longside () const
{
    int maxlen = length(0);
    for (int i = 1; i < BL_SPACEDIM; i++)
        if (length(i) > maxlen)
            maxlen = length(i);
    return maxlen;
}

int
Box::longside (int& dir) const
{
    int maxlen = length(0);
    dir = 0;
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        if (length(i) > maxlen)
        {
            maxlen = length(i);
            dir = i;
        }
    }
    return maxlen;
}

int
Box::shortside () const
{
    int minlen = length(0);
    for (int i = 1; i < BL_SPACEDIM; i++)
        if (length(i) < minlen)
            minlen = length(i);
    return minlen;
}

int
Box::shortside (int& dir) const
{
    int minlen = length(0);
    dir = 0;
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        if (length(i) < minlen)
        {
            minlen = length(i);
            dir = i;
        }
    }
    return minlen;
}

//
// Modified Box is low end, returned Box is high end.
// If CELL: chop_pnt included in high end.
// If NODE: chop_pnt included in both Boxes.
//

Box
Box::chop (int dir,
           int chop_pnt)
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

//
// Shift by half increments.
//

Box&
Box::shiftHalf (const IntVect& nz)
{
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
   {
       int nzones = nz[dir];
       int nbit = (nzones<0 ? -nzones : nzones)%2;
       int nshift = nzones/2;
       unsigned int bit_dir = btype[dir];
       //
       // Toggle btype bit if nzones is odd.
       //
       if (nbit)
           btype.flip(dir);
       if (nzones < 0)
           nshift -= (bit_dir ? nbit : 0);
       else
           nshift += (bit_dir ? 0 : nbit);
       smallend.shift(dir,nshift);
       bigend.shift(dir,nshift);
   }
   return *this;
}

Box&
Box::convert (int                  dir,
              IndexType::CellIndex typ)
{
   unsigned int bitval = btype[dir];
   int off = typ - bitval;
   bigend.shift(dir,off);
   //
   // Set dir'th bit to typ.
   //
   btype.setType(dir,typ);
   return *this;
}

Box&
Box::convert (IndexType t)
{
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
   {
      unsigned int typ = t[dir];
      unsigned int bitval = btype[dir];
      int off = typ - bitval;
      bigend.shift(dir,off);
      btype.setType(dir, (IndexType::CellIndex) typ);
   }
   return *this;
}

//
// Refinement functions.
//

Box
BoxLib::refine (const Box& b,
                int        refinement_ratio)
{
    IntVect small(b.smallEnd());
    small.scale(refinement_ratio);
    IntVect big(b.bigEnd());
    IntVect shft(IntVect::TheUnitVector());
    shft -= b.ixType().ixType();
    big += shft;        // Large end is not just multiply.
    big.scale(refinement_ratio);
    big -= shft;
    return Box(small,big,b.ixType());
}

Box
BoxLib::refine (const Box&     b,
                const IntVect& refinement_ratio)
{
    IntVect small(b.smallEnd());
    small *= refinement_ratio;
    IntVect big(b.bigEnd());
    IntVect shft(IntVect::TheUnitVector());
    shft -= b.ixType().ixType();
    big += shft;
    big *= refinement_ratio;
    big -= shft;
    return Box(small,big,b.ixType());
}

//
// Coarsening.
//

Box&
Box::coarsen (const IntVect& refinement_ratio)
{
    smallend.coarsen(refinement_ratio);

    if (btype.any())
    {
        IntVect off(IntVect::TheZeroVector());
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (btype[dir])
            {
                int b = bigend[dir];
                int r = refinement_ratio[dir];
                if (b%r)
                    off.setVal(dir,1);
            }
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
    {
        bigend.coarsen(refinement_ratio);
    }

    return *this;
}

Box
BoxLib::coarsen (const Box& b,
                 int        refinement_ratio)
{
    IntVect small(b.smallEnd());
    small.coarsen(refinement_ratio);
    IntVect big(b.bigEnd());

    if (b.ixType().any())
    {
        IntVect off(IntVect::TheZeroVector());
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (b.ixType()[dir])
                if (big[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.ixType());
}

Box&
Box::coarsen (int refinement_ratio)
{
    smallend.coarsen(refinement_ratio);
    if (btype.any())
    {
        IntVect off(IntVect::TheZeroVector());
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (btype[dir])
                if (bigend[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
    {
        bigend.coarsen(refinement_ratio);
    }

    return *this;
}

Box
BoxLib::coarsen (const Box&     b,
                 const IntVect& refinement_ratio)
{
    IntVect small(b.smallEnd());
    small.coarsen(refinement_ratio);
    IntVect big(b.bigEnd());

    if (b.ixType().any())
    {
        IntVect off(IntVect::TheZeroVector());
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (b.ixType()[dir])
                if (big[dir]%refinement_ratio[dir])
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.ixType());
}

//
// I/O functions.
//

std::ostream&
operator<< (std::ostream& os,
            const Box&    b)
{
    os << '('
       << b.smallEnd() << ' '
       << b.bigEnd()   << ' '
       << b.type()
       << ')';

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,Box&) failed");

    return os;
}

//
// Moved out of Utility.H
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            Box&          b)
{
    IntVect lo, hi, typ;

    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        is >> lo >> hi >> typ;

        is.ignore(BL_IGNORE_MAX,')');
    }
    else
    {
        BoxLib::Error("operator>>(istream&,Box&): expected \'(\'");
    }

    b = Box(lo,hi,typ);

    if (is.fail())
        BoxLib::Error("operator>>(istream&,Box&) failed");

    return is;
}

//
// Give smallest Box containing both Boxes.
//

Box&
Box::minBox (const Box &b)
{
    BL_ASSERT(b.ok() && ok());
    BL_ASSERT(sameType(b));
    smallend.min(b.smallend);
    bigend.max(b.bigend);
    return *this;
}

Box
BoxLib::minBox (const Box& b,
                const Box& o)
{
    BL_ASSERT(b.ok() && o.ok());
    BL_ASSERT(o.sameType(b));
    IntVect small = b.smallEnd();
    IntVect big   = b.bigEnd();
    small.min(o.smallEnd());
    big.max(o.bigEnd());
    return Box(small,big,b.ixType());
}

//
// Translation functions acting on relative vectors.
//

Box
BoxLib::bdryLo (const Box& b,
                int        dir,
                int        len)
{
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    int sm = low[dir];
    hi.setVal(dir,sm+len-1);
    //
    // set dir'th bit to 1 = IndexType::NODE.
    //
    IndexType typ(b.ixType());
    typ.set(dir);
    return Box(low,hi,typ);
}

Box
BoxLib::bdryHi (const Box& b,
                int        dir,
                int        len)
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

Box
BoxLib::bdryNode (const Box&         b,
                  const Orientation& face,
                  int                len)
{
    int dir = face.coordDir();
    IntVect low(b.smallEnd());
    IntVect hi(b.bigEnd());
    if (face.isLow())
    {
        int sm = low[dir];
        hi.setVal(dir,sm+len-1);
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

Box
BoxLib::adjCellLo (const Box& b,
                   int        dir,
                   int        len)
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

Box
BoxLib::adjCellHi (const Box& b,
                   int        dir,
                   int        len)
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

Box
BoxLib::adjCell (const Box&         b,
                 const Orientation& face,
                 int                len)
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
