
#include <iostream>
#include <limits>

#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Print.H>

namespace amrex {

const Box&
Box::TheUnitBox ()
{
    static const Box Unit(IntVect::TheZeroVector(), IntVect::TheZeroVector());
    return Unit;
}

Box::Box ()
    :
    smallend(IntVect::TheUnitVector()),
    bigend(IntVect::TheZeroVector()),
    btype()
{}

Box::Box (const IntVect& small,
          const int*     vec_len)
    :
    smallend(small),
    bigend(AMREX_D_DECL(small[0]+vec_len[0]-1,
                  small[1]+vec_len[1]-1,
                  small[2]+vec_len[2]-1))
{}

Box::Box (const IntVect& small,
          const IntVect& big,
          IndexType      t)
    :
    smallend(small),
    bigend(big),
    btype(t)
{}

Box::Box (const IntVect& small,
          const IntVect& big)
    :
    smallend(small),
    bigend(big)
{}

Box::Box (const IntVect& small,
          const IntVect& big,
          const IntVect& typ)
    :
    smallend(small),
    bigend(big),
    btype(typ)
{
    BL_ASSERT(typ.allGE(IntVect::TheZeroVector()) && typ.allLE(IntVect::TheUnitVector()));
}

Box&
Box::convert (const IntVect& typ)
{
    BL_ASSERT(typ.allGE(IntVect::TheZeroVector()) && typ.allLE(IntVect::TheUnitVector()));
    IntVect shft(typ - btype.ixType());
    bigend += shft;
    btype = IndexType(typ);
    return *this;
}

Box&
Box::convert (IndexType t)
{
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
   {
      const unsigned int typ = t[dir];
      const unsigned int bitval = btype[dir];
      const int off = typ - bitval;
      bigend.shift(dir,off);
      btype.setType(dir, (IndexType::CellIndex) typ);
   }
   return *this;
}

Box
convert (const Box& b, const IntVect& typ)
{
    Box bx(b);
    bx.convert(typ);
    return bx;
}

Box
convert (const Box& b, const IndexType& t)
{
    Box bx(b);
    bx.convert(t);
    return bx;
}

Box
surroundingNodes (const Box& b,
                          int        dir)
{
    Box bx(b);
    bx.surroundingNodes(dir);
    return bx;
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

Box
surroundingNodes (const Box& b)
{
    Box bx(b);
    bx.surroundingNodes();
    return bx;
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

Box
enclosedCells (const Box& b,
                       int        dir)
{
    Box bx(b);
    bx.enclosedCells(dir);
    return bx;
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
enclosedCells (const Box& b)
{
    Box bx(b);
    bx.enclosedCells();
    return bx;
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

Box
grow (const Box& b,
      int        i)
{
    Box result = b;
    result.grow(i);
    return result;
}

Box
grow (const Box&     b,
      const IntVect& v)
{
    Box result = b;
    result.grow(v);
    return result;
}

Box
grow (const Box& b, int idir, int n_cell)
{
    Box result = b;
    result.grow(idir, n_cell);
    return result;
}

Box
growLo (const Box& b, int idir, int n_cell)
{
    Box result = b;
    result.growLo(idir, n_cell);
    return result;
}

Box
growHi (const Box& b, int idir, int n_cell)
{
    Box result = b;
    result.growHi(idir, n_cell);
    return result;
}

Box&
Box::grow (Orientation face,
           int         n_cell)
{
    int idir = face.coordDir();
    if (face.isLow())
        smallend.shift(idir, -n_cell);
    else
        bigend.shift(idir,n_cell);
    return *this;
}

long
Box::numPts () const
{
    return AMREX_D_TERM( static_cast<long>(length(0)), 
                  *static_cast<long>(length(1)),
                  *static_cast<long>(length(2)));
}

double
Box::d_numPts () const
{
    BL_ASSERT(ok());

    return AMREX_D_TERM(double(length(0)), *double(length(1)), *double(length(2)));
}

long
Box::volume () const
{
    return AMREX_D_TERM( static_cast<long>(length(0)-btype[0]), 
                  *static_cast<long>(length(1)-btype[1]),
                  *static_cast<long>(length(2)-btype[2]));
}

Box&
Box::shiftHalf (int dir,
                int nzones)
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

Box&
Box::shiftHalf (const IntVect& nz)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
        shiftHalf(i,nz[i]);
   return *this;
}

void
Box::next (IntVect& p) const
{
    BL_ASSERT(contains(p));

    ++p[0];

#if (BL_SPACEDIM >= 2)
    if (p[0] > bigend[0])
    {
	p[0] = smallend[0];
	++p[1];
#if (BL_SPACEDIM == 3)
	if (p[1] > bigend[1])
	{
	    p[1] = smallend[1];
	    ++p[2];
	}
#endif
    }
#endif
}

Box
refine (const Box& b,
                int        ref_ratio)
{
    Box result = b;
    result.refine(IntVect(AMREX_D_DECL(ref_ratio,ref_ratio,ref_ratio)));
    return result;
}

Box&
Box::refine (int ref_ratio)
{
    return this->refine(IntVect(AMREX_D_DECL(ref_ratio,ref_ratio,ref_ratio)));
}

Box
refine (const Box&     b,
                const IntVect& ref_ratio)
{
    Box result = b;
    result.refine(ref_ratio);
    return result;
}

Box&
Box::refine (const IntVect& ref_ratio)
{
    if (ref_ratio != 1) {
        IntVect shft(IntVect::TheUnitVector());
        shft -= btype.ixType();
        smallend *= ref_ratio;
        bigend += shft;
        bigend *= ref_ratio;
        bigend -= shft;
    }
    return *this;
}

Box
shift (const Box& b, int dir, int nzones)
{
    Box result = b;
    result.shift(dir, nzones);
    return result;
}

//
// Define a macro which will compute an object's length vector from
// the smallend and bigend.  Do several versions according to dimension
// requires you to be in a member functions.

int
Box::longside () const
{
    int ignore = 0;
    return longside(ignore);
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
    int ignore = 0;
    return shortside(ignore);
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

Box
coarsen (const Box& b,
                 int        ref_ratio)
{
    Box result = b;
    result.coarsen(IntVect(AMREX_D_DECL(ref_ratio,ref_ratio,ref_ratio)));
    return result;
}

Box&
Box::coarsen (int ref_ratio)
{
    return this->coarsen(IntVect(AMREX_D_DECL(ref_ratio,ref_ratio,ref_ratio)));
}

Box
coarsen (const Box&     b,
                 const IntVect& ref_ratio)
{
    Box result = b;
    result.coarsen(ref_ratio);
    return result;
}

Box&
Box::coarsen (const IntVect& ref_ratio)
{
    smallend.coarsen(ref_ratio);

    if (btype.any())
    {
        IntVect off(IntVect::TheZeroVector());
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
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

    return *this;
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
        amrex::Error("operator<<(ostream&,Box&) failed");

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
        is >> lo >> hi;
	is >> c;
	// Read an optional IndexType
	is.putback(c);
	if ( c == '(' )
	{
	    is >> typ;
	}
        is.ignore(BL_IGNORE_MAX,')');
    }
    else if (c == '<')
    {
	is.putback(c);
        is >> lo >> hi;
	is >> c;
	// Read an optional IndexType
	is.putback(c);
	if ( c == '<' )
	{
	    is >> typ;
	}
        //is.ignore(BL_IGNORE_MAX,'>');
    }
    else
    {
        amrex::Error("operator>>(istream&,Box&): expected \'(\'");
    }

    b = Box(lo,hi,typ);

    if (is.fail())
        amrex::Error("operator>>(istream&,Box&) failed");

    return is;
}

Box
minBox (const Box& b,
                const Box& o)
{
    Box result = b;
    result.minBox(o);
    return result;
}

Box&
Box::minBox (const Box &b)
{
// BoxArray may call this with not ok boxes.  BL_ASSERT(b.ok() && ok());
    BL_ASSERT(sameType(b));
    smallend.min(b.smallend);
    bigend.max(b.bigend);
    return *this;
}

Box
bdryLo (const Box& b,
                int        dir,
                int        len)
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

Box
bdryHi (const Box& b,
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
bdryNode (const Box&  b,
                  Orientation face,
                  int         len)
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

Box
adjCellLo (const Box& b,
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
adjCellHi (const Box& b,
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
adjCell (const Box&  b,
                 Orientation face,
                 int         len)
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

bool
Box::sameSize (const Box& b) const
{
    BL_ASSERT(sameType(b));
    return AMREX_D_TERM(length(0) == b.length(0),
                  && length(1)==b.length(1),
                  && length(2)==b.length(2));
}

int
Box::operator[] (Orientation face) const
{
    const int dir = face.coordDir();
    return face.isLow() ? smallend[dir] : bigend[dir];
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

bool
Box::isSquare () const
{
    const IntVect& size = this->size();
#if BL_SPACEDIM==1
    return false; // can't build a square in 1-D
#elif BL_SPACEDIM==2
    return (size[0] == size[1]);
#elif BL_SPACEDIM==3
    return (size[0] == size[1] && (size[1] == size[2]));
#endif
}

Array<int> SerializeBox(const Box &b)
{
  int count(0);
  Array<int> retArray(BL_SPACEDIM * 3);
  for(int i(0); i < BL_SPACEDIM; ++i) {
    retArray[count] = b.smallEnd(i);
    ++count;
  }
  for(int i(0); i < BL_SPACEDIM; ++i) {
    retArray[count] = b.bigEnd(i);
    ++count;
  }
  IntVect ivType(b.type());
  for(int i(0); i < BL_SPACEDIM; ++i) {
    retArray[count] = ivType[i];
    ++count;
  }
  return retArray;
}


int SerializeBoxSize() {
  return (BL_SPACEDIM * 3);
}


Box UnSerializeBox(const Array<int> &serarray)
{
  BL_ASSERT(serarray.size() == (3 * BL_SPACEDIM));
  const int *iptr = serarray.dataPtr();
  return Box(IntVect(iptr),
             IntVect(iptr + BL_SPACEDIM),
	     IndexType(IntVect(iptr + (2 * BL_SPACEDIM))));
}

BoxCommHelper::BoxCommHelper (const Box& bx, int* p_)
    : p(p_)
{
    if (p == 0) {
	v.resize(3*BL_SPACEDIM);
	p = &v[0];
    }

    AMREX_D_EXPR(p[0]               = bx.smallend[0],
	   p[1]               = bx.smallend[1],
	   p[2]               = bx.smallend[2]);
    AMREX_D_EXPR(p[0+BL_SPACEDIM]   = bx.bigend[0],
	   p[1+BL_SPACEDIM]   = bx.bigend[1],
	   p[2+BL_SPACEDIM]   = bx.bigend[2]);
    const IntVect& typ = bx.btype.ixType();
    AMREX_D_EXPR(p[0+BL_SPACEDIM*2] = typ[0],
	   p[1+BL_SPACEDIM*2] = typ[1],
	   p[2+BL_SPACEDIM*2] = typ[2]);
}

}
