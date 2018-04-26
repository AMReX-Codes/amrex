
#include <iostream>
#include <limits>

#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {

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
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        shiftHalf(i,nz[i]);
   return *this;
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
bdryLo (const Box& b, int  dir, int len)
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

Vector<int> SerializeBox(const Box &b)
{
  int count(0);
  Vector<int> retArray(AMREX_SPACEDIM * 3);
  for(int i(0); i < AMREX_SPACEDIM; ++i) {
    retArray[count] = b.smallEnd(i);
    ++count;
  }
  for(int i(0); i < AMREX_SPACEDIM; ++i) {
    retArray[count] = b.bigEnd(i);
    ++count;
  }
  IntVect ivType(b.type());
  for(int i(0); i < AMREX_SPACEDIM; ++i) {
    retArray[count] = ivType[i];
    ++count;
  }
  return retArray;
}

Box UnSerializeBox(const Vector<int> &serarray)
{
  BL_ASSERT(serarray.size() == (3 * AMREX_SPACEDIM));
  const int *iptr = serarray.dataPtr();
  return Box(IntVect(iptr),
             IntVect(iptr + AMREX_SPACEDIM),
	     IndexType(IntVect(iptr + (2 * AMREX_SPACEDIM))));
}

BoxCommHelper::BoxCommHelper (const Box& bx, int* p_)
    : p(p_)
{
    if (p == 0) {
	v.resize(3*AMREX_SPACEDIM);
	p = &v[0];
    }

    AMREX_D_EXPR(p[0]               = bx.smallend[0],
	   p[1]               = bx.smallend[1],
	   p[2]               = bx.smallend[2]);
    AMREX_D_EXPR(p[0+AMREX_SPACEDIM]   = bx.bigend[0],
	   p[1+AMREX_SPACEDIM]   = bx.bigend[1],
	   p[2+AMREX_SPACEDIM]   = bx.bigend[2]);
    const IntVect& typ = bx.btype.ixType();
    AMREX_D_EXPR(p[0+AMREX_SPACEDIM*2] = typ[0],
	   p[1+AMREX_SPACEDIM*2] = typ[1],
	   p[2+AMREX_SPACEDIM*2] = typ[2]);
}

void
AllGatherBoxes (Vector<Box>& bxs)
{
#ifdef BL_USE_MPI
    // cell centered boxes only!
    const auto szof_bx = Box::linearSize();

    const long count = bxs.size() * static_cast<long>(szof_bx);
    const auto& countvec = ParallelDescriptor::Gather(count, ParallelDescriptor::IOProcessorNumber());
    
    long count_tot = 0L;
    Vector<long> offset(countvec.size(),0L);
    if (ParallelDescriptor::IOProcessor())
    {
        count_tot = countvec[0];
        for (int i = 1, N = offset.size(); i < N; ++i) {
            offset[i] = offset[i-1] + countvec[i-1];
            count_tot += countvec[i];
        }
    }

    ParallelDescriptor::Bcast(&count_tot, 1, ParallelDescriptor::IOProcessorNumber());

    if (count_tot == 0) return;

    Vector<char> send_buffer(count);
    char* psend = (count > 0) ? send_buffer.data() : nullptr;
    char* p = psend;
    for (const auto& b : bxs) {
        b.linearOut(p);
        p += szof_bx;
    }

    Vector<char> recv_buffer(count_tot);
    ParallelDescriptor::Gatherv(psend, count, recv_buffer.data(), countvec, offset, ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::Bcast(recv_buffer.data(), count_tot, ParallelDescriptor::IOProcessorNumber());

    const long nboxes_tot = count_tot/szof_bx;
    bxs.resize(nboxes_tot);

    p = recv_buffer.data();
    for (auto& b : bxs) {
        b.linearIn(p);
        p += szof_bx;
    }
#endif
}

}
