//
// $Id: Mask.cpp,v 1.1 1998-03-24 07:05:39 almgren Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <Looping.H>
#include <Mask.H>
#include <Utility.H>

const char NL = '\n';
const char SP = ' ';

// --------------------  CONSTRUCTORS AND I/O OPERATIONS
Mask::Mask(istream& is)
{
    readFrom(is);
}

ostream& operator << (ostream& os, const Mask& m)
{
    Box b(m.box());
    int ncomp = m.nComp();
    os << "(Mask: " << b << SP << ncomp << NL;
    IntVect sm = b.smallEnd();
    IntVect bg = b.bigEnd();
    for (IntVect p = sm; p <= bg; b.next(p)) {
        os << p;
        for( int k=0; k<ncomp; k++ ) os << "  " << m(p,k);
        os << NL;
    }
    os << ")\n";
    assert(os.good());
    return os;
}

istream& operator >> (istream& is, Mask& m)
{
    is.ignore(BL_IGNORE_MAX,':');
    Box b;
    int ncomp;
    is >> b >> ncomp;
    is.ignore(BL_IGNORE_MAX, '\n');
    m.resize(b,ncomp);
    IntVect sm = b.smallEnd();
    IntVect bg = b.bigEnd();
    IntVect q;
    for (IntVect p = sm; p <= bg; b.next(p)) {
        is >> q;
        assert( p == q);
        for( int k=0; k<ncomp; k++ ) is >> m(p,k);
        is.ignore(BL_IGNORE_MAX, '\n');
    }
    is.ignore(BL_IGNORE_MAX,'\n');
    assert(is.good());
    return is;
}

void
Mask::writeOn(ostream& os) const
{
    os << "(Mask: " << domain << SP << nvar << NL;
    const int* ptr = dataPtr();
    int len = domain.numPts();
    os.write( (char*) ptr, len*sizeof(int) );
    os << ")\n";
}

void
Mask::readFrom(istream& is)
{
    is.ignore(BL_IGNORE_MAX,':');
    Box b;
    int ncomp;
    is >> b >> ncomp;
    is.ignore(BL_IGNORE_MAX, '\n');
    resize(b,ncomp);
    int *ptr = dataPtr();
    int len = domain.numPts();
    is.read( (char*) ptr, len*sizeof(int) );
    is.ignore(BL_IGNORE_MAX, '\n');
}

// --------------------  LOGICAL AND  ---------------------------
Mask&
Mask::And(const Mask& src)
{
    ForAllThisXC(int,src) {
        thisR = (thisR ? srcR : 0);
    } EndForTX;
    return *this;
}

Mask&
Mask::And(const Mask& src, int srccomp, int destcomp, int numcomp)
{
    Box domain(box());
    ForAllThisBNNXC(int,domain,destcomp,numcomp,src,srccomp) {
        thisR = (thisR ? srcR : 0);
    } EndForTX;
    return *this;
}

Mask&
Mask::And(const Mask& src, const Box& subbox,
          int srccomp, int destcomp, int numcomp)
{
    ForAllThisBNNXC(int,subbox,destcomp,numcomp,src,srccomp) {
        thisR = (thisR ? srcR : 0);
    } EndForTX;
    return *this;
}

Mask&
Mask::And(const Mask& src, const Box& srcbox,
          const Box& destbox,
          int srccomp, int destcomp, int numcomp)
{
    ForAllThisBNNXCBN(int,destbox,destcomp,numcomp,src,srcbox,srccomp) {
        thisR = (thisR ? srcR : 0);
    } EndForTX;
    return *this;
}

// --------------------  LOGICAL OR  ---------------------------
Mask&
Mask::Or(const Mask& src)
{
    ForAllThisXC(int,src) {
        thisR = (thisR ? 1 : srcR);
    } EndForTX;
    return *this;
}

Mask&
Mask::Or(const Mask& src, int srccomp, int destcomp, int numcomp)
{
    Box domain(box());
    ForAllThisBNNXC(int,domain,destcomp,numcomp,src,srccomp) {
        thisR = (thisR ? 1 : srcR);
    } EndForTX;
    return *this;
}

Mask&
Mask::Or(const Mask& src, const Box& subbox,
         int srccomp, int destcomp, int numcomp)
{
    ForAllThisBNNXC(int,subbox,destcomp,numcomp,src,srccomp) {
        thisR = (thisR ? 1 : srcR);
    } EndForTX;
    return *this;
}

Mask&
Mask::Or(const Mask& src, const Box& srcbox,
         const Box& destbox,
         int srccomp, int destcomp, int numcomp)
{
    ForAllThisBNNXCBN(int,destbox,destcomp,numcomp,src,srcbox,srccomp) {
        thisR = (thisR ? 1 : srcR);
    } EndForTX;
    return *this;
}
