//
// $Id: IndexType.cpp,v 1.9 2001-07-17 23:02:22 lijewski Exp $
//

#include <iostream>
#include <iomanip>

#include <IndexType.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

int
IndexType::mask (int k)
{
    return 1<<k;
}

IndexType::IndexType ()
    : itype(0)
{}

IndexType::IndexType (const IndexType& bt)
    : itype(bt.itype)
{}

IndexType& IndexType::operator= (const IndexType& bt)
{
    itype = bt.itype;
    return *this;
}

IndexType::IndexType (const IntVect& iv)
{
    itype = D_TERM((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2));
}

IndexType::IndexType (D_DECL(CellIndex i, CellIndex j, CellIndex k))
{
    itype = D_TERM(i, | (j<<1), | (k<<2));
}

void
IndexType::set (int dir)
{
    itype |= mask(dir);
}

void
IndexType::unset (int dir)
{
    itype &= ~mask(dir);
}

bool
IndexType::test (int dir) const
{
    return (itype & mask(dir)) != 0;
}

void
IndexType::setall ()
{
    itype = (1 << BL_SPACEDIM) - 1;
}

void
IndexType::clear ()
{
    itype = 0;
}

bool
IndexType::any () const
{
    return itype != 0;
}

bool
IndexType::ok () const
{
    return itype < (1 << BL_SPACEDIM);
}

void
IndexType::flip (int i)
{
    itype ^= mask(i);
}

bool
IndexType::operator== (const IndexType& t) const
{
    return t.itype == itype;
}

bool
IndexType::operator!= (const IndexType& t) const
{
    return t.itype != itype;
}

bool
IndexType::cellCentered () const
{
    return itype == 0;
}

bool
IndexType::nodeCentered () const
{
    return itype == (1<<BL_SPACEDIM)-1;
}

void
IndexType::setType (int       dir,
                    CellIndex t)
{
    t == CELL ? unset(dir) : set(dir);
}

IndexType::CellIndex
IndexType::ixType (int dir) const
{
    return (CellIndex) ((itype & (1<<dir)) >> dir);
}

int
IndexType::operator[] (int dir) const
{
    return test(dir);
}

IntVect
IndexType::ixType () const
{
    return IntVect(D_DECL(itype&1, (itype>>1)&1, (itype>>2)&1));
}

IndexType
IndexType::TheCellType ()
{
    static const IndexType Cell(D_DECL(IndexType::CELL,
                                       IndexType::CELL,
                                       IndexType::CELL));
    return Cell;
}

IndexType
IndexType::TheNodeType ()
{
    static const IndexType Node(D_DECL(IndexType::NODE,
                                       IndexType::NODE,
                                       IndexType::NODE));
    return Node;
}

std::ostream&
operator<< (std::ostream&    os,
            const IndexType& it)
{
    os << '('
       << D_TERM( (it.test(0)?'N':'C'),
                  << ',' << (it.test(1)?'N':'C'),
                  << ',' << (it.test(2)?'N':'C')) << ')' << std::flush;

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,IndexType&) failed");

    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            IndexType&    it)
{
    char D_DECL(t0,t1,t2);

    D_EXPR( is.ignore(BL_IGNORE_MAX, '(') >> t0,
            is.ignore(BL_IGNORE_MAX, ',') >> t1,
            is.ignore(BL_IGNORE_MAX, ',') >> t2);
    is.ignore(BL_IGNORE_MAX, ')');
    D_TERM(
        BL_ASSERT(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
        BL_ASSERT(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
        BL_ASSERT(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2));

    if (is.fail())
        BoxLib::Error("operator>>(ostream&,IndexType&) failed");

    return is;
}



#ifdef BL_NAMESPACE
}
#endif

