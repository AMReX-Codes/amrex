//
// $Id: IndexType.cpp,v 1.12 2001-07-31 22:43:18 lijewski Exp $
//
#include <iostream>
#include <iomanip>

#include <IndexType.H>

IndexType::IndexType (const IntVect& iv)
{
    itype = D_TERM((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2));
}

IndexType::IndexType (D_DECL(CellIndex i, CellIndex j, CellIndex k))
{
    itype = D_TERM(i, | (j<<1), | (k<<2));
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
