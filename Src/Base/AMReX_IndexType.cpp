
#include <iostream>
#include <iomanip>

#include <AMReX_IndexType.H>

namespace amrex {

IndexType
IndexType::TheCellType ()
{
    static const IndexType Cell(AMREX_D_DECL(IndexType::CELL,
                                             IndexType::CELL,
                                             IndexType::CELL));
    return Cell;
}

IndexType
IndexType::TheNodeType ()
{
    static const IndexType Node(AMREX_D_DECL(IndexType::NODE,
                                             IndexType::NODE,
                                             IndexType::NODE));
    return Node;
}

std::ostream&
operator<< (std::ostream&    os,
            const IndexType& it)
{
    os << '('
       << AMREX_D_TERM( (it.test(0)?'N':'C'),
                  << ',' << (it.test(1)?'N':'C'),
                  << ',' << (it.test(2)?'N':'C')) << ')' << std::flush;

    if (os.fail())
        amrex::Error("operator<<(ostream&,IndexType&) failed");

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
    char AMREX_D_DECL(t0,t1,t2);

    AMREX_D_EXPR( is.ignore(BL_IGNORE_MAX, '(') >> t0,
            is.ignore(BL_IGNORE_MAX, ',') >> t1,
            is.ignore(BL_IGNORE_MAX, ',') >> t2);
    is.ignore(BL_IGNORE_MAX, ')');
    AMREX_D_TERM(
        BL_ASSERT(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
        BL_ASSERT(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
        BL_ASSERT(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2));

    if (is.fail())
        amrex::Error("operator>>(ostream&,IndexType&) failed");

    return is;
}

}
