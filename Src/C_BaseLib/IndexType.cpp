//BL_COPYRIGHT_NOTICE

//
// $Id: IndexType.cpp,v 1.2 1997-09-18 20:12:49 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <iomanip>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#include <iomanip.h>
#endif

#include <IndexType.H>
#include <Utility.H>

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

ostream&
operator<< (ostream&         os,
            const IndexType& it)
{
    os << '('
       << D_TERM( (it.test(0)?'N':'C'),
                  << ',' << (it.test(1)?'N':'C'),
                  << ',' << (it.test(2)?'N':'C')) << ')' << flush;

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,IndexType&) failed");

    return os;
}

istream&
operator>> (istream&   is,
            IndexType& it)
{
    char D_DECL(t0,t1,t2);

    D_EXPR( is.ignore(BL_IGNORE_MAX, '(') >> t0,
            is.ignore(BL_IGNORE_MAX, ',') >> t1,
            is.ignore(BL_IGNORE_MAX, ',') >> t2);
    is.ignore(BL_IGNORE_MAX, ')');
    D_TERM(
        assert(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
        assert(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
        assert(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2));

    if (is.fail())
        BoxLib::Error("operator>>(ostream&,IndexType&) failed");

    return is;
}



