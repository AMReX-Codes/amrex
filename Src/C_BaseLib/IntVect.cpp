//BL_COPYRIGHT_NOTICE

//
// $Id: IntVect.cpp,v 1.4 1998-04-16 20:31:47 lijewski Exp $
//

#include <Assert.H>
#include <BoxLib.H>
#include <Misc.H>
#include <IntVect.H>
#include <IndexType.H>

const IntVect&
IntVect::TheUnitVector ()
{
    static const IntVect Unit(D_DECL(1,1,1));
    return Unit;
}

const IntVect&
IntVect::TheZeroVector ()
{
    static const IntVect Zero(D_DECL(0,0,0));
    return Zero;
}

const IntVect&
IntVect::TheNodeVector ()
{
    static const IntVect Node(D_DECL(IndexType::NODE,IndexType::NODE,IndexType::NODE));
    return Node;
}

const IntVect&
IntVect::TheCellVector ()
{
    static const IntVect Cell(D_DECL(IndexType::CELL,IndexType::CELL,IndexType::CELL));
    return Cell;
}

//
// Returns IntVect which is the componentwise integer projection
// of IntVect p1 by IntVect p2.
//

ostream&
operator<< (ostream&       os,
            const IntVect& p)
{
    os << D_TERM( '(' << p[0] , <<
                  ',' << p[1] , <<
                  ',' << p[2])  << ')';
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

istream&
operator>> (istream& is,
            IntVect& p)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> p[0],
               is.ignore(BL_IGNORE_MAX, ',') >> p[1],
               is.ignore(BL_IGNORE_MAX, ',') >> p[2]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
        D_EXPR(is.ignore(BL_IGNORE_MAX, '<') >> p[0],
               is.ignore(BL_IGNORE_MAX, ',') >> p[1],
               is.ignore(BL_IGNORE_MAX, ',') >> p[2]);
        is.ignore(BL_IGNORE_MAX, '>');
    }
    else
        BoxLib::Error("operator>>(istream&,IntVect&): expected \'(\' or \'<\'");

    if (is.fail())
        BoxLib::Error("operator>>(istream&,IntVect&) failed");

    return is;
}

void
IntVect::printOn (ostream& os) const
{
    os << "IntVect: " << *this << '\n';
}

void
IntVect::dumpOn (ostream& os) const
{
    os << "IntVect(" << BoxLib::version << ")= " << *this << '\n';
}
