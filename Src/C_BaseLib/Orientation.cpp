//BL_COPYRIGHT_NOTICE

//
// $Id: Orientation.cpp,v 1.2 1997-12-17 23:05:20 lijewski Exp $
//

#include <BoxLib.H>
#include <Orientation.H>

ostream&
operator<< (ostream&           os,
            const Orientation& o)
{
    os << '('<< o.val << ')' ;
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,Orientation&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

istream&
operator>> (istream&     is,
            Orientation& o)
{
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        is.ignore(BL_IGNORE_MAX, '(');
        is >> o.val;
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
        BoxLib::Error("operator>>(istream&,Orientation&): expected \'(\'");

    if (is.fail())
        BoxLib::Error("operator>>(ostream&,Orientation&) failed");

    return is;
}
