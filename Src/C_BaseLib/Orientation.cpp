//
// $Id: Orientation.cpp,v 1.9 2001-07-31 22:43:19 lijewski Exp $
//
#include <iostream>

#include <BoxLib.H>
#include <Orientation.H>

std::ostream&
operator<< (std::ostream&      os,
            const Orientation& o)
{
    os << '('<< int(o) << ')' ;
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,Orientation&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            Orientation&  o)
{
    char c;
    is >> c;

    if (c == '(')
    {
        is >> o.val;
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        BoxLib::Error("operator>>(istream&,Orientation&): expected \'(\'");
    }

    if (is.fail())
        BoxLib::Error("operator>>(ostream&,Orientation&) failed");

    return is;
}

