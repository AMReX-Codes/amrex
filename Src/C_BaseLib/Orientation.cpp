//BL_COPYRIGHT_NOTICE

//
// $Id: Orientation.cpp,v 1.1 1997-09-12 18:00:13 lijewski Exp $
//

#include <BoxLib.H>
#include <Orientation.H>
#include <Utility.H>

ostream&
operator<< (ostream&           os,
            const Orientation& o)
{
    os << '('<< o.val << ')' ;
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,Orientation&) failed");
    return os;
}

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
