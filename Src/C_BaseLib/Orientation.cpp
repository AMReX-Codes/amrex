//BL_COPYRIGHT_NOTICE

//
// $Id: Orientation.cpp,v 1.3 2000-04-24 17:52:36 car Exp $
//

#include <BoxLib.H>
#include <Orientation.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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

#ifdef BL_NAMESPACE
}
#endif

