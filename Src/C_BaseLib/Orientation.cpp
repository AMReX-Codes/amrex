//
// $Id: Orientation.cpp,v 1.6 2001-07-19 16:57:34 lijewski Exp $
//

#include <BoxLib.H>
#include <Orientation.H>

Orientation::Orientation (int _val)
    : val(_val)
{}

Orientation::Orientation ()
    : val(-1)
{}

Orientation::Orientation (int  _dir,
                          Side _side)
    : val(BL_SPACEDIM*_side + _dir)
{
    BL_ASSERT(0 <= _dir && _dir < BL_SPACEDIM);
}

Orientation::Orientation (const Orientation& o)
    : val(o.val)
{}

Orientation&
Orientation::operator= (const Orientation& o)
{
    val = o.val;
    return *this;
}

bool
Orientation::operator== (const Orientation& o) const
{
    return val == o.val;
}

bool
Orientation::operator!= (const Orientation& o) const
{
    return val != o.val;
}

bool
Orientation::operator<  (const Orientation& o) const
{
    return val < o.val;
}

bool
Orientation::operator<= (const Orientation& o) const
{
    return val <= o.val;
}

bool
Orientation::operator>  (const Orientation& o) const
{
    return val > o.val;
}

bool
Orientation::operator>= (const Orientation& o) const
{
    return val >= o.val;
}

Orientation::operator int () const
{
    return val;
}

Orientation
Orientation::flip () const
{
    return Orientation(val < BL_SPACEDIM ? val+BL_SPACEDIM : val-BL_SPACEDIM);
}

int
Orientation::coordDir () const
{
    return val%BL_SPACEDIM;
}

Orientation::Side
Orientation::faceDir () const
{
    return Side(val/BL_SPACEDIM);
}

bool
Orientation::isLow () const
{
    return val < BL_SPACEDIM;
}

bool
Orientation::isHigh () const
{
    return val >= BL_SPACEDIM;
}

OrientationIter::OrientationIter (int _face)
    : face(_face)
{}

OrientationIter::OrientationIter ()
    : face(0)
{}

OrientationIter::OrientationIter (const Orientation& _face)
    : face(_face)
{}

bool
OrientationIter::ok () const
{
    return 0 <= face && face < 2*BL_SPACEDIM;
}

OrientationIter::OrientationIter (const OrientationIter& it)
{
    BL_ASSERT(it.ok());
    face = it.face;
}

OrientationIter&
OrientationIter::operator= (const OrientationIter& it)
{
    BL_ASSERT(it.ok());
    face = it.face;
    return *this;
}

void
OrientationIter::rewind ()
{
    face = 0;
}

Orientation
OrientationIter::operator() () const
{
    BL_ASSERT(ok());
    return Orientation(face);
}

OrientationIter::operator void* ()
{
    return 0 <= face && face < 2*BL_SPACEDIM ? this : 0;
}

OrientationIter&
OrientationIter::operator-- ()
{
    BL_ASSERT(ok());
    --face;
    return *this;
}

OrientationIter&
OrientationIter::operator++ ()
{
    BL_ASSERT(ok());
    ++face;
    return *this;
}

OrientationIter
OrientationIter::operator-- (int)
{
    BL_ASSERT(ok());
    return OrientationIter(face--);
}

OrientationIter
OrientationIter::operator++ (int)
{
    BL_ASSERT(ok());
    return OrientationIter(face++);
}

bool
OrientationIter::operator== (const OrientationIter& oi) const
{
    BL_ASSERT(ok() && oi.ok());
    return face == oi.face;
}

bool
OrientationIter::operator!= (const OrientationIter& oi) const
{
    BL_ASSERT(ok() && oi.ok());
    return face != oi.face;
}

std::ostream&
operator<< (std::ostream&      os,
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

std::istream&
operator>> (std::istream& is,
            Orientation&  o)
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

