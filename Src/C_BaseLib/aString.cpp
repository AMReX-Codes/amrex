//
// $Id: aString.cpp,v 1.17 2001-07-20 23:06:18 lijewski Exp $
//

#include <BLassert.H>
#include <aString.H>

aString::aString ()
    :
    std::string()
{}

aString::aString (const char* str) 
    :
    std::string(str)
{}

aString::aString (const std::string& str) 
    :
    std::string(str)
{}

char
aString::operator[] (size_type index) const
{
    BL_ASSERT(index < length());
    return std::string::operator[](index);
}

char&
aString::operator[] (size_type index)
{
    BL_ASSERT(index < size());
    return std::string::operator[](index);
}

bool
aString::isNull () const
{
    return size()==0;
}

std::istream&
aString::getline (std::istream& is)
{
    return std::getline(is, *this);
}


