//
// $Id: aString.cpp,v 1.18 2001-07-21 17:37:10 car Exp $
//

#include <string>
#include <BLassert.H>
#include <aString.H>
using std::string;

aString::aString ()
    :
    string()
{}

aString::aString (const char* str) 
    :
    string(str)
{}

aString::aString (const std::string& str) 
    :
    string(str)
{}

char
aString::operator[] (size_type index) const
{
    BL_ASSERT(index < length());
    const std::string& a= *this;
    return a[index];
}

char&
aString::operator[] (size_type index)
{
    BL_ASSERT(index < size());
    std::string& a = *this;
    return a[index];
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


