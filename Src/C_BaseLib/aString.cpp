//
// $Id: aString.cpp,v 1.15 2001-07-20 17:01:45 car Exp $
//

#include <BLassert.H>
#include <aString.H>

aString::aString()
    : std::string()
{
}

aString::aString(const char* str) 
    : std::string(str)
{
}

aString::aString(const std::string& str) 
    : std::string(str)
{
}

aString::size_type
aString::length() const
{
    return size();
}

char
aString::operator[] (size_type index) const
{
    BL_ASSERT(index >=0 && index < length());
    return std::string::operator[](index);
}

char&
aString::operator[] (size_type index)
{
    BL_ASSERT(index >= 0 && index < size());
    return std::string::operator[](index);
}

bool
aString::isNull() const
{
    return size()==0;
}

std::istream&
aString::getline(std::istream& is)
{
    return std::getline(is, *this);
}

int
aString::toInteger () const
{
    return size() == 0 ? 0 : atoi(c_str());
}

