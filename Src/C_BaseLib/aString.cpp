//
// $Id: aString.cpp,v 1.14 2001-07-20 04:46:22 car Exp $
//

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

int
aString::length() const
{
    return size();
}

char
aString::operator[] (int index) const
{
    BL_ASSERT(index >=0 && index < size());
    return std::string::operator[](index);
}

char&
aString::operator[] (int index)
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

