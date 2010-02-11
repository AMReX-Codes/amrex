//
// $Id: UseCount.cpp,v 1.3 2010-02-11 22:22:34 lijewski Exp $
//

#include <UseCount.H>

UseCount::UseCount ()
    :
    cnt(new unsigned int(1))
{}

void
UseCount::decrement ()
{
    if (unique())
    {
        delete cnt;
        cnt = 0;
    }
    else
    {
        --*cnt;
    }
}

UseCount&
UseCount::operator= (const UseCount& rhs)
{
    ++*rhs.cnt;
    decrement();
    cnt = rhs.cnt;
    return *this;
}

UseCount::~UseCount ()
{
    decrement();
}
