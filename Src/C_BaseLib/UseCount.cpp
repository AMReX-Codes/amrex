//
// $Id: UseCount.cpp,v 1.2 2001-07-31 22:43:19 lijewski Exp $
//

#include <UseCount.H>

UseCount::UseCount ()
    :
    cnt(new unsigned int(1))
{}

UseCount::UseCount (const UseCount& rhs)
    :
    cnt(rhs.cnt)
{
    ++*cnt;
}

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
