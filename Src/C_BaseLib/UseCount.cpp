//
// $Id: UseCount.cpp,v 1.1 2001-07-19 20:02:47 lijewski Exp $
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

bool
UseCount::unique () const
{
    return *cnt == 1;
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

int
UseCount::linkCount () const
{
    return *cnt;
}
