//BL_COPYRIGHT_NOTICE

//
// $Id: ErrorList.cpp,v 1.5 1997-12-11 23:27:50 lijewski Exp $
//

#include <ErrorList.H>

ErrorList::~ErrorList ()
{
   clear();
}

void
ErrorList::clear ()
{
    for (ListIterator<ErrorRec*> li(lst); li; ++li)
    {
        delete lst[li];
        lst[li] = 0;
    }
}

const ErrorRec*
ErrorList::operator[] (int k) const
{
    assert(k < length());

    ListIterator<ErrorRec*> li(lst);
    
    while (li && k > 0)
    {
        ++li;
        --k;
    }
    return lst[li];
}

void
ErrorList::add (const aString&      name,
                int                 nextra,
                ErrorRec::ErrorType typ,
                ErrorRec::ErrorFunc func)
{
    ErrorRec* er = new ErrorRec(name,nextra,typ,func);
    //
    // Keep list in order of definition, append().
    //
    lst.append(er);
}

static const char* err_name[] = { "Richardson", "Special" };

ostream&
operator << (ostream &os, const ErrorList& elst)
{
    for (ListIterator<ErrorRec*> li(elst.lst); li; ++li)
    {
        os << li()->name()
           << ' '
           << li()->nGrow()
           << ' '
           << err_name[li()->errType()]
           << '\n';
    }
    return os;
}
