//BL_COPYRIGHT_NOTICE

//
// $Id: ErrorList.cpp,v 1.3 1997-11-26 20:41:43 lijewski Exp $
//

#include <ErrorList.H>

// ------------------------------------------------------------------
ErrorRec::ErrorRec(const aString &nm, int ng, ErrorType etyp, ErrorFunc f)
    : ngrow(ng), err_type(etyp), err_func(f), derive_name(nm)
{
}

// ------------------------------------------------------------------
ErrorRec::~ErrorRec()
{
}

// ------------------------------------------------------------------
ErrorList::ErrorList()
    : lst()
{}

// ------------------------------------------------------------------
ErrorList::~ErrorList()
{
   clear();
}

// ------------------------------------------------------------------
void
ErrorList::clear()
{
    ListIterator<ErrorRec*> li(lst);
    while(li) {
	delete lst[li];
	lst[li] = 0;
	++li;
    }
}

// ------------------------------------------------------------------
const ErrorRec*
ErrorList::operator[] (int k) const
{
    assert( k < length() );
    ListIterator<ErrorRec*> li(lst);
    while(li && k > 0) {
	++li;
	--k;
    }
    return lst[li];
}

// ------------------------------------------------------------------
void
ErrorList::add(const aString &name, int nextra, ErrorType typ,
	       ErrorFunc func)
{
    ErrorRec* er = new ErrorRec(name,nextra,typ,func);
    if (er == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    lst.append(er);
      // keep list in order of definition, append
}

static aString err_name[] = {"Richardson", "Special"};

// ------------------------------------------------------------------
ostream&
operator << (ostream &os, const ErrorList& elst)
{
    for (ListIterator<ErrorRec*> li(elst.lst); li; ++li)
    {
	os << li()->derive_name;
	os << ' ' << li()->ngrow << ' ' << err_name[li()->err_type];
	os << '\n';
	++li;
    }
    return os;
}
