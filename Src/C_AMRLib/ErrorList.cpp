//
// $Id: ErrorList.cpp,v 1.10 2001-08-01 21:50:45 lijewski Exp $
//

#include <iostream>

#include <BLassert.H>
#include <ErrorList.H>

ErrorRec::ErrorRec (const std::string& nm,
                    int                ng,
                    ErrorType          etyp,
                    ErrorFunc          f)
    :
    derive_name(nm),
    ngrow(ng),
    err_func(f),
    err_type(etyp)
{}

const std::string&
ErrorRec::name () const
{
    return derive_name;
}

int
ErrorRec::nGrow () const
{
    return ngrow;
}

ErrorRec::ErrorType
ErrorRec::errType () const
{
    return err_type;
}

ErrorFunc
ErrorRec::errFunc () const
{
    return err_func;
}

int
ErrorList::size () const
{
    return vec.size();
}

void
ErrorList::add (const std::string&  name,
                int                 nextra, 
                ErrorRec::ErrorType typ,
                ErrorFunc           func)
{
    //
    // Keep list in order of definition, append().
    //
    vec.push_back(ErrorRec(name, nextra, typ, func));
}

const ErrorRec&
ErrorList::operator[] (int k) const
{
    BL_ASSERT(k < size());

    return vec[k];
}

static const char* err_name[] = { "Richardson", "Special" };

std::ostream&
operator << (std::ostream&    os,
             const ErrorList& elst)
{
    for (int i = 0; i < elst.size(); i++)
    {
        os << elst[i].name()
           << ' '
           << elst[i].nGrow()
           << ' '
           << err_name[elst[i].errType()]
           << '\n';
    }
    return os;
}
