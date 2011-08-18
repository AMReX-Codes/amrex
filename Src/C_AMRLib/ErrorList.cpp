
#include <iostream>
#include <BLassert.H>
#include <ErrorList.H>

ErrorRec::ErrorRec (const std::string& nm,
                    int                ng,
                    ErrorType          etyp,
                    ErrorFunc2         f2)
    :
    derive_name(nm),
    ngrow(ng),
    err_type(etyp),
    err_func(0),
    err_func2(f2)
{}

ErrorRec::ErrorRec (const std::string& nm,
                    int                ng,
                    ErrorType          etyp,
                    ErrorFunc          f)
    :
    derive_name(nm),
    ngrow(ng),
    err_type(etyp),
    err_func(f),
    err_func2(0)
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

ErrorFunc2
ErrorRec::errFunc2() const
{
    return err_func2;
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

void
ErrorList::add (const std::string&  name,
                int                 nextra,
                ErrorRec::ErrorType typ,
                ErrorFunc2          func2)
{
    //
    // Keep list in order of definition, append().
    //
    vec.push_back(ErrorRec(name, nextra, typ, func2));
}

const ErrorRec&
ErrorList::operator[] (int k) const
{
    BL_ASSERT(k < size());

    return vec[k];
}

static const char* err_name[] = { "Special", "Standard", "UseAverage" };

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
