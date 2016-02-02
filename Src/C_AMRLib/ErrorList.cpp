
#include <iostream>
#include <BLassert.H>
#include <ErrorList.H>
#include <SPACE.H>

ErrorRec::ErrorFunc::ErrorFunc ()
    :
    m_func(0),
    m_func3D(0)
{}

ErrorRec::ErrorFunc::ErrorFunc (ErrorFuncDefault inFunc)
    :
    m_func(inFunc),
    m_func3D(0)
{}

ErrorRec::ErrorFunc::ErrorFunc (ErrorFunc3DDefault inFunc)
    :
    m_func(0),
    m_func3D(inFunc)
{}

ErrorRec::ErrorFunc*
ErrorRec::ErrorFunc::clone () const
{
    return new ErrorFunc(*this);
}

ErrorRec::ErrorFunc::~ErrorFunc () {}

void
ErrorRec::ErrorFunc::operator () (int* tag, D_DECL(const int&tlo0,const int&tlo1,const int&tlo2), 
                                  D_DECL(const int&thi0,const int&thi1,const int&thi2), 
                                  const int* tagval, const int* clearval,
                                  Real* data, D_DECL(const int&dlo0,const int&dlo1,const int&dlo2), 
                                  D_DECL(const int&dhi0,const int&dhi1,const int&dhi2), 
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const
{
    BL_ASSERT(m_func != 0);

    m_func(tag,D_DECL(tlo0,tlo1,tlo2),D_DECL(thi0,thi1,thi2),
           tagval,clearval,data,D_DECL(dlo0,dlo1,dlo2),D_DECL(dhi0,dhi1,dhi2),lo,hi,nvar,
           domain_lo,domain_hi,dx,xlo,prob_lo,time,level);
}

void
ErrorRec::ErrorFunc::operator () (int* tag, const int* tlo, const int* thi, 
                                  const int* tagval, const int* clearval,
                                  Real* data, const int* dlo, const int* dhi,
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const
{
    BL_ASSERT(m_func3D != 0);

    m_func3D(tag,ARLIM_3D(tlo),ARLIM_3D(thi),
             tagval,clearval,data,ARLIM_3D(dlo),ARLIM_3D(dhi),ARLIM_3D(lo),ARLIM_3D(hi),nvar,
             ARLIM_3D(domain_lo),ARLIM_3D(domain_hi),ZFILL(dx),ZFILL(xlo),ZFILL(prob_lo),time,level);
}  


ErrorRec::ErrorFunc2::ErrorFunc2 ()
    :
    m_func(0)
{}

ErrorRec::ErrorFunc2::ErrorFunc2 (ErrorFunc2Default inFunc)
    :
    m_func(inFunc)
{}

ErrorRec::ErrorFunc2*
ErrorRec::ErrorFunc2::clone () const
{
    return new ErrorFunc2(*this);
}

ErrorRec::ErrorFunc2::~ErrorFunc2 () {}


void
ErrorRec::ErrorFunc2::operator () (int* tag, D_DECL(const int&tlo0,const int&tlo1,const int&tlo2), 
                                   D_DECL(const int&thi0,const int&thi1,const int&thi2), 
                                   const int* tagval, const int* clearval,
                                   Real* data, D_DECL(const int&dlo0,const int&dlo1,const int&dlo2), 
                                   D_DECL(const int&dhi0,const int&dhi1,const int&dhi2), 
                                   const int* lo, const int * hi, const int* nvar,
                                   const int* domain_lo, const int* domain_hi,
                                   const Real* dx, const int* level, const Real* avg) const
{
    BL_ASSERT(m_func != 0);

    m_func(tag,D_DECL(tlo0,tlo1,tlo2),D_DECL(thi0,thi1,thi2),
           tagval,clearval,data,D_DECL(dlo0,dlo1,dlo2),D_DECL(dhi0,dhi1,dhi2),lo,hi,nvar,
           domain_lo,domain_hi,dx,level,avg);
}


ErrorRec::ErrorRec (const std::string&          nm,
                    int                         ng,
                    ErrorRec::ErrorType         etyp,
                    const ErrorRec::ErrorFunc2& f2)
    :
    derive_name(nm),
    ngrow(ng),
    err_type(etyp),
    err_func(0),
    err_func2(f2.clone())
{}

ErrorRec::ErrorRec (const std::string&         nm,
                    int                        ng,
                    ErrorRec::ErrorType        etyp,
                    const ErrorRec::ErrorFunc& f)
    :
    derive_name(nm),
    ngrow(ng),
    err_type(etyp),
    err_func(f.clone()),
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

const ErrorRec::ErrorFunc&
ErrorRec::errFunc () const
{
    return *err_func;
}

const ErrorRec::ErrorFunc2&
ErrorRec::errFunc2() const
{
    return *err_func2;
}

ErrorRec::~ErrorRec()
{
    delete err_func;
    delete err_func2;
}

int
ErrorList::size () const
{
    return vec.size();
}

void
ErrorList::add (const std::string&         name,
                int                        nextra, 
                ErrorRec::ErrorType        typ,
                const ErrorRec::ErrorFunc& func)
{
    //
    // Keep list in order of definition, append().
    //
    int n = vec.size();
    vec.resize(n+1);
    vec.set(n,new ErrorRec(name, nextra, typ, func));
}

void
ErrorList::add (const std::string&          name,
                int                         nextra,
                ErrorRec::ErrorType         typ,
                const ErrorRec::ErrorFunc2& func2)
{
    //
    // Keep list in order of definition, append().
    //
    int n = vec.size();
    vec.resize(n+1);
    vec.set(n,new ErrorRec(name, nextra, typ, func2));
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
