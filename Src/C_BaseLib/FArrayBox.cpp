
#include <winstd.H>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>

#include <FArrayBox.H>
#include <FabConv.H>
#include <ParmParse.H>
#include <FabConv.H>
#include <FPC.H>

#include <BLassert.H>
#include <BoxLib.H>
#include <Looping.H>
#include <Utility.H>
#include <BL_CXX11.H>
#include <MemPool.H>

#if (__GNUC__ >= 6 || defined(BL_Darwin))
using std::isinf;
using std::isnan;
#endif

#if defined(DEBUG) || defined(BL_TESTING)
bool FArrayBox::do_initval = true;
bool FArrayBox::init_snan  = true;
#else
bool FArrayBox::do_initval = false;
bool FArrayBox::init_snan  = false;
#endif
Real FArrayBox::initval;

static const char sys_name[] = "IEEE";
//
// Default Ordering to Normal Order.
//
FABio::Ordering FArrayBox::ordering = FABio::FAB_NORMAL_ORDER;

//
// Our 8-bit FABio type.
//
class FABio_8bit
    :
    public FABio
{
public:
    virtual void read (std::istream& is,
                       FArrayBox&    fb) const BL_OVERRIDE;

    virtual void write (std::ostream&    os,
                        const FArrayBox& fb,
                        int              comp,
                        int              num_comp) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f,
		       int           nCompToSkip) const BL_OVERRIDE;
private:
    virtual void write_header (std::ostream&    os,
                               const FArrayBox& f,
                               int              nvar) const BL_OVERRIDE;
};

//
// Our ASCII FABio type.
//
class FABio_ascii
    :
    public FABio
{
public:
    virtual void read (std::istream&   is,
                       FArrayBox&      fb) const BL_OVERRIDE;

    virtual void write (std::ostream&    os,
                        const FArrayBox& fb,
                        int              comp,
                        int              num_comp) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f,
		       int           nCompToSkip) const BL_OVERRIDE;
private:
    virtual void write_header (std::ostream&    os,
                               const FArrayBox& f,
                               int              nvar) const BL_OVERRIDE;
};

//
// Our binary FABio type.
//
class FABio_binary
    :
    public FABio
{
public:
    FABio_binary (RealDescriptor* rd_);

    virtual void read (std::istream& is,
                       FArrayBox&    fb) const BL_OVERRIDE;

    virtual void write (std::ostream&    os,
                        const FArrayBox& fb,
                        int              comp,
                        int              num_comp) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f) const BL_OVERRIDE;

    virtual void skip (std::istream& is,
                       FArrayBox&    f,
		       int           nCompToSkip) const BL_OVERRIDE;

private:
    virtual void write_header (std::ostream&    os,
                               const FArrayBox& f,
                               int              nvar) const BL_OVERRIDE;

    CpClassPtr<RealDescriptor> rd;
};

//
// This isn't inlined as it's virtual.
//

FABio::~FABio () {}

void
FABio::write_header (std::ostream&    os,
                     const FArrayBox& f,
                     int              nvar) const
{
    BL_ASSERT(nvar <= f.nComp());
    BoxLib::StreamRetry sr(os, "FABio_write_header", 4);
    while(sr.TryOutput()) {
      os << f.box() << ' ' << nvar << '\n';
    }
}

FABio::Format FArrayBox::format;

FABio* FArrayBox::fabio = 0;

FArrayBox::FArrayBox ()
{
    if (fabio == 0) FArrayBox::Initialize();
}

FArrayBox::FArrayBox (const Box& b,
                      int        n,
		      bool       alloc,
		      bool       shared)
    :
    BaseFab<Real>(b,n,alloc,shared)
{
    if (fabio == 0) FArrayBox::Initialize();
    if (alloc) initVal();
}

FArrayBox&
FArrayBox::operator= (const Real& v)
{
    BaseFab<Real>::operator=(v);
    return *this;
}

void
FArrayBox::initVal ()
{
    if (init_snan) {
#ifdef BL_USE_DOUBLE
	array_init_snan(dataPtr(), truesize);
#endif
    } else if (do_initval) {
	setVal(initval);
    }
}

bool 
FArrayBox::contains_nan () const
{
#ifndef _CRAYC
    const Real* dp = dptr;
    for (int i = 0; i < numpts*nvar; i++)
        if (isnan(*dp++))
            return true;
#endif
    return false;
}

bool 
FArrayBox::contains_nan (const Box& bx, int scomp, int ncomp) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(scomp <  nComp());
    BL_ASSERT(ncomp <= nComp());
    BL_ASSERT(domain.contains(bx));

#ifndef _CRAYC
    for (int i = 0; i < ncomp; i++)
    {
        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (isnan(this->operator()(p,scomp+i)))
                return true;
        }
    }
#endif
    return false;
}

bool 
FArrayBox::contains_nan (IntVect& where) const
{
    return contains_nan(domain, 0, nComp(), where);
}

bool 
FArrayBox::contains_nan (const Box& bx, int scomp, int ncomp, IntVect& where) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(scomp <  nComp());
    BL_ASSERT(ncomp <= nComp());
    BL_ASSERT(domain.contains(bx));

#ifndef _CRAYC
    for (int i = 0; i < ncomp; i++)
    {
        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (isnan(this->operator()(p,scomp+i)))
            {
                where = p;

                return true;
            }
        }
    }
#endif
    return false;
}

bool 
FArrayBox::contains_inf () const
{
    const Real* dp = dptr;
    for (int i = 0; i < numpts*nvar; i++)
        if (isinf(*dp++))
            return true;
    return false;
}

bool 
FArrayBox::contains_inf (const Box& bx, int scomp, int ncomp) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(scomp <  nComp());
    BL_ASSERT(ncomp <= nComp());
    BL_ASSERT(domain.contains(bx));

    for (int i = 0; i < ncomp; i++)
    {
        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (isinf(this->operator()(p,scomp+i)))
                return true;
        }
    }
    return false;
}

bool 
FArrayBox::contains_inf (IntVect& where) const
{
    return contains_inf(domain,0,nComp(),where);
}

bool 
FArrayBox::contains_inf (const Box& bx, int scomp, int ncomp, IntVect& where) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(scomp <  nComp());
    BL_ASSERT(ncomp <= nComp());
    BL_ASSERT(domain.contains(bx));

    for (int i = 0; i < ncomp; i++)
    {
        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (isinf(this->operator()(p,scomp+i)))
            {
                where = p;

                return true;
            }
        }
    }
    return false;
}

void
FArrayBox::resize (const Box& b,
                   int        N)
{
    BaseFab<Real>::resize(b,N);
    initVal();
}

FABio::Format
FArrayBox::getFormat ()
{
    return format;
}

const FABio&
FArrayBox::getFABio ()
{
    return *fabio;
}

void
FArrayBox::setFABio (FABio* rd)
{
    BL_ASSERT(rd != 0);
    delete fabio;
    fabio = rd;
}

Real
FArrayBox::norm (int p,
                 int comp,
                 int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}

void
FArrayBox::writeOn (std::ostream& os) const
{
    writeOn(os,0,nComp());
}

void
FArrayBox::skipFAB (std::istream& is)
{
    int ignore = 0;
    skipFAB(is, ignore);
}

FArrayBox::~FArrayBox () {}

void
FArrayBox::setFormat (FABio::Format fmt)
{
    FABio* fio = 0;

    switch (fmt)
    {
    case FABio::FAB_ASCII:
        fio = new FABio_ascii;
        break;
    case FABio::FAB_8BIT:
        fio = new FABio_8bit;
        break;
    case FABio::FAB_NATIVE:
        fio = new FABio_binary(FPC::NativeRealDescriptor().clone());
        break;
    case FABio::FAB_IEEE:
        //BoxLib::Warning("FABio::FAB_IEEE has been deprecated");
        //
        // Fall through ...
        //
    case FABio::FAB_IEEE_32:
        fio = new FABio_binary(FPC::Ieee32NormalRealDescriptor().clone());
        break;
    case FABio::FAB_NATIVE_32:
        fio = new FABio_binary(FPC::Native32RealDescriptor().clone());
        break;
    default:
        std::cerr << "FArrayBox::setFormat(): Bad FABio::Format = " << fmt;
        BoxLib::Abort();
    }

    FArrayBox::format = fmt;

    setFABio(fio);
}

void
FArrayBox::setOrdering (FABio::Ordering ordering_)
{
    //BoxLib::Warning("FArrayBox::setOrdering() has been deprecated");
    ordering = ordering_;
}

FABio::Ordering
FArrayBox::getOrdering ()
{
    //BoxLib::Warning("FArrayBox::getOrdering() has been deprecated");
    return ordering;
}

void
FArrayBox::setPrecision (FABio::Precision)
{
    BoxLib::Warning("FArrayBox::setPrecision() has been deprecated");
}

FABio::Precision
FArrayBox::getPrecision ()
{
    BoxLib::Warning("FArrayBox::getPrecision() has been deprecated");
    return FABio::FAB_FLOAT;
}

bool
FArrayBox::set_do_initval (bool tf)
{
    bool o_tf = do_initval;
    do_initval = tf;
    return o_tf;
}

bool
FArrayBox::get_do_initval ()
{
    return do_initval;
}

Real
FArrayBox::set_initval (Real iv)
{
    Real o_iv = initval;
    initval = iv;
    return o_iv;
}

Real
FArrayBox::get_initval ()
{
    return initval;
}

void
FArrayBox::Initialize ()
{
    BL_ASSERT(fabio == 0);

    ParmParse pp("fab");

    std::string fmt;
    //
    // This block can legitimately set FAB output format.
    //
    if (pp.query("format", fmt))
    {
        FABio* fio = 0;

        if (fmt == "ASCII")
        {
            FArrayBox::format = FABio::FAB_ASCII;
            fio = new FABio_ascii;
        }
        else if (fmt == "8BIT")
        {
            FArrayBox::format = FABio::FAB_8BIT;
            fio = new FABio_8bit;
        }
        else if (fmt == "NATIVE")
        {
            FArrayBox::format = FABio::FAB_NATIVE;
            fio = new FABio_binary(FPC::NativeRealDescriptor().clone());
        }
        else if (fmt == "NATIVE_32")
        {
            FArrayBox::format = FABio::FAB_NATIVE_32;
            fio = new FABio_binary(FPC::Native32RealDescriptor().clone());
        }
        else if (fmt == "IEEE" || fmt == "IEEE32")
        {
            if (fmt == "IEEE")
            {
                FArrayBox::format = FABio::FAB_IEEE;
                //BoxLib::Warning("IEEE fmt in ParmParse files is deprecated");
            }
            else
            {
                FArrayBox::format = FABio::FAB_IEEE_32;
            }
            fio = new FABio_binary(FPC::Ieee32NormalRealDescriptor().clone());
        }
        else
        {
            std::cerr << "FArrayBox::init(): Bad FABio::Format = " << fmt;
            BoxLib::Abort();
        }

        setFABio(fio);
    }
    else
    {
        //
        // We default to "NATIVE" if nothing is specified.
        //
        FArrayBox::format = FABio::FAB_NATIVE;

        setFABio(new FABio_binary(FPC::NativeRealDescriptor().clone()));
    }
    //
    // This block sets ordering which doesn't affect output format.
    // It is only used when reading in an old FAB.
    //
    std::string ord;

    if (pp.query("ordering", ord))
    {
        if (ord == "NORMAL_ORDER")
            FArrayBox::setOrdering(FABio::FAB_NORMAL_ORDER);
        else if (ord == "REVERSE_ORDER")
            FArrayBox::setOrdering(FABio::FAB_REVERSE_ORDER);
        else if (ord == "REVERSE_ORDER_2")
            FArrayBox::setOrdering(FABio::FAB_REVERSE_ORDER_2);
        else
        {
            std::cerr << "FArrayBox::init(): Bad FABio::Ordering = " << ord;
            BoxLib::Abort();
        }
    }

    initval = std::numeric_limits<Real>::has_quiet_NaN
	    ? std::numeric_limits<Real>::quiet_NaN()
	    : std::numeric_limits<Real>::max();

    pp.query("initval",    initval);
    pp.query("do_initval", do_initval);
    pp.query("init_snan", init_snan);

    BoxLib::ExecOnFinalize(FArrayBox::Finalize);
}

void
FArrayBox::Finalize ()
{
    delete fabio;
    fabio = 0;
}

Real
FArrayBox::norm (const Box& subbox,
                 int        p,
                 int        comp,
                 int        ncomp) const
{
    BL_ASSERT(p >= 0);
    BL_ASSERT(comp >= 0 && comp+ncomp <= nComp());

    Real  nrm    = 0;
    Real* tmp    = 0;
    int   tmplen = 0;

    if (p == 0 || p == 1)
    {
        nrm = BaseFab<Real>::norm(subbox,p,comp,ncomp);
    }
    else if (p == 2)
    {
        ForAllThisCPencil(Real,subbox,comp,ncomp)
        {
            const Real* row = &thisR;
            if (tmp == 0)
            {
                tmp    = new Real[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = row[i]*row[i];
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += row[i]*row[i];
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
        nrm = std::sqrt(nrm);
    }
    else
    {
        Real pwr = Real(p);
        ForAllThisCPencil(Real,subbox,comp,ncomp)
        {
            const Real* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = std::pow(row[i],pwr);
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += std::pow(row[i],pwr);
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
        Real invpwr = 1.0/pwr;
        nrm = std::pow(nrm,invpwr);
    }

    delete [] tmp;

    return nrm;
}

//
// Copied from Utility.H.
//
#define BL_IGNORE_MAX 100000

//
// This is where lies all the smarts for recognizing FAB headers.
//

FABio*
FABio::read_header (std::istream& is,
                    FArrayBox&    f)
{
    int nvar;
    Box bx;
    FABio* fio = 0;
    RealDescriptor* rd = 0;
    char c;

    is >> c;
    if (c != 'F') BoxLib::Error("FABio::read_header(): expected \'F\'");
    is >> c;
    if (c != 'A') BoxLib::Error("FABio::read_header(): expected \'A\'");
    is >> c;
    if (c != 'B') BoxLib::Error("FABio::read_header(): expected \'B\'");

    is >> c;
    if (c == ':')
    {
        //
        // The "old" FAB format.
        //
        int typ_in, wrd_in;
        is >> typ_in;
        is >> wrd_in;

        char machine[128];
        is >> machine;
        is >> bx;
        is >> nvar;
        //
        // Set the FArrayBox to the appropriate size.
        //
        f.resize(bx,nvar);
        is.ignore(BL_IGNORE_MAX, '\n');
        switch (typ_in)
        {
        case FABio::FAB_ASCII: fio = new FABio_ascii; break;
        case FABio::FAB_8BIT:  fio = new FABio_8bit;  break;
        case FABio::FAB_NATIVE:
        case FABio::FAB_NATIVE_32:
        case FABio::FAB_IEEE:
            rd = RealDescriptor::newRealDescriptor(typ_in,
                                                   wrd_in,
                                                   machine,
                                                   FArrayBox::ordering);
            fio = new FABio_binary(rd);
            break;
        default:
            BoxLib::Error("FABio::read_header(): Unrecognized FABio header");
        }
    }
    else
    {
        //
        // The "new" FAB format.
        //
        is.putback(c);
        rd = new RealDescriptor;
        is >> *rd;
        is >> bx;
        is >> nvar;
        //
        // Set the FArrayBox to the appropriate size.
        //
        f.resize(bx,nvar);
        is.ignore(BL_IGNORE_MAX, '\n');
        fio = new FABio_binary(rd);
    }

    if (is.fail())
        BoxLib::Error("FABio::read_header() failed");

    return fio;
}


FABio*
FABio::read_header (std::istream& is,
                    FArrayBox&    f,
		    int           compIndex,
		    int&          nCompAvailable)
{
    int nvar;
    Box bx;
    FABio* fio = 0;
    RealDescriptor* rd = 0;
    char c;

    is >> c;
    if (c != 'F') BoxLib::Error("FABio::read_header(): expected \'F\'");
    is >> c;
    if (c != 'A') BoxLib::Error("FABio::read_header(): expected \'A\'");
    is >> c;
    if (c != 'B') BoxLib::Error("FABio::read_header(): expected \'B\'");

    is >> c;
    if (c == ':')
    {
        //
        // The "old" FAB format.
        //
        int typ_in, wrd_in;
        is >> typ_in;
        is >> wrd_in;

        char machine[128];
        is >> machine;
        is >> bx;
        is >> nvar;
	nCompAvailable = nvar;
	nvar = 1;    // make a single component fab
        //
        // Set the FArrayBox to the appropriate size.
        //
        f.resize(bx,nvar);
        is.ignore(BL_IGNORE_MAX, '\n');
        switch (typ_in)
        {
        case FABio::FAB_ASCII: fio = new FABio_ascii; break;
        case FABio::FAB_8BIT:  fio = new FABio_8bit;  break;
        case FABio::FAB_NATIVE:
        case FABio::FAB_NATIVE_32:
        case FABio::FAB_IEEE:
            rd = RealDescriptor::newRealDescriptor(typ_in,
                                                   wrd_in,
                                                   machine,
                                                   FArrayBox::ordering);
            fio = new FABio_binary(rd);
            break;
        default:
            BoxLib::Error("FABio::read_header(): Unrecognized FABio header");
        }
    }
    else
    {
        //
        // The "new" FAB format.
        //
        is.putback(c);
        rd = new RealDescriptor;
        is >> *rd;
        is >> bx;
        is >> nvar;
	nCompAvailable = nvar;
	nvar = 1;    // make a single component fab
        //
        // Set the FArrayBox to the appropriate size.
        //
        f.resize(bx,nvar);
        is.ignore(BL_IGNORE_MAX, '\n');
        fio = new FABio_binary(rd);
    }

    if (is.fail())
        BoxLib::Error("FABio::read_header() failed");

    return fio;
}

void
FArrayBox::writeOn (std::ostream& os,
                    int           comp,
                    int           num_comp) const
{
    BL_ASSERT(comp >= 0 && num_comp >= 1 && (comp+num_comp) <= nComp());
    fabio->write_header(os, *this, num_comp);
    os.flush();  // 2016-08-30: Titan requires this flush() (probably due to a bug).
    fabio->write(os, *this, comp, num_comp);
}

void
FArrayBox::readFrom (std::istream& is)
{
    FABio* fabrd = FABio::read_header(is, *this);
    fabrd->read(is, *this);
    delete fabrd;
}


int
FArrayBox::readFrom (std::istream& is, int compIndex)
{
    int nCompAvailable;
    FABio* fabrd = FABio::read_header(is, *this, compIndex, nCompAvailable);
    BL_ASSERT(compIndex >= 0 && compIndex < nCompAvailable);

    fabrd->skip(is, *this, compIndex);  // skip data up to the component we want
    fabrd->read(is, *this);
    int remainingComponents = nCompAvailable - compIndex - 1;
    fabrd->skip(is, *this, remainingComponents);  // skip to the end

    delete fabrd;
    return nCompAvailable;
}


Box
FArrayBox::skipFAB (std::istream& is,
                    int&          num_comp)
{
    FArrayBox f;
    FABio* fabrd = FABio::read_header(is, f);
    fabrd->skip(is, f);
    delete fabrd;
    num_comp = f.nComp();
    return f.box();
}

void
FABio_ascii::write (std::ostream&    os,
                    const FArrayBox& f,
                    int              comp,
                    int              num_comp) const
{
    BL_ASSERT(comp >= 0 && num_comp >= 1 && (comp+num_comp) <= f.nComp());

    const Box& bx = f.box();

    IntVect sm = bx.smallEnd();
    IntVect bg = bx.bigEnd();

    for (IntVect p = sm; p <= bg; bx.next(p))
    {
        os << p;
        for (int k=0; k < num_comp; k++)
            os << "  " << f(p,k+comp);
        os << '\n';
    }
    os << '\n';

    if (os.fail())
        BoxLib::Error("FABio_ascii::write() failed");
}

void
FABio_ascii::read (std::istream& is,
                   FArrayBox&    f) const
{
    const Box& bx = f.box();

    IntVect sm = bx.smallEnd();
    IntVect bg = bx.bigEnd();
    IntVect p, q;
    for (p = sm; p <= bg; bx.next(p))
    {
        is >> q;
        if (p != q)
        {
          std::cerr << "Error: read IntVect "
                    << q
                    << "  should be "
                    << p
                    << '\n';
          BoxLib::Error("FABio_ascii::read() bad IntVect");
        }
        for (int k = 0; k < f.nComp(); k++)
            is >> f(p, k);
    }

    if (is.fail())
        BoxLib::Error("FABio_ascii::read() failed");
}

void
FABio_ascii::skip (std::istream& is,
                   FArrayBox&    f) const
{
    FABio_ascii::read(is, f);
}

void
FABio_ascii::skip (std::istream& is,
                   FArrayBox&    f,
		   int           nCompToSkip) const
{
    BoxLib::Error("FABio_ascii::skip(..., int nCompToSkip) not implemented");
}

void
FABio_ascii::write_header (std::ostream&    os,
                           const FArrayBox& f,
                           int              nvar) const
{
    os << "FAB: "
       << FABio::FAB_ASCII
       << ' '
       << 0
       << ' '
       << sys_name
       << '\n';
    FABio::write_header(os, f, nvar);
}

void
FABio_8bit::write (std::ostream&    os,
                   const FArrayBox& f,
                   int              comp,
                   int              num_comp) const
{
    BL_ASSERT(comp >= 0 && num_comp >= 1 && (comp+num_comp) <= f.nComp());

    const Real eps = Real(1.0e-8); // FIXME - whats a better value?
    const long siz = f.box().numPts();

    unsigned char* c = new unsigned char[siz];

    for (int k = 0; k < num_comp; k++)
    {
        const Real mn   = f.min(k+comp);
        const Real mx   = f.max(k+comp);
        const Real* dat = f.dataPtr(k+comp);
        Real rng = std::fabs(mx-mn);
        rng = (rng < eps) ? 0.0 : 255.0/(mx-mn);
        for (long i = 0; i < siz; i++)
        {
            Real v = rng*(dat[i]-mn);
            int iv = (int) v;
            c[i]   = (unsigned char) iv;
        }
        os << mn << "  " << mx << '\n' << siz << '\n';
        os.write((char*)c,siz);
    }

    delete [] c;

    if (os.fail())
        BoxLib::Error("FABio_8bit::write() failed");
}

void
FABio_8bit::read (std::istream& is,
                  FArrayBox&    f) const
{
    long siz         = f.box().numPts();
    unsigned char* c = new unsigned char[siz];

    Real mn, mx;
    for (int nbytes, k = 0; k < f.nComp(); k++)
    {
        is >> mn >> mx >> nbytes;
        BL_ASSERT(nbytes == siz);
        while (is.get() != '\n')
            ;
        is.read((char*)c,siz);
        Real* dat       = f.dataPtr(k);
        const Real rng  = (mx-mn)/255.0;
        for (long i = 0; i < siz; i++)
        {
            int iv = (int) c[i];
            Real v = (Real) iv;
            dat[i] = mn + rng*v;
        }
    }
    if (is.fail())
        BoxLib::Error("FABio_8bit::read() failed");

    delete [] c;
}

void
FABio_8bit::skip (std::istream& is,
                  FArrayBox&    f) const
{
    const Box& bx = f.box();
    long siz      = bx.numPts();
    Real mn, mx;
    for (int nbytes, k = 0; k < f.nComp(); k++)
    {
        is >> mn >> mx >> nbytes;
        BL_ASSERT(nbytes == siz);
        while (is.get() != '\n')
            ;
        is.seekg(siz, std::ios::cur);
    }

    if (is.fail())
        BoxLib::Error("FABio_8bit::skip() failed");
}

void
FABio_8bit::skip (std::istream& is,
                  FArrayBox&    f,
		  int           nCompToSkip) const
{
    const Box& bx = f.box();
    long siz      = bx.numPts();
    Real mn, mx;
    for (int nbytes, k = 0; k < nCompToSkip; k++)
    {
        is >> mn >> mx >> nbytes;
        BL_ASSERT(nbytes == siz);
        while (is.get() != '\n')
            ;
        is.seekg(siz, std::ios::cur);
    }

    if (is.fail())
        BoxLib::Error("FABio_8bit::skip() failed");
}

void
FABio_8bit::write_header (std::ostream&    os,
                          const FArrayBox& f,
                          int              nvar) const
{
    os << "FAB: " << FABio::FAB_8BIT << ' ' << 0 << ' ' << sys_name << '\n';
    FABio::write_header(os, f, nvar);
}

FABio_binary::FABio_binary (RealDescriptor* rd_)
    :
    rd(rd_)
{}

void
FABio_binary::write_header (std::ostream&    os,
                            const FArrayBox& f,
                            int              nvar) const
{
    os << "FAB " << *rd;
    FABio::write_header(os, f, nvar);
}

void
FABio_binary::read (std::istream& is,
                    FArrayBox&    f) const
{
    const long base_siz = f.box().numPts();
    Real* comp_ptr      = f.dataPtr(0);
    const long siz      = base_siz*f.nComp();
    RealDescriptor::convertToNativeFormat(comp_ptr, siz, is, *rd);
    if (is.fail())
        BoxLib::Error("FABio_binary::read() failed");
}

void
FABio_binary::write (std::ostream&    os,
                     const FArrayBox& f,
                     int              comp,
                     int              num_comp) const
{
    BL_ASSERT(comp >= 0 && num_comp >= 1 && (comp+num_comp) <= f.nComp());

    const long base_siz  = f.box().numPts();
    const Real* comp_ptr = f.dataPtr(comp);
    const long siz       = base_siz*num_comp;
    RealDescriptor::convertFromNativeFormat(os, siz, comp_ptr, *rd);

    if (os.fail())
        BoxLib::Error("FABio_binary::write() failed");
}

void
FABio_binary::skip (std::istream& is,
                    FArrayBox&    f) const
{
    const Box& bx = f.box();
    long base_siz = bx.numPts();
    long siz      = base_siz * f.nComp();
    is.seekg(siz*rd->numBytes(), std::ios::cur);
    if (is.fail())
        BoxLib::Error("FABio_binary::skip() failed");
}

void
FABio_binary::skip (std::istream& is,
                    FArrayBox&    f,
		    int           nCompToSkip) const
{
    const Box& bx = f.box();
    long base_siz = bx.numPts();
    long siz      = base_siz * nCompToSkip;
    is.seekg(siz*rd->numBytes(), std::ios::cur);
    if (is.fail())
        BoxLib::Error("FABio_binary::skip(..., int nCompToSkip) failed");
}

std::ostream&
operator<< (std::ostream&    os,
            const FArrayBox& f)
{
    static FABio_ascii fabio_ascii;
    fabio_ascii.write(os,f,0,f.nComp());
    return os;
}

std::istream&
operator>> (std::istream& is,
            FArrayBox&    f)
{
    FABio* fabrd = FABio::read_header(is,f);
    fabrd->read(is,f);
    delete fabrd;
    return is;
}

