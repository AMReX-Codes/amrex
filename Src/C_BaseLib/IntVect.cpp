//
// $Id: IntVect.cpp,v 1.8 2001-07-17 23:02:24 lijewski Exp $
//

#include <BLassert.H>
#include <BoxLib.H>
#include <IntVect.H>
#include <IndexType.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

const IntVect&
IntVect::TheUnitVector ()
{
    static const IntVect Unit(D_DECL(1,1,1));
    return Unit;
}

const IntVect&
IntVect::TheZeroVector ()
{
    static const IntVect Zero(D_DECL(0,0,0));
    return Zero;
}

const IntVect&
IntVect::TheNodeVector ()
{
    static const IntVect Node(D_DECL(IndexType::NODE,IndexType::NODE,IndexType::NODE));
    return Node;
}

const IntVect&
IntVect::TheCellVector ()
{
    static const IntVect Cell(D_DECL(IndexType::CELL,IndexType::CELL,IndexType::CELL));
    return Cell;
}

IntVect::IntVect ()
{
    D_EXPR(vect[0] = 0, vect[1] = 0, vect[2] = 0);
}

IntVect::IntVect (D_DECL(int i, int j, int k))
{
    D_EXPR(vect[0] = i, vect[1] = j, vect[2] = k);
}

IntVect::IntVect (const int *a)
{
    D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

IntVect::IntVect (const Array<int> &a)
{
    BL_ASSERT(a.length() == BL_SPACEDIM);
    D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

IntVect::IntVect (const IntVect &iv)
{
    D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
}

IntVect&
IntVect::operator= (const IntVect &iv)
{
    D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
    return *this;
}

int&
IntVect::operator[] (int i)
{
    BL_ASSERT(i>=0 && i < BL_SPACEDIM);
    return vect[i];
}

int
IntVect::operator[] (int i) const
{
    BL_ASSERT(i>=0 && i < BL_SPACEDIM);
    return vect[i];
}

IntVect&
IntVect::setVal (int i,
                 int val)
{
    BL_ASSERT(i >=0 && i < BL_SPACEDIM);
    vect[i] = val;
    return *this;
}

const int*
IntVect::getVect () const
{
    return vect;
}

bool
IntVect::operator== (const IntVect& p) const
{
    return D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
}

bool
IntVect::operator!= (const IntVect& p) const
{
    return D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
}

bool
IntVect::operator< (const IntVect& p) const
{
    return D_TERM(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2]);
}

bool
IntVect::operator<= (const IntVect& p) const
{
    return D_TERM(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2]);
}

bool
IntVect::operator> (const IntVect& p) const
{
    return D_TERM(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2]);
}

bool
IntVect::operator>= (const IntVect& p) const
{
    return D_TERM(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2]);
}

bool
IntVect::lexLT (const IntVect &s) const
{
#define LLT0 (vect[0] < s[0])
#define LLT1 ((vect[0] == s[0]) && (vect[1] < s[1]))
#define LLT2 ((vect[1] == s[1]) && (vect[2] < s[2]))
#if   BL_SPACEDIM == 1
    return LLT0;
#elif BL_SPACEDIM == 2
    return LLT0 || LLT1;
#elif BL_SPACEDIM == 3
    return LLT0 || (vect[0]==s[0] && ((vect[1] < s[1] || LLT2)));
#endif
#undef LLT0
#undef LLT1
#undef LLT2
}

bool
IntVect::lexGT (const IntVect& s) const
{
#define LGT0 (vect[0] > s[0])
#define LGT1 ((vect[0] == s[0]) && (vect[1] > s[1]))
#define LGT2 ((vect[1] == s[1]) && (vect[2] > s[2]))
#if   BL_SPACEDIM == 1
   return LGT0;
#elif BL_SPACEDIM == 2
   return LGT0 || LGT1;
#elif BL_SPACEDIM == 3
   return LGT0 || (vect[0] == s[0] && ((vect[1] > s[1] || LGT2)));
#endif
#undef LGT0
#undef LGT1
#undef LGT2
}

IntVect
IntVect::operator+ () const
{
    return *this;
}

IntVect
IntVect::operator- () const
{
    return IntVect(D_DECL(-vect[0], -vect[1], -vect[2] ));
}

IntVect&
IntVect::operator+= (int s)
{
    D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
}

IntVect&
IntVect::operator+= (const IntVect& p)
{
    D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]);
    return *this;
}

IntVect&
IntVect::operator*= (int s)
{
    D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
    return *this;
}

IntVect&
IntVect::operator*= (const IntVect &p)
{
    D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]);
    return *this;
}

IntVect&
IntVect::operator/= (int s)
{
    D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s);
    return *this;
}

IntVect&
IntVect::operator/= (const IntVect& p)
{
    D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]);
    return *this;
}

IntVect&
IntVect::operator-= (int s)
{
    D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s);
    return *this;
}

IntVect&
IntVect::operator-= (const IntVect& p)
{
    D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]);
    return *this;
}

IntVect
IntVect::operator+ (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2]));
}

IntVect
IntVect::operator+ (int s) const
{
    return IntVect(D_DECL(vect[0] + s, vect[1] + s, vect[2] + s));
}

IntVect
IntVect::operator- (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2]));
}

IntVect
IntVect::operator- (int s) const
{
    return IntVect(D_DECL(vect[0] - s, vect[1] - s, vect[2] - s));
}

IntVect
IntVect::operator* (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2]));
}

IntVect
IntVect::operator* (int s) const
{
    return IntVect(D_DECL(vect[0] * s, vect[1] * s, vect[2] * s));
}

IntVect
IntVect::operator/ (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2]));
}

IntVect
IntVect::operator/ (int s) const
{
    return IntVect(D_DECL(vect[0] / s, vect[1] / s, vect[2] / s));
}

IntVect&
IntVect::min (const IntVect& p)
{
    D_EXPR(vect[0] = std::min(vect[0], p.vect[0]),
           vect[1] = std::min(vect[1], p.vect[1]),
           vect[2] = std::min(vect[2], p.vect[2]));
    return *this;
}

IntVect&
IntVect::max (const IntVect& p)
{
    D_EXPR(vect[0] = std::max(vect[0], p.vect[0]),
           vect[1] = std::max(vect[1], p.vect[1]),
           vect[2] = std::max(vect[2], p.vect[2]));
    return *this;
}

IntVect&
IntVect::scale (int s)
{
    D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
    return *this;
}

IntVect&
IntVect::reflect (int ref_ix,
                  int idir)
{
    BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
    vect[idir] = -vect[idir] + 2*ref_ix;
    return *this;
}

IntVect&
IntVect::shift (int coord,
                int s)
{
    BL_ASSERT(coord >= 0 && coord < BL_SPACEDIM);
    vect[coord] += s;
    return *this;
}

IntVect&
IntVect::shift (const IntVect& iv)
{
    *this += iv;
    return *this;
}

IntVect&
IntVect::diagShift (int s)
{
    D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
}

IntVect
operator+ (int            s,
           const IntVect& p)
{
    return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

IntVect
operator- (int            s,
           const IntVect& p)
{
    return IntVect(D_DECL(s - p[0], s - p[1], s - p[2]));
}

IntVect
operator* (int            s,
           const IntVect& p)
{
    return IntVect(D_DECL(s * p[0], s * p[1], s * p[2]));
}

IntVect
scale (const IntVect& p,
       int            s)
{
    return IntVect(D_DECL(s * p[0], s * p[1], s * p[2]));
}

IntVect
diagShift (const IntVect &p, int s)
{
    return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

IntVect
min (const IntVect& p1,
     const IntVect& p2)
{
    IntVect p(p1);
    return p.min(p2);
}

IntVect
max (const IntVect& p1,
     const IntVect& p2)
{
    IntVect p(p1);
    return p.max(p2);
}

IntVect
BASISV (int dir)
{
    BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM);
    IntVect tmp;
    tmp.vect[dir] = 1;
    return tmp;
}

IntVect
reflect (const IntVect& a,
         int            ref_ix,
         int            idir)
{
    BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
    IntVect b(a);
    b.vect[idir] = -b.vect[idir] + 2*ref_ix;
    return b;
}

IntVect
coarsen (const IntVect& p,
         int            s)
{
    BL_ASSERT(s > 0);
    return IntVect(
        D_DECL((p.vect[0]<0) ? -abs(p.vect[0]+1)/s-1 : p.vect[0]/s ,
               (p.vect[1]<0) ? -abs(p.vect[1]+1)/s-1 : p.vect[1]/s ,
               (p.vect[2]<0) ? -abs(p.vect[2]+1)/s-1 : p.vect[2]/s ));
}

IntVect
coarsen (const IntVect& p1,
         const IntVect& p2)
{
    BL_ASSERT(p2 > IntVect::TheZeroVector());
    return IntVect(
        D_DECL(
        (p1.vect[0]<0)?-abs(p1.vect[0]+1)/p2.vect[0]-1:p1.vect[0]/p2.vect[0],
        (p1.vect[1]<0)?-abs(p1.vect[1]+1)/p2.vect[1]-1:p1.vect[1]/p2.vect[1],
        (p1.vect[2]<0)?-abs(p1.vect[2]+1)/p2.vect[2]-1:p1.vect[2]/p2.vect[2])
    );
}

IntVect&
IntVect::coarsen (int s)
{
    BL_ASSERT(s > 0);
    for (int i = 0; i < BL_SPACEDIM; ++i)
        vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
    return *this;
}

IntVect&
IntVect::coarsen (const IntVect& p)
{
    BL_ASSERT(p > IntVect::TheZeroVector());
    for (int i = 0; i <BL_SPACEDIM; ++i)
    {
        const int s = p.vect[i];
        vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
    }
    return *this;
}

//
// Returns IntVect which is the componentwise integer projection
// of IntVect p1 by IntVect p2.
//

std::ostream&
operator<< (std::ostream&  os,
            const IntVect& p)
{
    os << D_TERM( '(' << p[0] , <<
                  ',' << p[1] , <<
                  ',' << p[2])  << ')';
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            IntVect&      p)
{
    is >> std::ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> p[0],
               is.ignore(BL_IGNORE_MAX, ',') >> p[1],
               is.ignore(BL_IGNORE_MAX, ',') >> p[2]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
        D_EXPR(is.ignore(BL_IGNORE_MAX, '<') >> p[0],
               is.ignore(BL_IGNORE_MAX, ',') >> p[1],
               is.ignore(BL_IGNORE_MAX, ',') >> p[2]);
        is.ignore(BL_IGNORE_MAX, '>');
    }
    else
        BoxLib::Error("operator>>(istream&,IntVect&): expected \'(\' or \'<\'");

    if (is.fail())
        BoxLib::Error("operator>>(istream&,IntVect&) failed");

    return is;
}

void
IntVect::printOn (std::ostream& os) const
{
    os << "IntVect: " << *this << '\n';
}

void
IntVect::dumpOn (std::ostream& os) const
{
    os << "IntVect(" << BoxLib::version << ")= " << *this << '\n';
}

#ifdef BL_NAMESPACE
}
#endif

