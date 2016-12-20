#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SPACE.H"
#include "MayDay.H"
#include "RealVect.H"
#include "Misc.H"
#include "IndexTM.H"

#include "NamespaceHeader.H"
#define CHOFFSET(object, member) (int)((char*)&(object.member) - (char*)&object)
using std::ostream;
using std::istream;
using std::ws;

RealVect tm;
size_t RealVect::io_offset = CHOFFSET(tm, vect);

const RealVect RealVect::Unit(D_DECL6(1.0,1.0,1.0,
                                      1.0,1.0,1.0));
const RealVect RealVect::Zero(D_DECL6(0.0,0.0,0.0,
                                      0.0,0.0,0.0));

const Real*
RealVect::dataPtr() const
{
  return vect;
}

Real*
RealVect::dataPtr()
{
  return vect;
}

RealVect::RealVect (D_DECL6(Real i, Real j, Real k,
                            Real l, Real m, Real n))
{
  D_EXPR6(vect[0] = i, vect[1] = j, vect[2] = k,
          vect[3] = l, vect[4] = m, vect[5] = n);
}

RealVect::RealVect (const Vector<Real>& vr )
{
  D_EXPR6(vect[0]=vr[0], vect[1]=vr[1], vect[2] = vr[2],
          vect[3]=vr[3], vect[4]=vr[4], vect[5] = vr[5]);
}

RealVect::RealVect ()
{
  D_EXPR6(vect[0]=0.0, vect[1]=0.0, vect[2] = 0.0,
          vect[3]=0.0, vect[4]=0.0, vect[5] = 0.0);
}

RealVect::RealVect(const IndexTM<Real, CH_SPACEDIM>& a_tm)
{
  D_EXPR6(vect[0] = a_tm[0], vect[1] = a_tm[1], vect[2] = a_tm[2],
          vect[3] = a_tm[3], vect[4] = a_tm[4], vect[5] = a_tm[5]);
}

RealVect&
RealVect::operator= (const RealVect &iv)
{
  D_EXPR6(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2],
          vect[3]=iv.vect[3], vect[4]=iv.vect[4], vect[5]=iv.vect[5]);
  return *this;
}

Real RealVect::dotProduct(const RealVect& a_rhs) const
{
  return D_TERM6(vect[0]*a_rhs.vect[0], +
                 vect[1]*a_rhs.vect[1], +
                 vect[2]*a_rhs.vect[2], +
                 vect[3]*a_rhs.vect[3], +
                 vect[4]*a_rhs.vect[4], +
                 vect[5]*a_rhs.vect[5]);
}

bool
RealVect::operator== (const RealVect& p) const
{
  return D_TERM6(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2], &&
                 vect[3] == p[3], && vect[4] == p[4], && vect[5] == p[5]);
}

bool
RealVect::operator!= (const RealVect& p) const
{
  return D_TERM6(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2], ||
                 vect[3] != p[3], || vect[4] != p[4], || vect[5] != p[5]);
}

RealVect&
RealVect::operator+= (Real s)
{
  D_EXPR6(vect[0] += s, vect[1] += s, vect[2] += s,
          vect[3] += s, vect[4] += s, vect[5] += s);
  return *this;
}

RealVect&
RealVect::operator+= (const RealVect& p)
{
  D_EXPR6(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2],
          vect[3] += p[3], vect[4] += p[4], vect[5] += p[5]);
  return *this;
}

RealVect&
RealVect::operator*= (Real s)
{
  D_EXPR6(vect[0] *= s, vect[1] *= s, vect[2] *= s,
          vect[3] *= s, vect[4] *= s, vect[5] *= s);
  return *this;
}

RealVect&
RealVect::operator*= (const RealVect &p)
{
  D_EXPR6(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2],
          vect[3] *= p[3], vect[4] *= p[4], vect[5] *= p[5]);
  return *this;
}

//XXXRealVect
//XXXRealVect::operator* (const RealVect& p) const
//XXX{
//XXX  RealVect v(D_DECL6(vect[0]*p[0], vect[1]*p[1], vect[2]*p[2]));
//XXX  return v;
//XXX}

RealVect
RealVect::operator* (Real s) const
{
  RealVect v(D_DECL6(vect[0]*s, vect[1]*s, vect[2]*s,
                     vect[3]*s, vect[4]*s, vect[5]*s));
  return v;
}

RealVect
RealVect::operator- (Real s) const
{
  RealVect v(D_DECL6(vect[0]-s, vect[1]-s, vect[2]-s,
                     vect[3]-s, vect[4]-s, vect[5]-s));
  return v;
}

RealVect
RealVect::operator+ (Real s) const
{
  RealVect v(D_DECL6(vect[0]+s, vect[1]+s, vect[2]+s,
                     vect[3]+s, vect[4]+s, vect[5]+s));
  return v;
}

RealVect&
RealVect::operator/= (Real s)
{
  D_EXPR6(vect[0] /= s, vect[1] /= s, vect[2] /= s,
          vect[3] /= s, vect[4] /= s, vect[5] /= s);
  return *this;
}

RealVect&
RealVect::operator/= (const RealVect& p)
{
  D_EXPR6(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2],
          vect[3] /= p[3], vect[4] /= p[4], vect[5] /= p[5]);
  return *this;
}

//XXXRealVect
//XXXRealVect::operator/ (const RealVect & p) const
//XXX{
//XXX  RealVect result( D_DECL6( vect[0] / p[0], vect[1] / p[1], vect[2] / p[2] ) );
//XXX  return result ;
//XXX}

RealVect
RealVect::operator/ (Real s) const
{
  RealVect result( D_DECL6( vect[0] / s, vect[1] / s, vect[2] / s,
                            vect[3] / s, vect[4] / s, vect[5] / s) );
  return result ;
}

int
RealVect::minDir(const bool& a_doAbs) const
{
  int mDir = 0;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      if (a_doAbs)
        {
          if (Abs(vect[idir]) < Abs(vect[mDir]))
            {
              mDir = idir;
            }
        }
      else
        {
          if (vect[idir] < vect[mDir])
            {
              mDir = idir;
            }
        }
    }
  return mDir;
}

int
RealVect::maxDir(const bool& a_doAbs) const
{
  int mDir = 0;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      if (a_doAbs)
        {
          if (Abs(vect[idir]) > Abs(vect[mDir]))
            {
              mDir = idir;
            }
        }
      else
        {
          if (vect[idir] > vect[mDir])
            {
              mDir = idir;
            }
        }
    }
  return mDir;
}

RealVect
BASISREALV (int dir)
{
  CH_assert(dir >= 0 && dir < SpaceDim);
  RealVect tmp = RealVect::Zero ;
  tmp.vect[dir] = 1;
  return tmp;
}

RealVect
operator/ (Real            s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s/p[0], s/p[1], s/p[2],
                          s/p[3], s/p[4], s/p[5]));
}
RealVect
operator+ (Real            s,
           const RealVect& p)
{
  return RealVect(D_DECL6(p[0] + s, p[1] + s, p[2] + s,
                          p[3] + s, p[4] + s, p[5] + s));
}

RealVect
operator- (Real            s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s - p[0], s - p[1], s - p[2],
                          s - p[3], s - p[4], s - p[5]));
}

RealVect
operator* (Real            s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s * p[0], s * p[1], s * p[2],
                          s * p[3], s * p[4], s * p[5]));
}

RealVect
operator/ (const RealVect& s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s[0] / p[0], s[1] /p[1], s[2] / p[2],
                          s[3] / p[3], s[4] /p[4], s[5] / p[5]));
}

RealVect
operator+ (const RealVect& s,
           const RealVect& p)
{
  return RealVect(D_DECL6(p[0] + s[0], p[1] +s[1], p[2] + s[2],
                          p[3] + s[3], p[4] +s[4], p[5] + s[5]));
}

RealVect
operator- (const RealVect& s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s[0] - p[0], s[1] - p[1], s[2] - p[2],
                          s[3] - p[3], s[4] - p[4], s[5] - p[5]));
}

RealVect
operator* (const RealVect& s,
           const RealVect& p)
{
  return RealVect(D_DECL6(p[0] * s[0], p[1] *s[1], p[2] * s[2],
                          p[3] * s[3], p[4] *s[4], p[5] * s[5]));
}

std::ostream&
operator<< (std::ostream& ostr, const RealVect& p)
{
  ostr << "(" << D_TERM6( p[0] ,<< "," << p[1], << "," << p[2], << "," <<
                          p[3] ,<< "," << p[4], << "," << p[5] ) << ")" ;
  if (ostr.fail()) MayDay::Error("operator<<(ostream&,RealVect&) failed");
  return ostr;
}

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

template < >
int linearSize(const RealVect& vindex)
{
  return sizeof(RealVect);
}

//RealVect specialization of linearIn
template < >
void linearIn(RealVect& a_outputT, const void* const inBuf)
{
  unsigned char* bob = (unsigned char*)inBuf;
  unsigned char* to = (unsigned char*)&a_outputT;
  memcpy(to + RealVect::io_offset, bob, SpaceDim*sizeof(Real));
}

//RealVect specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const RealVect& a_inputT)
{
  unsigned char* bob = (unsigned char*)a_outBuf;
  const unsigned char* from = (const unsigned char*)&a_inputT;
  memcpy(bob, from + RealVect::io_offset, SpaceDim*sizeof(Real));
}

#include "BaseNamespaceFooter.H"
