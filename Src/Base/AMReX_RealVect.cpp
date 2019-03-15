#include "AMReX_SPACE.H"
#include "AMReX_RealVect.H"
#include "AMReX_Utility.H"
using std::ostream;
using std::istream;
using std::ws;

namespace amrex
{

  const RealVect RealVect::Unit(AMREX_D_DECL(1.0,1.0,1.0));

  const RealVect RealVect::Zero(AMREX_D_DECL(0.0,0.0,0.0));

  const Real*
  RealVect::dataPtr() const noexcept
  {
    return vect;
  }

  Real*
  RealVect::dataPtr() noexcept
  {
    return vect;
  }

  RealVect::RealVect (AMREX_D_DECL(Real i, Real j, Real k)) noexcept
  {
    AMREX_D_EXPR(vect[0] = i, vect[1] = j, vect[2] = k);
  }

  RealVect::RealVect (const std::vector<Real>& vr ) noexcept
  {
    AMREX_D_EXPR(vect[0]=vr[0], vect[1]=vr[1], vect[2] = vr[2]);
  }

  RealVect::RealVect () noexcept
  {
    AMREX_D_EXPR(vect[0]=0.0, vect[1]=0.0, vect[2] = 0.0);
  }


  RealVect&
  RealVect::operator= (const RealVect &iv) noexcept
  {
    AMREX_D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
    return *this;
  }

  Real RealVect::dotProduct(const RealVect& a_rhs) const noexcept
  {
    return AMREX_D_TERM(vect[0]*a_rhs.vect[0], +
                  vect[1]*a_rhs.vect[1], +
                  vect[2]*a_rhs.vect[2]);
  }

  bool
  RealVect::operator== (const RealVect& p) const noexcept
  {
    return AMREX_D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
  }

  bool
  RealVect::operator!= (const RealVect& p) const noexcept
  {
    return AMREX_D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
  }

  RealVect&
  RealVect::operator+= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
  }

  RealVect&
  RealVect::operator+= (const RealVect& p) noexcept
  {
    AMREX_D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]);
    return *this;
  }

  RealVect&
  RealVect::operator*= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
    return *this;
  }

  RealVect&
  RealVect::operator*= (const RealVect &p) noexcept
  {
    AMREX_D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]);
    return *this;
  }

  RealVect
  RealVect::operator* (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]*s, vect[1]*s, vect[2]*s));
    return v;
  }

  RealVect
  RealVect::operator- (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]-s, vect[1]-s, vect[2]-s));
    return v;
  }

  RealVect
  RealVect::operator+ (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]+s, vect[1]+s, vect[2]+s));
    return v;
  }

  RealVect&
  RealVect::operator/= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s);
    return *this;
  }

  RealVect&
  RealVect::operator/= (const RealVect& p) noexcept
  {
    AMREX_D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]);
    return *this;
  }

  RealVect
  RealVect::operator/ (Real s) const noexcept
  {
    RealVect result( AMREX_D_DECL( vect[0] / s, vect[1] / s, vect[2] / s));
    return result ;
  }

  int
  RealVect::minDir(const bool& a_doAbs) const noexcept
  {
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
      {
        if (a_doAbs)
          {
            if (std::abs(vect[idir]) < std::abs(vect[mDir]))
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
  RealVect::maxDir(const bool& a_doAbs) const noexcept
  {
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
      {
        if (a_doAbs)
          {
            if (std::abs(vect[idir]) > std::abs(vect[mDir]))
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
  BASISREALV (int dir) noexcept
  {
    assert(dir >= 0 && dir < SpaceDim);
    RealVect tmp = RealVect::Zero ;
    tmp.vect[dir] = 1;
    return tmp;
  }

  RealVect
  operator/ (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s/p[0], s/p[1], s/p[2]));
  }
  RealVect
  operator+ (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
  }

  RealVect
  operator- (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s - p[0], s - p[1], s - p[2]));
  }

  RealVect
  operator* (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
  }

  RealVect
  operator/ (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s[0] / p[0], s[1] /p[1], s[2] / p[2]));
  }

  RealVect
  operator+ (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(p[0] + s[0], p[1] +s[1], p[2] + s[2]));
  }

  RealVect
  operator- (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s[0] - p[0], s[1] - p[1], s[2] - p[2]));
  }

  RealVect
  operator* (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(p[0] * s[0], p[1] *s[1], p[2] * s[2]));
  }

  std::ostream&
  operator<< (std::ostream& ostr, const RealVect& p)
  {
    ostr << "(" << AMREX_D_TERM ( p[0] ,<< "," << p[1], << "," << p[2]) << ")" ;
    return ostr;
  }

  std::istream&
  operator>> (std::istream& is,
              RealVect&      iv)
  {
    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        AMREX_D_EXPR(is >> iv[0],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[1],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[2]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,IntVect&): expected \'(\'");
    }

    if (is.fail())
        amrex::Error("operator>>(istream&,IntVect&) failed");

    return is;
}
} //namespace amrex
