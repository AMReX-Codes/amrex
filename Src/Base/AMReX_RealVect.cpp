#include "AMReX_SPACE.H"
#include "AMReX_RealVect.H"
#include "AMReX_Utility.H"

using std::ostream;
using std::istream;
using std::ws;

namespace amrex
{

  const RealVect RealVect::Zero = RealVect::TheZeroVector();
  const RealVect RealVect::Unit = RealVect::TheUnitVector();

  AMREX_GPU_HOST_DEVICE
  const Real*
  RealVect::dataPtr() const noexcept
  {
    return vect;
  }

  AMREX_GPU_HOST_DEVICE
  Real*
  RealVect::dataPtr() noexcept
  {
    return vect;
  }

  AMREX_GPU_HOST_DEVICE
  Real RealVect::dotProduct(const RealVect& a_rhs) const noexcept
  {
    return AMREX_D_TERM(vect[0]*a_rhs.vect[0], +
                  vect[1]*a_rhs.vect[1], +
                  vect[2]*a_rhs.vect[2]);
  }

#if (AMREX_SPACEDIM == 3)
  AMREX_GPU_HOST_DEVICE
  RealVect RealVect::crossProduct(const RealVect& a_rhs) const noexcept
  {
    RealVect tmp(vect[1]*a_rhs[2] - vect[2]*a_rhs[1],
                 vect[2]*a_rhs[0] - vect[0]*a_rhs[2],
                 vect[0]*a_rhs[1] - vect[1]*a_rhs[0]);
    return tmp;
  }
#endif

  AMREX_GPU_HOST_DEVICE
  bool
  RealVect::operator== (const RealVect& p) const noexcept
  {
    return AMREX_D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
  }

  AMREX_GPU_HOST_DEVICE
  bool
  RealVect::operator!= (const RealVect& p) const noexcept
  {
    return AMREX_D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator+= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator+= (const RealVect& p) noexcept
  {
    AMREX_D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator*= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator*= (const RealVect &p) noexcept
  {
    AMREX_D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  RealVect::operator* (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]*s, vect[1]*s, vect[2]*s));
    return v;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  RealVect::operator- (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]-s, vect[1]-s, vect[2]-s));
    return v;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  RealVect::operator+ (Real s) const noexcept
  {
    RealVect v(AMREX_D_DECL(vect[0]+s, vect[1]+s, vect[2]+s));
    return v;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator/= (Real s) noexcept
  {
    AMREX_D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect&
  RealVect::operator/= (const RealVect& p) noexcept
  {
    AMREX_D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]);
    return *this;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  RealVect::operator/ (Real s) const noexcept
  {
    RealVect result( AMREX_D_DECL( vect[0] / s, vect[1] / s, vect[2] / s));
    return result ;
  }

  AMREX_GPU_HOST_DEVICE
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

  AMREX_GPU_HOST_DEVICE
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

  AMREX_GPU_HOST_DEVICE
  RealVect
  BASISREALV (int dir) noexcept
  {
    AMREX_ASSERT(dir >= 0 && dir < SpaceDim);
    RealVect tmp(AMREX_D_DECL(0.,0.,0.));
    tmp.vect[dir] = 1;
    return tmp;
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator/ (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s/p[0], s/p[1], s/p[2]));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator+ (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator- (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s - p[0], s - p[1], s - p[2]));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator* (Real            s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator/ (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s[0] / p[0], s[1] /p[1], s[2] / p[2]));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator+ (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(p[0] + s[0], p[1] +s[1], p[2] + s[2]));
  }

  AMREX_GPU_HOST_DEVICE
  RealVect
  operator- (const RealVect& s,
             const RealVect& p) noexcept
  {
    return RealVect(AMREX_D_DECL(s[0] - p[0], s[1] - p[1], s[2] - p[2]));
  }

  AMREX_GPU_HOST_DEVICE
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
