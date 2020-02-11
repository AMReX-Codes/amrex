#include <AMReX_RealVect.H>

namespace amrex
{

  const RealVect RealVect::Zero = RealVect::TheZeroVector();
  const RealVect RealVect::Unit = RealVect::TheUnitVector();

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
