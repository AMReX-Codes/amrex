#include <AMReX_GeomIntersectUtils.H>

namespace amrex
{
  Edge::Edge(const IntVect& lhs, const IntVect& rhs)
  {
    ID = -1;
    if (lhs < rhs) {
      IV_l = lhs;
      IV_r = rhs;
    } else {
      IV_l = rhs;
      IV_r = lhs;
    }
  }

  bool
  Edge::operator== (const Edge& rhs) const
  {
    return IV_l==rhs.IV_l && IV_r==rhs.IV_r;
  }

  bool
  Edge::operator< (const Edge& rhs) const
  {
    if (IV_l==rhs.IV_l) {
      return IV_r < rhs.IV_r;
    }
    return IV_l < rhs.IV_l;
  }
  
  std::ostream& operator<<(std::ostream& os, const Edge& e)
  {
    os << e.IV_l << " " << e.IV_r << " " << e.ID;
    return os;
  }
}
