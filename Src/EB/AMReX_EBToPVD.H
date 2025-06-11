#ifndef AMReX_EBToPVD_H_
#define AMReX_EBToPVD_H_
#include <AMReX_Config.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>

#include <vector>
#include <array>
#include <iosfwd>

namespace amrex {

class EBToPVD {

public:
   EBToPVD() = default;

   void EBToPolygon(const Real* problo, const Real* dx,
         const Box & bx, Array4<EBCellFlag const> const& flag,
         Array4<Real const> const& bcent,
         Array4<Real const> const& apx, Array4<Real const> const& apy, Array4<Real const> const& apz);

   void WriteEBVTP(int myID) const;
   static void WritePVTP(int nProcs) ;

   void EBGridCoverage(int myID, const Real* problo, const Real* dx,
         const Box &bx, Array4<EBCellFlag const> const& flag);

private:
   void reorder_polygon(const std::vector<std::array<Real,3>>& lpoints, std::array<int,7>& lconnect,
         const std::array<Real,3>& lnormal);

   // Calculates the Hesse Normal FOrm corresponding to normal and centroid
   static void calc_hesse(Real& distance, std::array<Real,3>& n0, Real& p,
         const std::array<Real,3>& normal, const std::array<Real,3>& centroid) ;

   // Fills the alpha vector
   static void calc_alpha(std::array<Real,12>& alpha,
         const std::array<Real,3>& n0, Real p,
         const std::array<std::array<Real,3>,8>& vertex,
         const Real* dx) ;

   // Fills count and flags selecting the alphas which are in (0,1)
   static void calc_intersects(int& int_count, std::array<bool,12>& intersects_flags,
         const std::array<Real,12>& alpha) ;

   void print_points(std::ofstream& myfile) const;
   void print_connectivity(std::ofstream& myfile) const;
   void print_grids(std::ofstream& myfile) const;

   std::vector<std::array<Real,3>> m_points;
   std::vector<std::array<int,7>> m_connectivity;
   int m_grid{0};

};

}

#endif
