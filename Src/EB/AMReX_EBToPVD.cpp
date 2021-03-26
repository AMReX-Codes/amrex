#include <AMReX_EBToPVD.H>
#include <AMReX_BLassert.H>
#include <AMReX_Dim3.H>

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

namespace {
amrex::Real dot_product(const std::array<amrex::Real,3>& a, const std::array<amrex::Real,3>& b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

bool intersects(amrex::Real val)
{
   return (val > 0.0 && val < 1.0);
}

}

namespace amrex {

void EBToPVD::EBToPolygon(const Real* problo, const Real* dx,
      const Box & bx, Array4<EBCellFlag const> const& flag,
      Array4<Real const> const& bcent,
      Array4<Real const> const& apx, Array4<Real const> const& apy, Array4<Real const> const& apz)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);

   for(int k = lo.z; k <= hi.z; ++k) {
      for(int j = lo.y; j <= hi.y; ++j) {
         for(int i = lo.x; i <= hi.x; ++i) {
            // NOTE: do not skip fully enclosed cells (is_covered_cell), as this seems
            // to skip thin walls in the domain:
            // if(.not.is_regular_cell(flag(i,j,k)) .and. &
            // .not.is_covered_cell(flag(i,j,k))) then

            // Instead only look for EBs
            // if( .not.is_regular_cell(flag(i,j,k))) then

            // If covered cells are accounted for in this loop, a FPE arises
            // since apnorm is zero.

            if(flag(i,j,k).isSingleValued()) {
               Real axm = apx(i  ,j  ,k  );
               Real axp = apx(i+1,j  ,k  );
               Real aym = apy(i  ,j  ,k  );
               Real ayp = apy(i  ,j+1,k  );
               Real azm = apz(i  ,j  ,k  );
               Real azp = apz(i  ,j  ,k+1);

               Real apnorm = std::sqrt((axm-axp)*(axm-axp) + (aym-ayp)*(aym-ayp) + (azm-azp)*(azm-azp));
               Real apnorminv = 1.0/apnorm;

               std::array<Real,3> normal, centroid;
               std::array<std::array<Real,3>,8> vertex;

               normal[0] = (axp-axm) * apnorminv;
               normal[1] = (ayp-aym) * apnorminv;
               normal[2] = (azp-azm) * apnorminv;

               // convert bcent to global coordinate system centered at plo
               centroid[0] = problo[0] + bcent(i,j,k,0)*dx[0] + (i + 0.5)*dx[0];
               centroid[1] = problo[1] + bcent(i,j,k,1)*dx[1] + (j + 0.5)*dx[1];
               centroid[2] = problo[2] + bcent(i,j,k,2)*dx[2] + (k + 0.5)*dx[2];

               // vertices of bounding cell (i,j,k)
               vertex[0] = {problo[0] + (i  )*dx[0], problo[1] + (j  )*dx[1], problo[2] + (k  )*dx[2]};
               vertex[1] = {problo[0] + (i+1)*dx[0], problo[1] + (j  )*dx[1], problo[2] + (k  )*dx[2]};
               vertex[2] = {problo[0] + (i  )*dx[0], problo[1] + (j+1)*dx[1], problo[2] + (k  )*dx[2]};
               vertex[3] = {problo[0] + (i+1)*dx[0], problo[1] + (j+1)*dx[1], problo[2] + (k  )*dx[2]};
               vertex[4] = {problo[0] + (i  )*dx[0], problo[1] + (j  )*dx[1], problo[2] + (k+1)*dx[2]};
               vertex[5] = {problo[0] + (i+1)*dx[0], problo[1] + (j  )*dx[1], problo[2] + (k+1)*dx[2]};
               vertex[6] = {problo[0] + (i  )*dx[0], problo[1] + (j+1)*dx[1], problo[2] + (k+1)*dx[2]};
               vertex[7] = {problo[0] + (i+1)*dx[0], problo[1] + (j+1)*dx[1], problo[2] + (k+1)*dx[2]};

               // NOTE: this seems to be unncessary:
               // skip cells that have a tiny intersection and cells that have
               // the centroid on a face/edge/corner
               // if(apnorm > stol    .and. &
               //   vertex(1,1) < centroid(1) .and. centroid(1) < vertex(8,1) .and. &
               //   vertex(1,2) < centroid(2) .and. centroid(2) < vertex(8,2) .and. &
               //   vertex(1,3) < centroid(3) .and. centroid(3) < vertex(8,3)) then

               // Compute EB facets for current cell
               int count;
               Real distance, p;
               std::array<Real,3> n0;
               std::array<Real,12> alpha;
               std::array<bool,12> alpha_intersect;

               calc_hesse(distance, n0, p, normal, centroid);
               calc_alpha(alpha, n0, p, vertex, dx);
               calc_intersects(count, alpha_intersect, alpha);

               // If the number of facet "contained" in does not describe a facet:
               // ... I.e. there's less than 3 (not even a triangle) or more than 6
               // ... (I have no idea what that is):
               //   => Move the centroid a little back and forth along the normal
               //      to see if that makes a difference:

               if((count < 3) || (count > 6)) {
                  int count_d;
                  Real p_d;
                  std::array<Real,3> n0_d;
                  std::array<Real,12> alpha_d;
                  std::array<bool,12> alpha_d_intersect;

                  Real tol = std::min({dx[0], dx[1], dx[2]})/100; // bit of a fudge factor

                  std::array<Real,3> centroid_d;
                  for(int idim = 0; idim < 3; ++idim) {
                     centroid_d[idim] = centroid[idim] + tol*normal[idim];
                  }

                  calc_hesse(distance, n0_d, p_d, normal, centroid_d);
                  calc_alpha(alpha_d, n0_d, p_d, vertex, dx);
                  calc_intersects(count_d, alpha_d_intersect, alpha_d);
                  if((count_d >= 3) && (count_d <= 6)) {
                     count = count_d;
                     alpha_intersect = alpha_d_intersect;
                  }

                  for(int idim = 0; idim < 3; ++idim) {
                     centroid_d[idim] = centroid[idim] - tol*normal[idim];
                  }

                  calc_hesse(distance, n0_d, p_d, normal, centroid_d);
                  calc_alpha(alpha_d, n0_d, p_d, vertex, dx);
                  calc_intersects(count_d, alpha_d_intersect, alpha_d);
                  if((count_d >= 3) && (count_d <= 6)) {
                     count = count_d;
                     alpha_intersect = alpha_d_intersect;
                  }
               }
               // I know this was a bit of a hack, but it's the only way I prevent
               // missing facets...

               if((count >=3) && (count <=6)) {
                  m_connectivity.push_back({0,0,0,0,0,0,0});

                  // calculate intersection points.
                  std::array<std::array<Real,3>,12> apoints;

                  std::array<Real,3> ihat = {1, 0, 0};
                  std::array<Real,3> jhat = {0, 1, 0};
                  std::array<Real,3> khat = {0, 0, 1};

                  for(int idim = 0; idim < 3; ++idim) {
                     apoints[ 0][idim] = vertex[0][idim] + ihat[idim]*dx[0]*alpha[ 0];
                     apoints[ 1][idim] = vertex[1][idim] + jhat[idim]*dx[1]*alpha[ 1];
                     apoints[ 2][idim] = vertex[2][idim] + ihat[idim]*dx[0]*alpha[ 2];
                     apoints[ 3][idim] = vertex[0][idim] + jhat[idim]*dx[1]*alpha[ 3];
                     apoints[ 4][idim] = vertex[0][idim] + khat[idim]*dx[2]*alpha[ 4];
                     apoints[ 5][idim] = vertex[1][idim] + khat[idim]*dx[2]*alpha[ 5];
                     apoints[ 6][idim] = vertex[3][idim] + khat[idim]*dx[2]*alpha[ 6];
                     apoints[ 7][idim] = vertex[2][idim] + khat[idim]*dx[2]*alpha[ 7];
                     apoints[ 8][idim] = vertex[4][idim] + ihat[idim]*dx[0]*alpha[ 8];
                     apoints[ 9][idim] = vertex[5][idim] + jhat[idim]*dx[1]*alpha[ 9];
                     apoints[10][idim] = vertex[6][idim] + ihat[idim]*dx[0]*alpha[10];
                     apoints[11][idim] = vertex[4][idim] + jhat[idim]*dx[1]*alpha[11];
                  }

                  // store intersections with grid cell alpha in [0,1]
                  for(int lc1 = 0; lc1 < 12; ++lc1) {
                     if(alpha_intersect[lc1]) {
                        m_points.push_back(apoints[lc1]);
                        int lc2 = m_connectivity.back()[0]+1;
                        m_connectivity.back()[0] = lc2;
                        m_connectivity.back()[lc2] = m_points.size()-1;
                     }
                  }

                  reorder_polygon(m_points, m_connectivity.back(), n0);
               }
            }
         }
      }
   };
}

void EBToPVD::WriteEBVTP(const int myID) const
{
   std::stringstream ss;
   ss << std::setw(8) << std::setfill('0') << myID;
   std::string cID = "eb_" + ss.str() + ".vtp";

   std::ofstream myfile(cID);
   if(myfile.is_open()) {
      myfile.precision(6);
      myfile << "<?xml version=\"1.0\"?>\n";
      myfile << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      myfile << "<PolyData>\n";
      myfile << "<Piece NumberOfPoints=\"" << m_points.size() << "\" NumberOfVerts=\"0\" "
         << "NumberOfLines=\"0\" NumberOfString=\"0\" NumberOfPolys=\" "
         << m_connectivity.size() << "\">\n";
      print_points(myfile);
      print_connectivity(myfile);
      myfile << "<PointData></PointData>\n";
      myfile << "<CellData></CellData>\n";
      myfile << "</Piece>\n";
      myfile << "</PolyData>\n";
      myfile << "</VTKFile>\n";

      myfile.close();
   }
}

void EBToPVD::WritePVTP(const int nProcs) const
{
   std::ofstream myfile("eb.pvtp");

   if(myfile.is_open()) {
      myfile << "<?xml version=\"1.0\"?>\n";
      myfile << "<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      myfile << "<PPolyData GhostLevel=\"0\">\n";
      myfile << "<PPointData/>\n";
      myfile << "<PCellData/>\n";
      myfile << "<PPoints>\n";
      myfile << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
      myfile << "</PPoints>\n";

      for(int lc1 = 0; lc1 < nProcs; ++lc1) {
         std::stringstream ss;
         ss << std::setw(8) << std::setfill('0') << lc1;
         std::string clc1 = "eb_" + ss.str() + ".vtp";
         myfile << "<Piece Source=\"" << clc1 << "\"/>\n";
      }

      myfile << "</PPolyData>\n";
      myfile << "</VTKFile>\n";
      myfile.close();
   }
}
void EBToPVD::reorder_polygon(const std::vector<std::array<Real,3>>& lpoints,
      std::array<int,7>& lconnect,
      const std::array<Real,3>& lnormal)
{
   std::array<Real,3> center;
   center.fill(0.0);

   int longest = 2;
   if(Math::abs(lnormal[0]) > Math::abs(lnormal[1])) {
      if(Math::abs(lnormal[0]) > Math::abs(lnormal[2]))
         longest = 0;
   }
   else {
      if(Math::abs(lnormal[1]) > Math::abs(lnormal[2]))
         longest = 1;
   }

   for(int i = 1; i <= lconnect[0]; ++i) {
      center[0] += m_points[lconnect[i]][0];
      center[1] += m_points[lconnect[i]][1];
      center[2] += m_points[lconnect[i]][2];
   }
   center = {center[0]/lconnect[0], center[1]/lconnect[0], center[2]/lconnect[0]};

   int pi, pk;
   Real ref_angle, angle;
   if(longest == 0)
   {
      for(int i = 1; i <= lconnect[0]-1; ++i) {
         pi = lconnect[i];
         ref_angle = std::atan2(lpoints[pi][2]-center[2], lpoints[pi][1]-center[1]);
         for(int k = i+1; k <= lconnect[0]; ++k) {
            pk = lconnect[k];
            angle = std::atan2(lpoints[pk][2]-center[2], lpoints[pk][1]-center[1]);
            if(angle < ref_angle) {
               ref_angle = angle;
               lconnect[k] = pi;
               lconnect[i] = pk;
               pi = pk;
            }
         }
      }
   }
   else if(longest == 1) {
      for(int i = 1; i <= lconnect[0]-1; ++i) {
         pi = lconnect[i];
         ref_angle = std::atan2(lpoints[pi][0]-center[0], lpoints[pi][2]-center[2]);
         for(int k = i+1; k <= lconnect[0]; ++k) {
            pk = lconnect[k];
            angle = std::atan2(lpoints[pk][0]-center[0], lpoints[pk][2]-center[2]);
            if(angle < ref_angle) {
               ref_angle = angle;
               lconnect[k] = pi;
               lconnect[i] = pk;
               pi = pk;
            }
         }
      }
   }
   else if(longest == 2) {
      for(int i = 1; i <= lconnect[0]-1; ++i) {
         pi = lconnect[i];
         ref_angle = std::atan2(lpoints[pi][1]-center[1], lpoints[pi][0]-center[0]);
         for(int k = i+1; k <= lconnect[0]; ++k) {
            pk = lconnect[k];
            angle = std::atan2(lpoints[pk][1]-center[1], lpoints[pk][0]-center[0]);
            if(angle < ref_angle) {
               ref_angle = angle;
               lconnect[k] = pi;
               lconnect[i] = pk;
               pi = pk;
            }
         }
      }
   }
}

void EBToPVD::calc_hesse(Real& distance, std::array<Real,3>& n0, Real& p,
      const std::array<Real,3>& normal, const std::array<Real,3>& centroid) const
{
   Real sign_of_dist;

   // General equation of a plane: Ax + By + Cz + D = 0
   // here D := distance
   distance = -dot_product(normal, centroid);

   // Get the sign of the distance
   sign_of_dist = -distance / Math::abs(distance);

   // Get the Hessian form
   Real fac = sign_of_dist/dot_product(normal, normal);
   for(int idim = 0; idim < 3; ++idim) {
      n0[idim] = fac*normal[idim];
   }
   p = sign_of_dist*(-distance);
}

void EBToPVD::calc_alpha(std::array<Real,12>& alpha,
      const std::array<Real,3>& n0, Real p,
      const std::array<std::array<Real,3>,8>& vertex,
      const Real* dx) const
{
   // default (large) value
   std::fill(alpha.begin(), alpha.end(), 10.0);

   // Ray-xAxis intersection
   if(Math::abs(n0[0]) > std::numeric_limits<Real>::epsilon()) {
      alpha[0]  = (p - dot_product(n0,vertex[0]))/(n0[0]*dx[0]);
      alpha[2]  = (p - dot_product(n0,vertex[2]))/(n0[0]*dx[0]);
      alpha[8]  = (p - dot_product(n0,vertex[4]))/(n0[0]*dx[0]);
      alpha[10] = (p - dot_product(n0,vertex[6]))/(n0[0]*dx[0]);
   }

   // Ray-yAxis intersection
   if(Math::abs(n0[1]) > std::numeric_limits<Real>::epsilon()) {
      alpha[1]  = (p - dot_product(n0,vertex[1]))/(n0[1]*dx[1]);
      alpha[3]  = (p - dot_product(n0,vertex[0]))/(n0[1]*dx[1]);
      alpha[9]  = (p - dot_product(n0,vertex[5]))/(n0[1]*dx[1]);
      alpha[11] = (p - dot_product(n0,vertex[4]))/(n0[1]*dx[1]);
   }

   // Ray-zAxis intersection
   if(Math::abs(n0[2]) > std::numeric_limits<Real>::epsilon()) {
      alpha[4] = (p - dot_product(n0,vertex[0]))/(n0[2]*dx[2]);
      alpha[5] = (p - dot_product(n0,vertex[1]))/(n0[2]*dx[2]);
      alpha[6] = (p - dot_product(n0,vertex[3]))/(n0[2]*dx[2]);
      alpha[7] = (p - dot_product(n0,vertex[2]))/(n0[2]*dx[2]);
   }
}

void EBToPVD::calc_intersects(int& int_count, std::array<bool,12>& intersects_flags,
      const std::array<Real,12>& alpha) const
{
   int_count = 0;
   std::fill(intersects_flags.begin(), intersects_flags.end(), false);

   for(int lc1 = 0; lc1 < 12; ++lc1) {
      if(intersects(alpha[lc1])) {
         ++int_count;
         intersects_flags[lc1] = true;
      }
   }
}

void EBToPVD::print_points(std::ofstream& myfile) const
{
   myfile << "<Points>\n";
   myfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

   for(size_t lc1 = 0; lc1 < m_points.size(); ++lc1) {
      myfile << std::fixed << std::scientific
         << m_points[lc1][0] << " " << m_points[lc1][1] << " " << m_points[lc1][2] << "\n";
   }

   myfile << "</DataArray>\n";
   myfile << "</Points>\n";
}

void EBToPVD::print_connectivity(std::ofstream& myfile) const
{
   myfile << "<Polys>\n";
   myfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
   for(size_t lc1 = 0; lc1 < m_connectivity.size(); ++lc1) {
      for(int lc2 = 1; lc2 <= m_connectivity[lc1][0]; ++lc2) {
         myfile << " " << m_connectivity[lc1][lc2];
      }
      myfile << "\n";
   }
   myfile << "</DataArray>\n";

   myfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
   int lc2 = 0;
   for(size_t lc1 = 0; lc1 < m_connectivity.size(); ++lc1) {
      lc2 = lc2 + m_connectivity[lc1][0];
      myfile << " " << lc2;
   }
   myfile << "\n";
   myfile << "</DataArray>\n";

   myfile << "</Polys>\n";
}

void EBToPVD::EBGridCoverage(const int myID, const Real* problo, const Real* dx,
      const Box &bx, Array4<EBCellFlag const> const& flag)
{
   int lc1 = 0;

   const auto lo = lbound(bx);
   const auto hi = ubound(bx);

   std::array<int,3> nodes = {hi.x-lo.x + 1, hi.y-lo.y + 1, hi.z-lo.z + 1};
   std::array<int,3> low = {lo.x, lo.y, lo.z};

   for(int k = lo.z; k <= hi.z; ++k) {
      for(int j = lo.y; j <= hi.y; ++j) {
         for(int i = lo.x; i <= hi.x; ++i)
         {
            if(flag(i,j,k).isSingleValued())
               lc1 = lc1 + 1;
         }
      }
   };

   ++m_grid;
   if(lc1 == 0) return;

   std::stringstream ss;
   ss << std::setw(4) << std::setfill('0') << myID;
   std::string cID = ss.str();

   ss.str("");
   ss.clear();
   ss << std::setw(4) << std::setfill('0') << m_grid;
   std::string cgrid = ss.str();
   std::string fname = "eb_grid_" + cID + "_" + cgrid + ".vtr";

   std::ofstream myfile(fname);
   if(myfile.is_open()) {
      myfile.precision(6);
      myfile << "<?xml version=\"1.0\"?>\n";
      myfile << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      myfile << "<RectilinearGrid WholeExtent=\" 0 "
         << nodes[0] << " 0 " << nodes[1] << " 0 " << nodes[2] << "\">\n";
      myfile << "<Piece Extent=\" 0 "
         << nodes[0] << " 0 " << nodes[1] << " 0 " << nodes[2] << "\">\n";
      myfile << "<Coordinates>\n";

      for(int idim = 0; idim < 3; ++idim) {
         std::vector<Real> lines(nodes[idim]+1);
         Real grid_start = problo[idim] + low[idim]*dx[idim];
         for(int llc = 0; llc <= nodes[idim]; ++llc) {
            lines[llc] = grid_start + llc*dx[idim];
         }

         myfile << "<DataArray type=\"Float32\" format=\"ascii\" RangeMin=\""
            << std::fixed
            << lines[0] << "\" RangeMax=\"" << lines[nodes[idim]] << "\">\n";

         for(size_t llc = 0; llc < lines.size(); ++llc) {
            myfile << " " << lines[llc];
         }
         myfile << "\n";

         myfile << "</DataArray>\n";
      }

      myfile << "</Coordinates>\n";
      myfile << "</Piece>\n";
      myfile << "</RectilinearGrid>\n";
      myfile << "</VTKFile>\n";
   }

   myfile.close();
}

}
