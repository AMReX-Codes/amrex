#include <AMReX_EBToPVD.H>
#include <AMReX_BLassert.H>
#include <AMReX_Dim3.H>

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

namespace amrex {

namespace {
constexpr Real intersection_tol = Real(16)*std::numeric_limits<Real>::epsilon();

amrex::Real dot_product(const std::array<amrex::Real,3>& a, const std::array<amrex::Real,3>& b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

}

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

               Real adx = (axm-axp) * dx[1] * dx[2];
               Real ady = (aym-ayp) * dx[0] * dx[2];
               Real adz = (azm-azp) * dx[0] * dx[1];

               Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);
               Real apnorminv = Real(1.0)/apnorm;

               std::array<Real,3> normal, centroid;
               std::array<std::array<Real,3>,8> vertex;

               normal[0] = adx * apnorminv;
               normal[1] = ady * apnorminv;
               normal[2] = adz * apnorminv;

               // Work relative to vertex 0 to avoid cancellation in plane distances.
               centroid[0] = (bcent(i,j,k,0) + Real(0.5))*dx[0];
               centroid[1] = (bcent(i,j,k,1) + Real(0.5))*dx[1];
               centroid[2] = (bcent(i,j,k,2) + Real(0.5))*dx[2];

               vertex[0] = {0,     0,     0    };
               vertex[1] = {dx[0], 0,     0    };
               vertex[2] = {0,     dx[1], 0    };
               vertex[3] = {dx[0], dx[1], 0    };
               vertex[4] = {0,     0,     dx[2]};
               vertex[5] = {dx[0], 0,     dx[2]};
               vertex[6] = {0,     dx[1], dx[2]};
               vertex[7] = {dx[0], dx[1], dx[2]};

               // NOTE: this seems to be unnecessary:
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
               std::array<std::array<Real,3>,12> apoints{};

               calc_hesse(distance, n0, p, normal, centroid);
               calc_intersects(count, apoints, n0, p, vertex, dx);

               // If the number of facet "contained" in does not describe a facet:
               // ... I.e. there's less than 3 (not even a triangle) or more than 6
               // ... (I have no idea what that is):
               //   => Move the centroid a little back and forth along the normal
               //      to see if that makes a difference:

               if((count < 3) || (count > 6)) {
                  int count_d;
                  Real p_d;
                  std::array<Real,3> n0_d;
                  std::array<std::array<Real,3>,12> apoints_d{};

                  Real tol = std::min({dx[0], dx[1], dx[2]})/Real(100.); // bit of a fudge factor

                  std::array<Real,3> centroid_d;
                  for(int idim = 0; idim < 3; ++idim) {
                     centroid_d[idim] = centroid[idim] + tol*normal[idim];
                  }

                  calc_hesse(distance, n0_d, p_d, normal, centroid_d);
                  calc_intersects(count_d, apoints_d, n0_d, p_d, vertex, dx);
                  if((count_d >= 3) && (count_d <= 6)) {
                     count = count_d;
                     apoints = apoints_d;
                  }

                  for(int idim = 0; idim < 3; ++idim) {
                     centroid_d[idim] = centroid[idim] - tol*normal[idim];
                  }

                  calc_hesse(distance, n0_d, p_d, normal, centroid_d);
                  calc_intersects(count_d, apoints_d, n0_d, p_d, vertex, dx);
                  if((count_d >= 3) && (count_d <= 6)) {
                     count = count_d;
                     apoints = apoints_d;
                  }
               }
               // I know this was a bit of a hack, but it's the only way I prevent
               // missing facets...

               if((count >=3) && (count <=6)) {
                  m_connectivity.push_back({0,0,0,0,0,0,0});

                  for(int n = 0; n < count; ++n) {
                     m_points.push_back(apoints[n]);
                     m_connectivity.back()[0] = n+1;
                     m_connectivity.back()[n+1] = static_cast<int>(m_points.size()-1);
                  }

                  reorder_polygon(m_points, m_connectivity.back(), n0);

                  // If the plane lies on a cell face, the polygon covers the whole
                  // face, but only part of it is wall.  Trim it with a straight line
                  // on the side away from bcent so that its area matches the EB area.
                  auto& conn = m_connectivity.back();
                  int const np = conn[0];
                  std::array<std::array<Real,3>,8> poly;
                  std::array<Real,3> pcent = {0, 0, 0};
                  for(int n = 0; n < np; ++n) {
                     poly[n] = m_points[conn[n+1]];
                     for(int idim = 0; idim < 3; ++idim) { pcent[idim] += poly[n][idim]/Real(np); }
                  }
                  auto area = [&] (std::array<std::array<Real,3>,8> const& pts, int n) {
                     std::array<Real,3> s = {0, 0, 0};
                     for(int a = 0; a < n; ++a) {
                        auto const& u = pts[a];
                        auto const& v = pts[(a+1)%n];
                        s[0] += u[1]*v[2] - u[2]*v[1];
                        s[1] += u[2]*v[0] - u[0]*v[2];
                        s[2] += u[0]*v[1] - u[1]*v[0];
                     }
                     return Real(0.5)*std::abs(dot_product(s, n0));
                  };
                  std::array<Real,3> m;
                  Real const mc = dot_product(n0, centroid) - dot_product(n0, pcent);
                  for(int idim = 0; idim < 3; ++idim) {
                     m[idim] = centroid[idim] - pcent[idim] - mc*n0[idim];
                  }
                  Real const mnorm = std::sqrt(dot_product(m, m));
                  if(area(poly, np) > apnorm*Real(1.01) &&
                     mnorm > Real(1.e-3)*std::min({dx[0], dx[1], dx[2]})) {
                     for(auto& x : m) { x /= mnorm; }
                     // keep the part with dot(m,x) >= t
                     std::array<std::array<Real,3>,8> clipped;
                     int nc = 0;
                     auto clip = [&] (Real t) {
                        nc = 0;
                        for(int a = 0; a < np; ++a) {
                           auto const& u = poly[a];
                           auto const& v = poly[(a+1)%np];
                           Real const du = dot_product(m, u) - t;
                           Real const dv = dot_product(m, v) - t;
                           if(du >= 0) { clipped[nc++] = u; }
                           if((du >= 0) != (dv >= 0)) {
                              Real const f = du/(du-dv);
                              for(int idim = 0; idim < 3; ++idim) {
                                 clipped[nc][idim] = u[idim] + f*(v[idim]-u[idim]);
                              }
                              ++nc;
                           }
                        }
                        return (nc >= 3) ? area(clipped, nc) : Real(0);
                     };
                     Real tlo = std::numeric_limits<Real>::max();
                     Real thi = std::numeric_limits<Real>::lowest();
                     for(int n = 0; n < np; ++n) {
                        tlo = std::min(tlo, dot_product(m, poly[n]));
                        thi = std::max(thi, dot_product(m, poly[n]));
                     }
                     for(int iter = 0; iter < 60; ++iter) {
                        Real const t = Real(0.5)*(tlo+thi);
                        if(clip(t) > apnorm) { tlo = t; } else { thi = t; }
                     }
                     clip(tlo);
                     if(nc >= 3 && nc <= 6) {
                        for(int n = 0; n < np; ++n) { m_points.pop_back(); }
                        conn = {0,0,0,0,0,0,0};
                        for(int n = 0; n < nc; ++n) {
                           m_points.push_back(clipped[n]);
                           conn[0] = n+1;
                           conn[n+1] = static_cast<int>(m_points.size()-1);
                        }
                     }
                  }

                  // Translate the completed polygon to physical coordinates.
                  std::array<Real,3> const origin = {
                     problo[0] + static_cast<Real>(i)*dx[0],
                     problo[1] + static_cast<Real>(j)*dx[1],
                     problo[2] + static_cast<Real>(k)*dx[2]};
                  for(int n = 1; n <= conn[0]; ++n) {
                     for(int idim = 0; idim < 3; ++idim) {
                        m_points[conn[n]][idim] += origin[idim];
                     }
                  }
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
      myfile << "<Piece NumberOfPoints=\"" << m_points.size() << "\" NumberOfVerts=\"0\" " // NOLINT
         << "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\" " // NOLINT
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

void EBToPVD::WritePVTP(const int nProcs)
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
   if(std::abs(lnormal[0]) > std::abs(lnormal[1])) {
      if(std::abs(lnormal[0]) > std::abs(lnormal[2])) {
         longest = 0;
      }
   }
   else {
      if(std::abs(lnormal[1]) > std::abs(lnormal[2])) {
         longest = 1;
      }
   }

   for(int i = 1; i <= lconnect[0]; ++i) {
      center[0] += m_points[lconnect[i]][0];
      center[1] += m_points[lconnect[i]][1];
      center[2] += m_points[lconnect[i]][2];
   }
   center = {center[0]/static_cast<Real>(lconnect[0]),
             center[1]/static_cast<Real>(lconnect[0]),
             center[2]/static_cast<Real>(lconnect[0])};

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
      const std::array<Real,3>& normal, const std::array<Real,3>& centroid)
{
   Real sign_of_dist;

   // General equation of a plane: Ax + By + Cz + D = 0
   // here D := distance
   distance = -dot_product(normal, centroid);

   // Get the sign of the distance
   sign_of_dist = std::copysign(Real(1), -distance);

   // normal is already normalized; apply the same sign to it and the offset.
   for(int idim = 0; idim < 3; ++idim) {
      n0[idim] = sign_of_dist*normal[idim];
   }
   p = sign_of_dist*(-distance);
}

void EBToPVD::calc_intersects(int& int_count,
      std::array<std::array<Real,3>,12>& points,
      const std::array<Real,3>& n0, Real p,
      const std::array<std::array<Real,3>,8>& vertex, const Real* dx)
{
   int_count = 0;
   std::array<Real,8> distance;
   // Bound roundoff in the plane equation before dividing by an edge's
   // normal component, which can amplify it for nearly parallel edges.
   Real const tol = intersection_tol *
      (std::abs(n0[0])*dx[0] + std::abs(n0[1])*dx[1] + std::abs(n0[2])*dx[2]);
   for(int v = 0; v < 8; ++v) {
      distance[v] = dot_product(n0, vertex[v]) - p;
      if(std::abs(distance[v]) <= tol) {
         distance[v] = Real(0);
         points[int_count++] = vertex[v];
      }
   }

   constexpr std::array<std::array<int,2>,12> edge_vertices = {{
      {0,1}, {1,3}, {2,3}, {0,2}, {0,4}, {1,5},
      {3,7}, {2,6}, {4,5}, {5,7}, {6,7}, {4,6}
   }};
   for(auto const& edge : edge_vertices) {
      int const a = edge[0];
      int const b = edge[1];
      // On-plane corners were already added once, including the endpoints
      // of edges lying in the plane. Only strict sign changes add new points.
      if((distance[a] < 0 && distance[b] > 0) ||
         (distance[a] > 0 && distance[b] < 0)) {
         Real const alpha = distance[a]/(distance[a] - distance[b]);
         for(int idim = 0; idim < 3; ++idim) {
            points[int_count][idim] = vertex[a][idim] +
               alpha*(vertex[b][idim] - vertex[a][idim]);
         }
         ++int_count;
      }
   }
}

void EBToPVD::print_points(std::ofstream& myfile) const
{
   myfile << "<Points>\n";
   myfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

   for(const auto & m_point : m_points) {
      myfile << std::fixed << std::scientific
         << m_point[0] << " " << m_point[1] << " " << m_point[2] << "\n";
   }

   myfile << "</DataArray>\n";
   myfile << "</Points>\n";
}

void EBToPVD::print_connectivity(std::ofstream& myfile) const
{
   myfile << "<Polys>\n";
   myfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
   for(const auto & lc1 : m_connectivity) {
      for(int lc2 = 1; lc2 <= lc1[0]; ++lc2) {
         myfile << " " << lc1[lc2];
      }
      myfile << "\n";
   }
   myfile << "</DataArray>\n";

   myfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
   int lc2 = 0;
   for(const auto & lc1 : m_connectivity) {
      lc2 = lc2 + lc1[0];
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
            if(flag(i,j,k).isSingleValued()) {
               lc1 = lc1 + 1;
            }
         }
      }
   };

   ++m_grid;
   if(lc1 == 0) { return; }

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
         Real grid_start = problo[idim] + static_cast<Real>(low[idim])*dx[idim];
         for(int llc = 0; llc <= nodes[idim]; ++llc) {
            lines[llc] = grid_start + static_cast<Real>(llc)*dx[idim];
         }

         myfile << "<DataArray type=\"Float32\" format=\"ascii\" RangeMin=\"" // NOLINT
            << std::fixed
            << lines[0] << "\" RangeMax=\"" << lines[nodes[idim]] << "\">\n";

         for(auto line : lines) {
            myfile << " " << line;
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
