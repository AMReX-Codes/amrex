
#include <cmath>
#include <cstdio>
#include <iostream>

#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_ParmParse.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_PlaneIF.H"
#include "AMReX_SphereIF.H"
#include "AMReX_IntersectionIF.H"
#include "AMReX_UnionIF.H"
#include "CommonCode.H"
#include "DebugDump.H"
#include "AMReX_LatheIF.H"
#include "AMReX_ComplementIF.H"
#include "AMReX_TransformIF.H"
#include "WriteEBPlotFile.H"
#include "AMReX_RealVect.H"

namespace amrex
{
  UnionIF* makeCrossSection(const Vector<Vector<RealVect> >& a_polygons)
  {
    // The final result
    UnionIF* retval;
        
    // Get the number of polygons and make this inside of the domain
    // the inside of the polygons
    int numPolys = a_polygons.size();
    bool inside = true;
        
    // A list of all the polygons as implicit functions
    Vector<BaseIF*> polytopes;
    polytopes.resize(0);
        
    // Process each polygon
    for (int p = 0; p < numPolys; p++)
    {
      // All the half planes/spaces used to make a polygon
      Vector<BaseIF*> planes;
      planes.resize(0);
            
      // Get the current polygon (as a vector of points)
      const Vector<RealVect>& polygon = a_polygons[p];
            
      // Get the number of points in the polygon
      int numPts = polygon.size();
            
      // Process each pair of points
      for (int n = 0; n < numPts; n++)
      {
        // The normal and point is space used to specify each half plane/space
        RealVect normal(RealVect::Zero);
        RealVect point;
                
        // Set the normal remembering that the last point connects to the first
        // point.
        normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
        normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);
                
        point = polygon[n];
                
        // Generate the appropriate half plane/space (as an implicit function)
        PlaneIF* plane;
        plane = new PlaneIF(normal,point,inside);
                
        // Save the result
        planes.push_back(plane);
      }
            
      // Intersect all the half planes/spaces to create an implicit function
      // that represents the polygon
      IntersectionIF* polygonIF = new IntersectionIF(planes);
            
      polytopes.push_back(polygonIF);
    }
        
    // Union all the polygon implicit functions to get the implicit function
    // returned
    retval = new UnionIF(polytopes);
        
    return retval;
  }

  void
  makeEBIS(const Box&       a_domain,
           const RealVect&  a_dx)
  {
    amrex::Print() << "creating geometry from polygon surfaces of revolution" << endl;
    bool insideRegular = false;
      
    std::unique_ptr<BaseIF> impfunc;
    // Data for polygons making up nozzle
    Vector<Vector<RealVect> > polygons;
      
      
    // For building each polygon
      
    int num_poly;
    RealVect translation;
      
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int finesize = a_domain.size()[idir];
      translation[idir] = 0.5*finesize*a_dx[idir];
    }
    Real scale = a_domain.size()[0]*a_dx[0];
    ParmParse pp;
    pp.get("num_poly", num_poly);
    // Nothing initially
    polygons.resize(num_poly);
    amrex::Print() << "num poly = " << num_poly << endl;
    for(int ipoly = 0; ipoly < num_poly; ipoly++)
    {
      string nptsstr = "poly_" + std::to_string(ipoly) + "_num_pts";
      int num_pts;
      pp.get(nptsstr.c_str(), num_pts);
      Vector<RealVect> polygon(num_pts);
      for(int ipt = 0; ipt < num_pts; ipt++)
      {
        RealVect point(RealVect::Zero);
        string    pointstr = "poly_" + std::to_string(ipoly) + "_point_" + std::to_string(ipt);
        Vector<Real> vecpt;
        pp.getarr(pointstr.c_str(), vecpt,  0, SpaceDim);
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          point[idir] = vecpt[idir] ;
        }
        //now scale by the size of the domain
        point *= scale;
        polygon[ipt] = point;
      }
            
      amrex::Print() << "scaled poly" << ipoly << " = " << endl;
      for(int ipt = 0; ipt < num_pts; ipt++)
      {
        amrex::Print() << polygon[ipt] << " ";
      }
      amrex::Print() << endl;
      polygons[ipoly] = polygon;
    }
      
        
    // Make the vector of (convex) polygons (vectors of points) into a union
    // of convex polygons, each made from the intersection of a set of half
    // planes/spaces - all represented by implicit functions.
    UnionIF* crossSection = makeCrossSection(polygons);
    //SmoothUnion* crossSection = makeSmoothCrossSection(polygons, fine_dx);
      
    if (SpaceDim == 2)
    {
      // In 2D use "as is"
      
      // Complement if necessary
      if(!insideRegular)
      {
        ComplementIF insideOut(*crossSection);
        impfunc.reset(insideOut.newImplicitFunction());
      }
      else
      {
        impfunc.reset(crossSection->newImplicitFunction());
      }

    }
    else
    {
      // In 3D rotate about the z-axis and complement if necessary
      LatheIF lathe(*crossSection,insideRegular);
      //we are starting around the y axis so we need to translate
      //over to the center 
      
      translation[2] = 0;
      TransformIF implicit(lathe);
      implicit.translate(translation);
      impfunc.reset(implicit.newImplicitFunction());
    }
 
    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int biggridsize;

    pp.get("max_grid_size", biggridsize);
    GeometryShop workshop(*impfunc);
    int ebmaxcoarsen = 0;
    RealVect origin = RealVect::Zero;
    ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, ebmaxcoarsen);
  }
/************/
  void
  makeAndOutputEBIS()
  {
    //make layouts == domain
    Box domainBox;
    RealVect dx;
    getFinestDomain(domainBox, dx);

    BoxArray ba(domainBox);
    DistributionMapping dm(ba);

    makeEBIS(domainBox,  dx);
    EBLevelGrid eblg(ba, dm, domainBox, 2);
    
    MultiFab data(ba, dm, 1, 0);
    setDataToSomething(data, eblg);
    
    std::string filename = string("viva_la_revolucion.") + convertInt(SpaceDim) + "d.plt";
    WriteSingleLevelEBPlotFile(filename, data, eblg, Vector<string>());

  }
}
/************/
/************/
int
main(int argc, char* argv[])
{
  int retval = 0;
  amrex::Initialize(argc,argv);

  amrex::makeAndOutputEBIS();

  amrex::Finalize();
  return retval;
}
