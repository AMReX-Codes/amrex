
#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_AnisotropicDxPlaneIF.H>
#include <AMReX_AnisotropicIF.H>

#include <AMReX_ParmParse.H>

#include "MyTest.H"

using namespace amrex;

namespace {
    bool eb_initialized = false;

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

}

void
MyTest::initializeEBIS ()
{
    BL_PROFILE("initialize_EBIS");

    if (!eb_initialized)
    {
        ParmParse pp;

        const Box& finest_domain = geom.back().Domain();
        const Real* dx = geom.back().CellSize();
        const RealVect dxVec {AMREX_D_DECL(dx[0], dx[1], dx[2])};
        RealVect origin = RealVect::Zero;
        
        amrex::Print() << "defining EBIS with ... ";
        
        if (geom_type == "all_regular")
        {
            //allregular
            amrex::Print() << "all regular geometry\n";
            AllRegularService regserv;
            AMReX_EBIS::instance()->define(finest_domain, origin, dxVec[0], regserv,
                                           max_grid_size, max_level);
        }
        else
        {
            std::unique_ptr<BaseIF> impfunc;
            
            if (geom_type == "ramp_normal_point")
            {
                amrex::Print() << "ramp geometry using normal and point directly \n";
                Vector<Real> pointvec {AMREX_D_DECL(0.,0.,0.)};
                Vector<Real> normalvec{AMREX_D_DECL(-0.1, 1.0, 0.0)};
                int inside = 1;
                pp.queryarr("ramp_normal", normalvec, 0, SpaceDim);
                pp.queryarr("ramp_point" ,  pointvec, 0, SpaceDim);
                pp.query("ramp_inside", inside);
                bool normalInside = (inside == 0);
                RealVect normal{AMREX_D_DECL(normalvec[0],normalvec[1],normalvec[2])};
                RealVect point {AMREX_D_DECL( pointvec[0], pointvec[1], pointvec[2])};

                impfunc.reset(static_cast<BaseIF*>(new PlaneIF(normal,point,normalInside)));
            }
            else
            {
                amrex::Print() << " bogus geom_type input = " << geom_type << "\n";
                amrex::Error("invalid inputs");
            }
        
            bool eb_verbosity = true;
            pp.query("eb_verbosity", eb_verbosity);
            GeometryShop gshop(*impfunc, eb_verbosity);
            AMReX_EBIS::instance()->define(finest_domain, origin, dxVec[0], gshop,
                                           max_grid_size, max_level+100);
        }

        eb_initialized = true;

    }
}
