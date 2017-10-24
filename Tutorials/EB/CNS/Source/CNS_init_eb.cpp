
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

#include <AMReX_ParmParse.H>

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
initialize_EBIS(const int max_level)
{
    BL_PROFILE("initialize_EBIS");

    if (!eb_initialized)
    {
        ParmParse pp;
        ParmParse ppa("amr");
        ParmParse ppg("geometry");
        Vector<int> n_cell;
        ppa.getarr("n_cell", n_cell, 0, SpaceDim);

        int max_grid_size = 0;
        ppa.query("max_grid_size", max_grid_size);

        // const Geometry& geom = parent->Geom(max_level);
        // const Box& finest_domain = geom.Domain();
        // const Real fine_dx = geom.CellSize()[0];
        // 
        //   ParmParse pp;

        Vector<int> ref_ratio;
        ppa.getarr("ref_ratio", ref_ratio, 0, max_level);
      
        IntVect ivlo = IntVect::Zero;
        IntVect ivhi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          ivhi[idir] = n_cell[idir]-1;
        }
        Box finest_domain(ivlo, ivhi);
        // < maxlev because there is one less refinement than number of levels
        for(int ilev = 0; ilev < max_level; ilev++)
        {
          finest_domain.refine(ref_ratio[ilev]);
        }
      
        Vector<Real> prob_lo, prob_hi;
        ppg.getarr("prob_lo", prob_lo, 0, SpaceDim);
        ppg.getarr("prob_hi", prob_hi, 0, SpaceDim);
        Real fine_dx = (prob_hi[0]-prob_lo[0])/finest_domain.size()[0];
        RealVect dxVec;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          dxVec[idir] = (prob_hi[idir]-prob_lo[idir])/finest_domain.size()[idir];
        }
        RealVect origin = RealVect::Zero;
  

        amrex::Print() << "defining EBIS with..." << endl;
        amrex::Print() << "finest domain = " << finest_domain << endl;
        amrex::Print() << "finest dx     = " << fine_dx       << endl;

        
        std::string geom_type;
        pp.get("geom_type", geom_type);
        if (geom_type == "all_regular")
        {
          //allregular
          amrex::Print() << "all regular geometry\n";
          AllRegularService regserv;
          AMReX_EBIS::instance()->define(finest_domain, origin, fine_dx, regserv, max_grid_size, max_level);
        }
        else if (geom_type == "flat_plate") 
        {
          amrex::Print() << "flat plate  geometry\n";
          Vector<Real>  platelovec(SpaceDim);
          Vector<Real>  platehivec(SpaceDim);

          int normalDir;
          Real plateLoc;
    
          pp.getarr("plate_lo", platelovec, 0, SpaceDim);
          pp.getarr("plate_hi", platehivec, 0, SpaceDim);
          pp.get("plate_location", plateLoc);
          pp.get("plate_normal", normalDir);

          RealVect plateLo, plateHi;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            plateLo[idir] = platelovec[idir];
            plateHi[idir] = platehivec[idir];
          }
          FlatPlateGeom flat_plate(normalDir, plateLoc, plateLo, plateHi);
          AMReX_EBIS::instance()->define(finest_domain, origin, fine_dx, flat_plate, max_grid_size, max_level);
        }
        else
        {
          std::unique_ptr<BaseIF> impfunc;
          
          if (geom_type == "ramp")
          {
            amrex::Print() << "ramp geometry\n";
            int upDir;
            int indepVar;
            Real startPt;
            Real slope;
            pp.get("up_dir",upDir);
            pp.get("indep_var",indepVar);
            pp.get("start_pt", startPt);
            pp.get("ramp_slope", slope);

            RealVect normal = RealVect::Zero;
            normal[upDir] = 1.0;
            normal[indepVar] = -slope;

            RealVect point = RealVect::Zero;
            point[upDir] = -slope*startPt;

            bool normalInside = true;

            impfunc.reset(static_cast<BaseIF*>(new PlaneIF(normal,point,normalInside)));
          }
          else if (geom_type == "anisotropic_ramp")
          {
            amrex::Print() << "anisotropic ramp geometry\n";
            int upDir;
            int indepVar;
            Real startPt;
            Real slope;
            pp.get("up_dir",upDir);
            pp.get("indep_var",indepVar);
            pp.get("start_pt", startPt);
            pp.get("ramp_slope", slope);

            RealVect normal = RealVect::Zero;
            normal[upDir] = 1.0;
            normal[indepVar] = -slope;

            RealVect point = RealVect::Zero;
            point[upDir] = -slope*startPt;

            bool normalInside = true;

            impfunc.reset(static_cast<BaseIF*>(new AnisotropicDxPlaneIF(normal,point,normalInside,dxVec)));
          }
          else if (geom_type == "sphere")
          {
            amrex::Print() << "sphere geometry\n";
            Vector<Real> centervec(SpaceDim);
            Real radius;
            pp.get(   "sphere_radius", radius);
            pp.getarr("sphere_center", centervec, 0, SpaceDim);
            RealVect center;
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }
            bool insideRegular = false;
            impfunc.reset(static_cast<BaseIF*>(new SphereIF(radius, center, insideRegular)));

          }
          else if (geom_type == "parabola")
          {

            amrex::Print() << "parabola geometry\n";
            Vector<PolyTerm> poly;

            PolyTerm mono;
            Real coef;
            IntVect powers;
            Real amplitude = 1;

            // y^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[1] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

#if BL_SPACEDIM==3
            // z^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[2] = 2;
            mono.coef   = coef;
            mono.powers = powers;
            poly.push_back(mono);
#endif
            // x term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[0] = 1;
            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            PolynomialIF mirror(poly,false);
            RealVect translation;
      
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              int finesize = finest_domain.size()[idir];
              translation[idir] = 0.5*finesize*fine_dx;
            }
            translation[0] = 0;

            TransformIF implicit(mirror);
            implicit.translate(translation);
            impfunc.reset(implicit.newImplicitFunction());

          }

          else if (geom_type == "parabola_and_sphere")
          {

            amrex::Print() << "parabola + sphere geometry\n";
            Vector<PolyTerm> poly;

            PolyTerm mono;
            Real coef;
            IntVect powers;
            Real amplitude = 1;

            // y^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[1] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

#if BL_SPACEDIM==3
            // z^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[2] = 2;
            mono.coef   = coef;
            mono.powers = powers;
            poly.push_back(mono);
#endif
            // x term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[0] = 1;
            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            PolynomialIF mirror(poly,false);
            RealVect translation;
      
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              int finesize = finest_domain.size()[idir];
              translation[idir] = 0.5*finesize*fine_dx;
            }
            RealVect center = translation;
            translation[0] = 0;

            TransformIF transform(mirror);
            transform.translate(translation);

            Real radius = 0.2*center[0];
            SphereIF sphere(radius, center, true);
            Vector<BaseIF*> funcs(2);
            funcs[0] = &transform;
            funcs[1] = &sphere;
            UnionIF implicit(funcs);
            impfunc.reset(implicit.newImplicitFunction());

          }
          else if (geom_type == "interior_box")
          {
            amrex::Print() << "box within a box" << endl;
            bool inside = true;  
            
            // A list of all the polygons as implicit functions
            Vector<BaseIF*> planes;
            // Process each polygon
            for(int idir =0; idir < SpaceDim; idir++)
            {
              Real domlen = fine_dx*finest_domain.size()[idir];
              RealVect pointLo = 1.5*fine_dx*BASISREALV(idir);
              RealVect pointHi = (domlen - 1.5*fine_dx)*BASISREALV(idir);
              RealVect normalLo = BASISREALV(idir);
              RealVect normalHi = -BASISREALV(idir);
              PlaneIF* planeLo = new PlaneIF(normalLo,pointLo,inside);
              PlaneIF* planeHi = new PlaneIF(normalHi,pointHi,inside);

              planes.push_back(planeLo);
              planes.push_back(planeHi);
            }
            
        

            IntersectionIF* intersectif = new IntersectionIF(planes);
            UnionIF* unionif = new UnionIF(planes);

            impfunc.reset(intersectif);
        
          }
      
          else if (geom_type == "polygon_revolution")
          {
            amrex::Print() << "creating geometry from polygon surfaces of revolution" << endl;
            bool insideRegular = false;
      
            // Data for polygons making up nozzle
            Vector<Vector<RealVect> > polygons;
      
      
            // For building each polygon
      
            int num_poly;
            RealVect translation;
      
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              int finesize = finest_domain.size()[idir];
              translation[idir] = 0.5*finesize*fine_dx;
            }
            Real scale = finest_domain.size()[0]*fine_dx;
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
      
          }
          else
          {
            //bogus which_geom
            amrex::Print() << " bogus geom_type input = " << geom_type << endl;
            amrex::Error("invalid inputs");
          }
        
          bool eb_verbosity = false;
          pp.query("eb_verbosity", eb_verbosity);
          GeometryShop gshop(*impfunc, eb_verbosity);
          AMReX_EBIS::instance()->define(finest_domain, origin, fine_dx, gshop, max_grid_size, max_level);
        }

        eb_initialized = true;

    } // End of finest level EB init - below gets done for all levels

}
