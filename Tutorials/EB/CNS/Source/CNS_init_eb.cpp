
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
initialize_EBIS(const Geometry& geom, const int max_level)
{
    BL_PROFILE("initialize_EBIS");

    if (!eb_initialized)
    {
        ParmParse pp;
        ParmParse ppa("amr");
        ParmParse ppg("geometry");
        ParmParse ppeb2("eb2");
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
        ppeb2.get("geom_type", geom_type);
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
    
          ppeb2.getarr("plate_lo", platelovec, 0, SpaceDim);
          ppeb2.getarr("plate_hi", platehivec, 0, SpaceDim);
          ppeb2.get("plate_location", plateLoc);
          ppeb2.get("plate_normal", normalDir);

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
            ppeb2.get("up_dir",upDir);
            ppeb2.get("indep_var",indepVar);
            ppeb2.get("start_pt", startPt);
            ppeb2.get("ramp_slope", slope);

            RealVect normal = RealVect::Zero;
            normal[upDir] = 1.0;
            normal[indepVar] = -slope;

            RealVect point = RealVect::Zero;
            point[upDir] = -slope*startPt;

            bool normalInside = true;

            impfunc.reset(static_cast<BaseIF*>(new PlaneIF(normal,point,normalInside)));
          }
          else if (geom_type == "plane")//"ramp_normal_point")
          {
            amrex::Print() << "ramp geometry using normal and point directly \n";
            RealVect normal;
            RealVect point;
            Vector<Real> pointvec; 
            Vector<Real> normalvec;
            int inside;
            ppeb2.getarr("ramp_normal", normalvec, 0, SpaceDim);
            ppeb2.getarr("ramp_point" ,  pointvec, 0, SpaceDim);
            ppeb2.get("ramp_inside", inside);
            bool normalInside = (inside == 0);
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              point[idir]  =  pointvec[idir];
              normal[idir] = normalvec[idir];
            }


            impfunc.reset(static_cast<BaseIF*>(new PlaneIF(normal,point,normalInside)));
          }
          else if (geom_type == "anisotropic_ramp")
          {
            amrex::Print() << "anisotropic ramp geometry\n";
            int upDir;
            int indepVar;
            Real startPt;
            Real slope;
            ppeb2.get("up_dir",upDir);
            ppeb2.get("indep_var",indepVar);
            ppeb2.get("start_pt", startPt);
            ppeb2.get("ramp_slope", slope);

            RealVect normal = RealVect::Zero;
            normal[upDir] = 1.0;
            normal[indepVar] = -slope;

            RealVect point = RealVect::Zero;
            point[upDir] = -slope*startPt;

            bool normalInside = true;

            impfunc.reset(static_cast<BaseIF*>(new AnisotropicDxPlaneIF(normal,point,normalInside,dxVec)));
          }

          else if (geom_type == "anisotropic_sphere")
          {
            amrex::Print() << "anisotropic sphere geometry\n";
            Vector<Real> centervec(SpaceDim);
            Real radius;
            ppeb2.get(   "sphere_radius", radius);
            ppeb2.getarr("sphere_center", centervec, 0, SpaceDim);
            RealVect center;
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }
            bool insideRegular = false;
            shared_ptr<BaseIF> baseif(new SphereIF(radius, center, insideRegular));

            impfunc.reset(static_cast<BaseIF*>(new AnisotropicIF(baseif,dxVec)));
          }
          else if (geom_type == "sphere")
          {
            amrex::Print() << "sphere geometry\n";
            Vector<Real> centervec(SpaceDim);
            Real radius;
            ppeb2.get(   "sphere_radius", radius);
            ppeb2.getarr("sphere_center", centervec, 0, SpaceDim);
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

#if AMREX_SPACEDIM==3
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

#if AMREX_SPACEDIM==3
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
//            UnionIF* unionif = new UnionIF(planes);

            impfunc.reset(intersectif);
        
          }
      
          else if (geom_type == "double_ramp")
          {
            amrex::Print() << "box within a box" << endl;
            bool inside = true;  
            
            // A list of all the polygons as implicit functions
            Vector<BaseIF*> planes;
            // Process each polygon
            int idir = 1;
            {
//              Real domlen = fine_dx*finest_domain.size()[idir];
              RealVect pointLo = 10.0*BASISREALV(idir);
              RealVect pointHi = 19.0*BASISREALV(idir);
              RealVect normalLo = BASISREALV(idir);
              RealVect normalHi = -BASISREALV(idir);
              PlaneIF* planeLo = new PlaneIF(normalLo,pointLo,inside);
              PlaneIF* planeHi = new PlaneIF(normalHi,pointHi,inside);

              planes.push_back(planeLo);
              planes.push_back(planeHi);
            }
            
        

            IntersectionIF* intersectif = new IntersectionIF(planes);
//            UnionIF* unionif = new UnionIF(planes);

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
	  else if(geom_type=="combustor")
	  {
		amrex::Print() << "Direct Point input for combustor! \n"; 
		bool insideRegular = false; 
		Real hix = geom.ProbHi(0); 
		Real hiy = geom.ProbHi(1); 
		Real loy = geom.ProbLo(1); 

		//Data for constructing combustor
		Vector<Vector<RealVect> > polygons; 
		polygons.resize(4); //3 polygons "Pipe, Ramp, Flattened Corner and Far Wall"  

		//=============== Far Wall =====================================
		//User gives location, from this generate plane polygon 
		Vector<RealVect> farwall(4); //four points 
		Real fwl; 
		ppeb2.get("far_wall_loc", fwl);
		//gets x location of far wall  
		RealVect fwpnt0 = RealVect::Zero; 
		RealVect fwpnt1 = RealVect::Zero; 
		RealVect fwpnt2 = RealVect::Zero; 
		RealVect fwpnt3 = RealVect::Zero; 
		fwpnt0[0] = fwl; //Should be physical location 
		fwpnt0[1] = loy; //Scaled for convenience 
		fwpnt1[0] = 1.0*hix; 
		fwpnt1[1] = loy; 
		fwpnt2[0] = 1.0*hix; 
		fwpnt2[1] = 1.0*hiy; 
		fwpnt3[0] = fwl; 
		fwpnt3[1] = 1.0*hiy; 
		farwall[0] = fwpnt0; 
		farwall[1] = fwpnt1; 
		farwall[2] = fwpnt2; 
		farwall[3] = fwpnt3; 
		polygons[0] = farwall; 

		//=============== Pipe Section ====================================
		//Require 4 Points, much like the Farwall. However, user passes in "high" and "low"
		//Utilizing high and low, a plane polygon is derived
		Vector<RealVect> pipe(4);  
                Vector<Real> pipelovec, pipehivec;
                ppeb2.getarr("pipe_lo", pipelovec, 0, SpaceDim);
                ppeb2.getarr("pipe_hi", pipehivec, 0, SpaceDim);
		RealVect pipelo1(RealVect::Zero), pipelo2(RealVect::Zero); 
		RealVect pipehi1(RealVect::Zero), pipehi2(RealVect::Zero); 
		pipelo1[0] = pipelovec[0]; 
		pipelo1[1] = pipelovec[1]; 
		pipelo2[0] = pipehivec[0]; 
		pipelo2[1] = pipelovec[1]; 
	
		pipehi1[0] = pipehivec[0]; 
		pipehi1[1] = pipehivec[1]; 
		pipehi2[0] = pipelovec[0]; 
		pipehi2[2] = pipehivec[1]; 
		
		pipe[0] = pipelo1; 
		pipe[1] = pipelo2; 
		pipe[2] = pipehi1; 
		pipe[3] = pipehi2; 


		//================ Ramp Section ===============================
		//Ramp contains 3 planes and a flattened corner, we will need 6 points to 
		//construct it as a polygon. The user supplies the 3 plane points, and one normal vec
		Vector<RealVect> ramp(5); 
				 
                Vector<Real> pl1ptvec, pl2ptvec, pl2nrmlvec, pl3ptvec;
                RealVect ramp0(RealVect::Zero), ramp1(RealVect::Zero); 
		RealVect ramp2(RealVect::Zero), ramp3(RealVect::Zero);
		RealVect ramp4(RealVect::Zero), ramp5(RealVect::Zero);
                ppeb2.getarr("ramp_plane1_point", pl1ptvec, 0, SpaceDim);
                ppeb2.getarr("ramp_plane2_point", pl2ptvec, 0, SpaceDim);
                ppeb2.getarr("ramp_plane2_normal", pl2nrmlvec, 0, SpaceDim);
                ppeb2.getarr("ramp_plane3_point", pl3ptvec, 0, SpaceDim);
		Real rampy = -pl2nrmlvec[0]/pl2nrmlvec[1]*(pipehi1[0] - pl2ptvec[0]) + pl2ptvec[1]; 


		ramp0[0] = pipehi1[0]; 
		ramp0[1] = pipelo1[1];
		ramp1[0] = fwl; 
		ramp1[1] = pipelo1[1]; 
		ramp2[0] = fwl;  
		ramp2[1] = pl2ptvec[1]; 
		ramp3[0] = pl2ptvec[0]; 
		ramp3[1] = pl2ptvec[1]; 
		ramp4[0] = pipehi1[0]; 
		ramp4[1] = rampy; 
		
		ramp[0] = ramp0; 
		ramp[1] = ramp1; 
		ramp[2] = ramp2; 
		ramp[3] = ramp3; 
		ramp[4] = ramp4; 
		
		polygons[1] = ramp; 

		//================ Flattened Corner ==============================
		// This flattens the corner between the pipe and ramp, allowing for better 
		// Cell cuts. 

		Vector<RealVect> flatcnr(4); 
		Real dx = std::max(geom.CellSize()[0],std::max(geom.CellSize()[1], geom.CellSize()[2]));
                Real ydx;

                //Derive location for flat corner from plane2/pipe intersection. 
                if(dx < ((geom.ProbHi(0) - geom.ProbLo(0))/32.))
                        ydx = -pl2nrmlvec[0]/pl2nrmlvec[1]*(pipehi1[0]+ 4./3.*dx
			      - pl2ptvec[0]) + pl2ptvec[1];
                else
                        ydx = pipehi1[1];

		RealVect fltpt0(RealVect::Zero), fltpt1(RealVect::Zero); 
		RealVect fltpt2(RealVect::Zero), fltpt3(RealVect::Zero); 

		fltpt0 = pipelo2; 
		fltpt1[0] = pipelo2[0] + dx; 
		fltpt1[1] = pipelo2[1]; 
		fltpt2[0] = pipelo2[0] + dx; 
		fltpt2[1] = ydx; 
		fltpt3[0] = pipelo2[0]; 
		fltpt3[1] = ydx; 

		flatcnr[0] = fltpt0; 
		flatcnr[1] = fltpt1; 
		flatcnr[2] = fltpt2; 
		flatcnr[3] = fltpt3; 


		polygons[2] = flatcnr; 
		polygons[3] = pipe; //Polygons in order of farthest to nearest (to the origin)

		UnionIF* crossSection = makeCrossSection(polygons); 

		if(SpaceDim == 2) 
		{
		  //In 2D use "as is" 
		  impfunc.reset(crossSection->newImplicitFunction()); 
		} 
		else{
		  //In 3D rotate about the z-axis
			LatheIF lathe(*crossSection, insideRegular);  
		  //We are starting around the y-axis so we need to translate 
		  //Over to the center 

			Vector<Real> trans_vec; 
		 	ppeb2.getarr("lathe_trans_vec", trans_vec, 0, SpaceDim); 
		 	RealVect translation; 
		 	for(int idir = 0; idir < SpaceDim; idir++) 
				translation[idir] = trans_vec[idir]; 
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
