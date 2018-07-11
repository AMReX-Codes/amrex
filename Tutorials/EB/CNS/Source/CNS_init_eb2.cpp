
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>

#include <AMReX_ParmParse.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
                const int max_coarsening_level)
{
    BL_PROFILE("initializeEB2");

    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);
    
    if (geom_type == "combustor")
    {
	        Real fwl; 
		pp.get("far_wall_loc",fwl);
 
	        EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
        	                     {AMREX_D_DECL(1.  ,0.,0.)});

		Vector<Real> pl1ptvec, pl2ptvec, pl2nrmlvec, pl3ptvec; 
		RealArray pl1pt, pl2pt, pl2nrm, pl3pt; 
		pp.getarr("ramp_plane1_point", pl1ptvec, 0, SpaceDim);
		pp.getarr("ramp_plane2_point", pl2ptvec, 0, SpaceDim); 
		pp.getarr("ramp_plane2_normal", pl2nrmlvec, 0, SpaceDim); 
		pp.getarr("ramp_plane3_point", pl3ptvec, 0, SpaceDim); 
		for(int idir=0; idir<SpaceDim; idir++)
		{
			pl1pt[idir]  = pl1ptvec[idir]; 
			pl2pt[idir]  = pl2ptvec[idir]; 
			pl2nrm[idir] = pl2nrmlvec[idir]; 
			pl3pt[idir]  = pl3ptvec[idir];
		}

	        auto ramp = EB2::makeIntersection(EB2::PlaneIF(pl1pt,
                                                 {AMREX_D_DECL(0.  , -1. , 0.)}),
	                                         EB2::PlaneIF(pl2pt,pl2nrm),
                                                 EB2::PlaneIF(pl3pt,
                                                 {AMREX_D_DECL(1.  , 0.  , 0.)}));

		Vector<Real> pipelovec, pipehivec; 
		pp.getarr("pipe_lo", pipelovec, 0, SpaceDim); 
		pp.getarr("pipe_hi", pipehivec, 0, SpaceDim); 
		RealArray pipe_lo, pipe_hi; 
		for(int idir=0; idir<SpaceDim; idir++)
		{
			pipe_lo[idir]  = pipelovec[idir]; 
			pipe_hi[idir]  = pipehivec[idir]; 
		}


	        EB2::BoxIF pipe(pipe_lo, pipe_hi, false);
		Real dx = std::max(geom.CellSize()[0],std::max(geom.CellSize()[1], geom.CellSize()[2])); 
		Real ydx; 

		//Derive location for flat corner from ramp intersection. 
		if(dx < ((geom.ProbHi(0) - geom.ProbLo(0))/32.))
			ydx = -pl2nrm[0]/pl2nrm[1]*(pipe_hi[0]+ 4./3.*dx - pl2pt[0]) + pl2pt[1]; 
		else 
			ydx = pipe_hi[1];

	        EB2::BoxIF flat_corner({AMREX_D_DECL(pipe_lo[0]-0.0001, -geom.ProbHi(1), -100.0)},
                               {AMREX_D_DECL(geom.ProbHi(0), ydx, 100.0)}, false);

		Vector<Real> lathe_trans_vec; 
		RealArray lathe_trans; 

		pp.getarr("lathe_trans_vec", lathe_trans_vec, 0, SpaceDim); 
		for(int idir=0; idir<SpaceDim; idir++) 
			lathe_trans[idir] = lathe_trans_vec[idir]; 


	        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);
	        auto pr = EB2::translate(EB2::lathe(polys), lathe_trans);

	        auto gshop = EB2::makeShop(pr);
	        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
    }
    else if(geom_type == "ramp") 
    {
	amrex::Print() << " Ramp using slope!\n"; 
	int upDir; 
	int indepVar; 
	Real startPt; 
	Real slope; 
	pp.get("up_dir", upDir); 
	pp.get("indep_var", indepVar);
	pp.get("start_pt", startPt); 
	pp.get("ramp_slope", slope); 
	RealArray normal; 
	RealArray point; 
	for(int idir = 0; idir < SpaceDim; idir++)
	{
		normal[idir] = 0.; 
		point[idir] = 0.; 		
	}
	normal[upDir] = -1.0; 
	normal[indepVar] = slope; 
	
	point[upDir] = -startPt*slope;
	EB2::PlaneIF myplane(point, normal); 
	
	auto gshop = EB2::makeShop(myplane); 
	EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level); 
	 
    }
    else if(geom_type == "ramp_normal_point")
    {
	amrex::Print() << " Ramp Geometry using normal and point directly \n"; 
	RealArray normal; 
	RealArray point; 
	Vector<Real> pointvec; 
	Vector<Real> normalvec; 
	pp.getarr("ramp_normal", normalvec, 0, SpaceDim); 
	pp.getarr("ramp_point", pointvec, 0, SpaceDim); 
	for(int idir = 0; idir < SpaceDim; idir++)
	{
		point[idir] = pointvec[idir]; 
		normal[idir] = normalvec[idir]; 
	}
	EB2::PlaneIF myplane(point,normal); 
	auto gshop = EB2::makeShop(myplane); 
	EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level); 
    }
    else
    {
        EB2::Build(geom, max_coarsening_level, max_coarsening_level);
    }
}


