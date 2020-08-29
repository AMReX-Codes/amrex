#ifndef AMREX_EB_GEOMETRY_H_
#define AMREX_EB_GEOMETRY_H_


#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>



namespace amrex {
namespace geom {



RealVect facets_nearest_pt (const IntVect  & ind_pt, const IntVect  & ind_loop,
                            const RealVect & r_vec,  const RealVect & eb_normal,
                            const RealVect & eb_p0,  const RealVect & dx );


void closest_dist (Real & min_dist, bool & proj_valid,
                   const Vector<Real> & eb_data,
                   const RealVect & dx_eb, const RealVect & pos);

}
}


#endif
