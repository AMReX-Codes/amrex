#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>

#include <AMReX_GeometryShop.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBTower.H>
#include <AMReX_EBFArrayBox.H>

#include <AMReX_BaseIF.H>
#include <AMReX_EB_levelset.H>

std::unique_ptr<amrex::BaseIF>
make_cylinder_geom(int dir, amrex::Real radius, amrex::Real length, const amrex::RealVect & translation,
                   int lev, const amrex::Geometry & geom, const amrex::DistributionMapping & dm, amrex::LSFactory * level_set);

std::unique_ptr<amrex::BaseIF>
make_poly_geom(int lev, int max_order, std::string field_prefix);

void
make_my_eb(int lev, const amrex::BoxArray& grids, const amrex::DistributionMapping& dmap,
           const amrex::Geometry& geom, amrex::LSFactory * level_set);
