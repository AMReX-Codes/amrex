#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EBFabFactory.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

using namespace amrex;

using WallsIF = typename
    EB2::TranslationIF<
        EB2::UnionIF<
            EB2::PlaneIF,
            EB2::PlaneIF
        >
    >;

using CylinderIF = typename
    EB2::TranslationIF<
        EB2::PolynomialIF
    >;

using CappedCylinderIF = typename
    EB2::UnionIF<
        WallsIF,
        CylinderIF
    >;

void
make_my_eb2(int lev, const BoxArray & grids, const DistributionMapping & dmap,
            const Geometry & geom, LSFactory * level_set);


std::unique_ptr<CappedCylinderIF>
make_cylinder_eb2_geom(int dir, Real radius, Real length, const RealVect & translation,
                   int lev, const Geometry & geom, const DistributionMapping & dm,
                   LSFactory * level_set);

std::unique_ptr<CylinderIF>
make_poly_eb2_geom(int lev, int max_order, std::string field_prefix);
