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

#include <AMReX_EB_LSCoreBase.H>

using namespace amrex;


/********************************************************************************
 * ALIAS for the otherwise loooooooooong template parameter names:              *
 *       WallsIF          : Implicit function for two planar walls              *
 *       CylinderIF       : Implicit for a polynomial (which has been           *
 *                          translated representing the cylinder walls          *
 *       CappedCylinderIF : Implicit for a polynomial (e.g. cylinder walls)     *
 *                          with planar caps on either end                      *
 *******************************************************************************/

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


void make_my_eb2(int n_lev, const BoxArray & grids, const DistributionMapping & dmap,
                 const Geometry & geom, LSFactory * level_set, LSCoreBase *& ls_core);


std::unique_ptr<CappedCylinderIF>
make_cylinder_eb2_geom(int dir, Real radius, Real length, const RealVect & translation,
                       int lev, const Geometry & geom, const DistributionMapping & dm,
                       LSFactory * level_set, LSCoreBase *& ls_core);


std::unique_ptr<CylinderIF>
make_poly_eb2_geom(int lev, int max_order, std::string field_prefix);
