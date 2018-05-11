#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_EB2_GeometryShop.H>

namespace amrex { namespace EB2 {

template <typename G>
class ISLevel
{
public:

    ISLevel (G const& gshop, Geometry const& geom, int max_grid_size);

private:

    Geometry m_geom;    
    BoxArray m_grids;
};

template <typename G>
ISLevel<G>::ISLevel (G const& gshop, Geometry const& geom, int max_grid_size)
    : m_geom(geom)
{
    
}


}}

#endif
