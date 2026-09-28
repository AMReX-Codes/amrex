#include <AMReX_EB2_Level.H>

#include <limits>

namespace amrex::EB2 {

namespace detail {

// These kernels are free functions because CUDA does not allow extended
// device lambdas in protected member functions.

// Cut cells get -1, which build_cell_fractions overwrites.
void prefill_volume_fractions (Box const& bx, Array4<Real const> const& sdf,
                               Array4<Real> const& vfrac)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int nfluid = 0;
        nfluid += sdf(i  ,j  ,k  ) > 0.0_rt;
        nfluid += sdf(i+1,j  ,k  ) > 0.0_rt;
        nfluid += sdf(i  ,j+1,k  ) > 0.0_rt;
        nfluid += sdf(i+1,j+1,k  ) > 0.0_rt;
        nfluid += sdf(i  ,j  ,k+1) > 0.0_rt;
        nfluid += sdf(i+1,j  ,k+1) > 0.0_rt;
        nfluid += sdf(i  ,j+1,k+1) > 0.0_rt;
        nfluid += sdf(i+1,j+1,k+1) > 0.0_rt;
        vfrac(i,j,k) = nfluid == 0 ? 0.0_rt : (nfluid == 8 ? 1.0_rt : -1.0_rt);
    });
}

void negate_levelset (Box const& bx, Array4<Real> const& phi)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        phi(i,j,k) = -phi(i,j,k);
    });
}

} // namespace detail

void
Level::assert_marching_cubes_supported (Geometry const& geom)
{
    auto const cell_size = geom.CellSizeArray();
    Real const max_cell_size = amrex::max(cell_size[0], amrex::max(cell_size[1], cell_size[2]));
    Real const cubic_tolerance = 16.0_rt * std::numeric_limits<Real>::epsilon() * max_cell_size;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.Coord() == CoordSys::cartesian &&
            std::abs(cell_size[0] - cell_size[1]) <= cubic_tolerance &&
            std::abs(cell_size[0] - cell_size[2]) <= cubic_tolerance,
        "Marching-cubes EB construction requires a 3D Cartesian grid with dx == "
        "dy == dz");
}

} // namespace amrex::EB2
