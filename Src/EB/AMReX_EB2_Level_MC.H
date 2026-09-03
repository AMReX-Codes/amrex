#ifndef AMREX_EB2_LEVEL_MC_H_
#define AMREX_EB2_LEVEL_MC_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_BaseFab.H>
#include <AMReX_EB2_Graph.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>

/**
 * \file AMReX_EB2_Level_MC.H
 *
 * Adapters that produce the two geometry-source inputs of the marching-cubes
 * EB builder (see AMReX_MarchingCubes.H, "Geometry-source contract"):
 *
 * - the nodal MC level set (positive in fluid), and
 * - the exact Cartesian-edge crossings cached in MC::MCFab.
 *
 * A geometry source that provides its own MC sampling (STLtools) is used as is.
 * Every other source, i.e. any EB2::GeometryShop implicit function, is served
 * by the generic adapters below, which reuse the source's `fillFab` and
 * `getIntercept` (Brent root finding on cut edges).
 *
 * These are implementation details of GShopLevel; they are not part of the
 * public EB2 interface.
 */

namespace amrex::EB2::detail {

//! Geometry sources with dedicated marching-cubes sampling (e.g. STLtools).
template <class G>
concept HasMarchingCubesSampling = requires (G const& g, MultiFab& mf, IntVect const& ng,
                                             Geometry const& geom, FArrayBox const& fab,
                                             MC::MCFab& mc_fab)
{
    g.fillMarchingCubesLevelSet(mf, ng, geom);
    g.fillMarchingCubesEdgeIntersections(fab, geom, mc_fab);
};

//! Execution target used by the generic marching-cubes adapters for \p gshop.
template <class G>
RunOn marching_cubes_run_on (G const& gshop) noexcept
{
    return (Gpu::inLaunchRegion() && gshop.isGPUable()) ? RunOn::Gpu : RunOn::Cpu;
}

/**
 * \brief Fill the nodal marching-cubes level set from geometry source \p gshop.
 *
 * On return every node of \p sdf (including ghost nodes) holds a value that is
 * positive in fluid and negative in the covered region.  Every node is a
 * deterministic function of the geometry alone, so shared and ghost nodes
 * agree across FABs without communication.  Implicit functions are sampled
 * without clamping to the domain because the builder applies its own
 * domain-face extension; as in the legacy generator, ghost nodes across a
 * periodic boundary sample the function at their unwrapped coordinates.
 */
template <class G>
void fill_marching_cubes_levelset (G const& gshop, MultiFab& sdf, Geometry const& geom)
{
    if constexpr (HasMarchingCubesSampling<G>) {
        gshop.fillMarchingCubesLevelSet(sdf, sdf.nGrowVect(), geom);
    } else {
        BL_PROFILE("EB2::fill_marching_cubes_levelset");
        RunOn const run_on = marching_cubes_run_on(gshop);

        // Roundoff can leave a node that lies on the surface at a tiny nonzero
        // value instead of exactly zero, which would create a sliver cut cell
        // that the builder cannot represent.  Snap such nodes to the exact-zero
        // ON state.  The threshold is relative to the local magnitude of the
        // implicit function, i.e. it bounds the crossing position on the
        // incident edges as a fraction of the cell size, so it is independent
        // of the function's units; it matches the sliver tolerance below which
        // the marching-cubes moments reject a cut.
#ifdef AMREX_USE_FLOAT
        constexpr Real snap_factor = 1.e-4_rt;
#else
        constexpr Real snap_factor = 1.e-11_rt;
#endif
        for (MFIter mfi(sdf); mfi.isValid(); ++mfi) {
            Box const nbx = sdf[mfi].box();
            // Sample one more ring of nodes so that the snap stencil is
            // complete on every node of the FAB, ghost nodes included.
            FArrayBox raw(amrex::grow(nbx,1), 1, The_Arena());
            Elixir raw_eli = raw.elixir();
            gshop.fillFab(raw, geom, run_on, raw.box());
            auto const f = raw.const_array();
            auto const a = sdf.array(mfi);
            ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const v = f(i,j,k);
                Real scale = std::abs(v);
                scale = amrex::max(scale, std::abs(f(i-1,j,k)), std::abs(f(i+1,j,k)));
                scale = amrex::max(scale, std::abs(f(i,j-1,k)), std::abs(f(i,j+1,k)));
                scale = amrex::max(scale, std::abs(f(i,j,k-1)), std::abs(f(i,j,k+1)));
                // EB2 implicit functions are positive inside the body; marching
                // cubes extracts the positive material, i.e. the fluid.
                a(i,j,k) = (std::abs(v) <= snap_factor*scale) ? 0.0_rt : -v;
            });
        }
        Gpu::streamSynchronize();
    }
}

/**
 * \brief Cache exact edge crossings for one FAB from geometry source \p gshop.
 *
 * \p sdf is the FAB's nodal marching-cubes level set (positive in fluid).  For
 * every edge whose end points differ in sign, the crossing position from the
 * low node, normalized by the cell size, is stored in
 * \p mc_fab.m_edge_intersections; other edges keep MC::invalid_edge_intersection.
 */
template <class G>
void fill_marching_cubes_edge_intersections (G const& gshop, FArrayBox const& sdf_fab,
                                             Geometry const& geom, MC::MCFab& mc_fab)
{
    if constexpr (HasMarchingCubesSampling<G>) {
        gshop.fillMarchingCubesEdgeIntersections(sdf_fab, geom, mc_fab);
    } else {
        BL_PROFILE("EB2::fill_marching_cubes_edge_intersections");
        RunOn const run_on = marching_cubes_run_on(gshop);
        auto const sdf = sdf_fab.const_array();
        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        Array<BaseFab<Type_t>,AMREX_SPACEDIM> type_fabs;
        Array<BaseFab<Real>,AMREX_SPACEDIM> inter_fabs;
        Array<Array4<Type_t const>,AMREX_SPACEDIM> type_arrs;
        Array<Array4<Real>,AMREX_SPACEDIM> inter_arrs;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Box const ebx = mc_fab.m_edge_intersections[idim].box();
            AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(amrex::surroundingNodes(ebx)));
            type_fabs[idim].resize(ebx, 1);
            inter_fabs[idim].resize(ebx, 1);
            IntVect const shift = IntVect::TheDimensionVector(idim);
            auto const type = type_fabs[idim].array();
            ParallelFor(ebx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const lo = sdf(i,j,k);
                Real const hi = sdf(i+shift[0],j+shift[1],k+shift[2]);
                // An exact-zero (ON) end point already is the crossing; the
                // builder places it at that node without an exact intercept.
                // Root finding on such an edge is not needed and could fail
                // when the zero came from snapping a roundoff-level value.
                // Away from snapped nodes the field is the negated implicit
                // function, so a strict sign change brackets the root for
                // getIntercept.
                bool const strict_sign_change = (lo > 0.0_rt && hi < 0.0_rt)
                                             || (lo < 0.0_rt && hi > 0.0_rt);
                type(i,j,k) = strict_sign_change ? Type::irregular : Type::regular;
            });
            type_arrs[idim] = type_fabs[idim].const_array();
            inter_arrs[idim] = inter_fabs[idim].array();
        }
        Gpu::streamSynchronize();

        // Physical crossing coordinates on irregular edges; NaN elsewhere.
        // GeometryShop ignores the level-set argument (it is only a hint for
        // STL ray casting), so the MC-convention values are passed unchanged.
        gshop.getIntercept(inter_arrs, type_arrs, sdf, geom, run_on, sdf_fab.box());

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Box const ebx = mc_fab.m_edge_intersections[idim].box();
            auto const type = type_fabs[idim].const_array();
            auto const inter = inter_fabs[idim].const_array();
            auto const crossing = mc_fab.m_edge_intersections[idim].array();
            Real const problo_d = problo[idim];
            Real const dx_d = dx[idim];
            ParallelFor(ebx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (type(i,j,k) == Type::irregular) {
                    Real const x = inter(i,j,k);
                    if (amrex::isnan(x)) {
                        crossing(i,j,k) = MC::invalid_edge_intersection;
                    } else {
                        int const iv[3] = {i,j,k};
                        Real const lo = problo_d + static_cast<Real>(iv[idim])*dx_d;
                        crossing(i,j,k) = amrex::Clamp((x-lo)/dx_d, 0.0_rt, 1.0_rt);
                    }
                }
            });
        }
        Gpu::streamSynchronize();
    }
}

} // namespace amrex::EB2::detail

#endif
