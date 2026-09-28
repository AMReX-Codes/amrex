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
 * Adapters that produce the marching-cubes level set and edge crossings from
 * a geometry source.  Sources with their own marching-cubes sampling (STLtools)
 * are used directly; GeometryShop implicit functions use fillFab and
 * getIntercept.
 */

namespace amrex::EB2::detail {

void prefill_volume_fractions (Box const& bx, Array4<Real const> const& sdf,
                               Array4<Real> const& vfrac);
void negate_levelset (Box const& bx, Array4<Real> const& phi);

//! Geometry sources with dedicated marching-cubes sampling (e.g. STLtools).
template <class G>
concept HasMarchingCubesSampling = requires (G const& g, MultiFab& mf, IntVect const& ng,
                                             Geometry const& geom, FArrayBox const& fab,
                                             Array<Array4<Real>,AMREX_SPACEDIM> const& edges)
{
    g.fillMarchingCubesLevelSet(mf, ng, geom);
    g.fillMarchingCubesEdgeIntersections(fab, geom, edges);
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
 * All nodes of \p sdf, including ghost nodes, are set to a value that is
 * positive in fluid and negative in the covered region.  Each node depends
 * only on its position, so shared and ghost nodes agree across FABs without
 * communication.
 */
template <class G>
void fill_marching_cubes_levelset (G const& gshop, MultiFab& sdf, Geometry const& geom)
{
    if constexpr (HasMarchingCubesSampling<G>) {
        gshop.fillMarchingCubesLevelSet(sdf, sdf.nGrowVect(), geom);
    } else {
        BL_PROFILE("EB2::fill_marching_cubes_levelset");
        RunOn const run_on = marching_cubes_run_on(gshop);

        // Snap nodes within roundoff of the surface to exactly zero to avoid
        // sliver cut cells.  The threshold is relative to the neighboring
        // values, so it is independent of the implicit function's units.
#ifdef AMREX_USE_FLOAT
        constexpr Real snap_factor = 1.e-4_rt;
#else
        constexpr Real snap_factor = 1.e-11_rt;
#endif
        for (MFIter mfi(sdf); mfi.isValid(); ++mfi) {
            Box const nbx = sdf[mfi].box();
            // One extra ring for the snap stencil.
            FArrayBox raw(amrex::grow(nbx,1), 1,
                          run_on == RunOn::Gpu ? The_Async_Arena() : The_Arena());
            Elixir raw_eli;
            if (run_on == RunOn::Cpu) {
                raw_eli = raw.elixir();
            }
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
                // EB2 implicit functions are negative in fluid.
                a(i,j,k) = (std::abs(v) <= snap_factor*scale) ? 0.0_rt : -v;
            });
        }
    }
}

/**
 * \brief Cache exact edge crossings for one FAB from geometry source \p gshop.
 *
 * For every edge of \p sdf_fab whose end points differ in sign, the crossing
 * position from the low node, normalized by the cell size, is stored in
 * \p edges.  The caller must initialize \p edges to
 * MC::invalid_edge_intersection.
 */
template <class G>
void fill_marching_cubes_edge_intersections (G const& gshop, FArrayBox const& sdf_fab,
                                             Geometry const& geom,
                                             Array<Array4<Real>,AMREX_SPACEDIM> const& edges)
{
    if constexpr (HasMarchingCubesSampling<G>) {
        gshop.fillMarchingCubesEdgeIntersections(sdf_fab, geom, edges);
    } else {
        BL_PROFILE("EB2::fill_marching_cubes_edge_intersections");
        RunOn const run_on = marching_cubes_run_on(gshop);
        auto const sdf = sdf_fab.const_array();
        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        // Host-only getIntercept uses The_Arena() to decide whether staging
        // is needed, so its buffers must retain that arena.
        Arena* arena = (run_on == RunOn::Gpu) ? The_Async_Arena() : The_Arena();
        Array<BaseFab<Type_t>,AMREX_SPACEDIM> type_fabs;
        Array<BaseFab<Real>,AMREX_SPACEDIM> inter_fabs;
        Array<Elixir,AMREX_SPACEDIM> type_elixirs;
        Array<Elixir,AMREX_SPACEDIM> inter_elixirs;
        Array<Array4<Type_t const>,AMREX_SPACEDIM> type_arrs;
        Array<Array4<Real>,AMREX_SPACEDIM> inter_arrs;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Box const ebx = amrex::enclosedCells(sdf_fab.box(), idim);
            AMREX_ALWAYS_ASSERT(Box(edges[idim]).setType(ebx.ixType()).contains(ebx));
            AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(amrex::surroundingNodes(ebx)));
            type_fabs[idim].resize(ebx, 1, arena);
            inter_fabs[idim].resize(ebx, 1, arena);
            if (run_on == RunOn::Cpu) {
                // The kernels below may still run on the device.
                type_elixirs[idim] = type_fabs[idim].elixir();
                inter_elixirs[idim] = inter_fabs[idim].elixir();
            }
            IntVect const shift = IntVect::TheDimensionVector(idim);
            auto const type = type_fabs[idim].array();
            ParallelFor(ebx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const lo = sdf(i,j,k);
                Real const hi = sdf(i+shift[0],j+shift[1],k+shift[2]);
                // A zero end point is already the crossing, and root finding
                // could fail if it came from snapping.
                bool const strict_sign_change = (lo > 0.0_rt && hi < 0.0_rt)
                                             || (lo < 0.0_rt && hi > 0.0_rt);
                type(i,j,k) = strict_sign_change ? Type::irregular : Type::regular;
            });
            type_arrs[idim] = type_fabs[idim].const_array();
            inter_arrs[idim] = inter_fabs[idim].array();
        }
        // CPU root finding can read managed edge types directly. GPU root
        // finding and explicit device-to-host copies are stream ordered.
        if (run_on == RunOn::Cpu && The_Arena()->isHostAccessible()) {
            Gpu::streamSynchronize();
        }

        // GeometryShop does not use the level-set argument, so its sign
        // convention does not matter here.
        gshop.getIntercept(inter_arrs, type_arrs, sdf, geom, run_on, sdf_fab.box());

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Box const ebx = amrex::enclosedCells(sdf_fab.box(), idim);
            auto const type = type_fabs[idim].const_array();
            auto const inter = inter_fabs[idim].const_array();
            auto const crossing = edges[idim];
            Real const problo_d = problo[idim];
            Real const dx_d = dx[idim];
            ParallelFor(ebx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (type(i,j,k) == Type::irregular) {
                    Real const x = inter(i,j,k);
                    if (x == no_intercept) {
                        crossing(i,j,k) = MC::invalid_edge_intersection;
                    } else {
                        int const iv[3] = {i,j,k};
                        Real const lo = problo_d + static_cast<Real>(iv[idim])*dx_d;
                        crossing(i,j,k) = amrex::Clamp((x-lo)/dx_d, 0.0_rt, 1.0_rt);
                    }
                }
            });
        }
    }
}

} // namespace amrex::EB2::detail

#endif
