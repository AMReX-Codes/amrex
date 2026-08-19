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
 * positive in fluid and negative in the covered region, and shared nodes are
 * identical across FABs.  Implicit functions are sampled without clamping to
 * the domain because the builder applies its own domain-face extension.
 */
template <class G>
void fill_marching_cubes_levelset (G const& gshop, MultiFab& sdf, Geometry const& geom)
{
    if constexpr (HasMarchingCubesSampling<G>) {
        gshop.fillMarchingCubesLevelSet(sdf, sdf.nGrowVect(), geom);
    } else {
        BL_PROFILE("EB2::fill_marching_cubes_levelset");
        RunOn const run_on = marching_cubes_run_on(gshop);
        for (MFIter mfi(sdf); mfi.isValid(); ++mfi) {
            auto& fab = sdf[mfi];
            gshop.fillFab(fab, geom, run_on, fab.box());
            // EB2 implicit functions are positive inside the body; marching
            // cubes extracts the positive material, i.e. the fluid.
            auto const a = fab.array();
            ParallelFor(fab.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i,j,k) = -a(i,j,k);
            });
        }
        Gpu::streamSynchronize();
        sdf.FillBoundary(geom.periodicity());

        // Roundoff can leave a node that lies on the surface at a tiny nonzero
        // value instead of exactly zero, which would create a sliver cut cell
        // that the builder cannot represent.  Snap such nodes to the exact-zero
        // ON state.  The threshold is relative to the local magnitude of the
        // implicit function, i.e. it bounds the crossing position on the
        // incident edges as a fraction of the cell size, so it is independent
        // of the function's units; it matches the sliver tolerance below which
        // the marching-cubes moments reject a cut.  The decision reads the
        // unsnapped values from a copy of the FAB, so it does not depend on
        // the order in which nodes are visited; valid nodes are updated and
        // OverrideSync then makes shared nodes identical across FABs.
#ifdef AMREX_USE_FLOAT
        constexpr Real snap_factor = 1.e-4_rt;
#else
        constexpr Real snap_factor = 1.e-11_rt;
#endif
        for (MFIter mfi(sdf); mfi.isValid(); ++mfi) {
            Box const vbx = mfi.validbox();
            AMREX_ALWAYS_ASSERT(sdf[mfi].box().contains(amrex::grow(vbx,1)));
            FArrayBox unsnapped(amrex::grow(vbx,1), 1, The_Arena());
            Elixir unsnapped_eli = unsnapped.elixir();
            unsnapped.copy<RunOn::Device>(sdf[mfi], amrex::grow(vbx,1));
            auto const src = unsnapped.const_array();
            auto const a = sdf.array(mfi);
            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const v = src(i,j,k);
                Real scale = std::abs(v);
                scale = amrex::max(scale, std::abs(src(i-1,j,k)), std::abs(src(i+1,j,k)));
                scale = amrex::max(scale, std::abs(src(i,j-1,k)), std::abs(src(i,j+1,k)));
                scale = amrex::max(scale, std::abs(src(i,j,k-1)), std::abs(src(i,j,k+1)));
                if (std::abs(v) <= snap_factor*scale) {
                    a(i,j,k) = 0.0_rt;
                }
            });
        }
        Gpu::streamSynchronize();
        sdf.OverrideSync(geom.periodicity());
        sdf.FillBoundary(geom.periodicity());
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

        // The nodal field has been synchronized across FABs and periodic
        // images, so its ghost nodes need not agree with the implicit function
        // evaluated at their physical coordinates (a non-periodic function in
        // a periodic domain).  Root finding brackets with the function itself,
        // so decide which edges to bracket from a raw sample of the function
        // (EB2 convention: positive in the body).
        FArrayBox raw(sdf_fab.box(), 1, The_Arena());
        Elixir raw_eli = raw.elixir();
        gshop.fillFab(raw, geom, run_on, raw.box());
        auto const raw_if = raw.const_array();

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
                bool const strict_sign_change = (lo > 0.0_rt && hi < 0.0_rt)
                                             || (lo < 0.0_rt && hi > 0.0_rt);
                Real const raw_lo = raw_if(i,j,k);
                Real const raw_hi = raw_if(i+shift[0],j+shift[1],k+shift[2]);
                bool const bracketed = (raw_lo > 0.0_rt && raw_hi < 0.0_rt)
                                    || (raw_lo < 0.0_rt && raw_hi > 0.0_rt);
                // Edges where the synchronized field changes sign but the
                // function does not keep the linear fallback.
                type(i,j,k) = (strict_sign_change && bracketed) ? Type::irregular
                                                                : Type::regular;
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
