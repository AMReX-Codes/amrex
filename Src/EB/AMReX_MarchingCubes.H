#ifndef AMREX_MARCHING_CUBES_H_
#define AMREX_MARCHING_CUBES_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>

#include <string>

/**
 * \file AMReX_MarchingCubes.H
 *
 * GPU-enabled marching-cubes routines used to extract EB surfaces.
 *
 * A geometry source provides two inputs to the marching-cubes EB builder:
 *
 * 1. A nodal level set that is positive in fluid and negative in the covered
 *    region.  An exact zero lies on the boundary and counts as covered.  This
 *    is the opposite of the EB2 implicit-function convention.
 *
 * 2. Optional exact crossings on Cartesian edges, stored as the position in
 *    [0,1] from the low node.  Edges without one hold
 *    invalid_edge_intersection and use linear interpolation of the level set.
 *
 * See AMReX_EB2_Level_MC.H for the adapters that produce these inputs.
 */

namespace amrex::MC {

//! Sentinel for an edge without an exact crossing.
inline constexpr Real invalid_edge_intersection = -1.0_rt;

enum RejectionReason : int {
    invalid_topology = 1,
    small_volume = 2
};

//! Components stored in MCFab::m_cell_data.
enum CellDataComponent : int {
    triangle_count = 0,
    triangle_offset,
    interior_vertex_count,
    interior_vertex_offset,
    face_decision_valid_mask,
    face_fluid_connected_mask,
    num_cell_data_components
};

/**
 * \brief Slots of the per-FAB counter block that the functions below
 *        accumulate into.
 *
 * The caller provides num_fab_counters ints per FAB in device-accessible
 * memory and copies them to the host when needed; the functions do not
 * synchronize with the host.
 */
enum FabCounter : int {
    counter_invalid_triangles = 0,   //!< marching_cubes: triangle without a crossing vertex
    counter_face_decision_errors,    //!< build_face_fractions
    counter_degenerate_faces,        //!< build_face_fractions
    counter_closure_errors,          //!< build_cell_fractions
    counter_volume_errors,           //!< build_cell_fractions
    counter_centroid_errors,         //!< build_cell_fractions
    counter_area_vector_errors,      //!< build_cell_fractions
    counter_face_rejections,         //!< mark_faces_for_cleanup
    counter_topology_rejections,     //!< mark_cells_for_cleanup
    counter_small_cell_rejections,   //!< mark_cells_for_cleanup
    counter_changed_nodes,           //!< zero_nodes_for_cleanup
    counter_extended_nodes,          //!< extend_domain_face_levelset
    num_fab_counters
};

//! Initialize internal lookup tables and device buffers for marching cubes.
void Initialize ();
//! Destroy any static storage allocated by Initialize().
void Finalize ();

//! Vertex buffer storing coordinates for extracted surfaces.
struct Vertex
{
    Gpu::DeviceVector<Real> x, y, z;

    //! Resize the vertex buffers to hold \p n entries.
    void resize (int n);

    //! Raw device pointers used when launching kernels.
    GpuArray<Real*,3> dataPtrs () { return {x.data(), y.data(), z.data()}; }
};

//! Triangle index buffer referencing Vertex entries.
struct Triangle
{
    Gpu::DeviceVector<int> v1, v2, v3; // vertex index

    //! Resize the triangle buffer to hold \p n faces.
    void resize (int n);

    //! Raw device pointers for kernel launches.
    GpuArray<int*,3> dataPtrs () { return {v1.data(), v2.data(), v3.data()}; }
};

//! Minimal container holding marching-cubes output for one FAB.
struct MCFab
{
    BaseFab<int> m_cell_data;
    Triangle m_triangles;
    Vertex m_vertices;
};

/**
 * \brief Run marching cubes on the nodal level set \p sdf_fab.
 *
 * \p counters is the FAB's counter block (see FabCounter).  \p edgecent holds
 * EB edge centroids; empty views use linear interpolation of \p sdf_fab.
 */
void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab, int* counters,
                     Array<Array4<Real const>,AMREX_SPACEDIM> const& edgecent = {});

/**
 * \brief Build face apertures and face centroids for all faces of the
 *        cell-centered box \p bx.
 *
 * Faces with a degenerate fluid polygon get zero aperture and are marked in
 * \p rejected_x, \p rejected_y and \p rejected_z for nodal repair.
 * \p edgecent must match the centroids passed to marching_cubes.
 */
void build_face_fractions (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& fcx, FArrayBox& fcy, FArrayBox& fcz,
    IArrayBox& rejected_x, IArrayBox& rejected_y, IArrayBox& rejected_z,
    int* counters, Array<Array4<Real const>,AMREX_SPACEDIM> const& edgecent = {});

/**
 * \brief Convert the edge arrays \p edges to EB edge centroids.
 *
 * Regular and covered edges get 1 and -1.  On the first pass \p edges holds
 * exact crossings (or invalid_edge_intersection); later passes update only
 * the edges touching repaired nodes.
 */
void build_edge_centroids (
    FArrayBox const& sdf, Array<Array4<Real>,AMREX_SPACEDIM> const& edges, bool first_pass);

/**
 * \brief Compute volume fraction, centroid and boundary moments of cut cells
 *        from the face apertures and the marching-cubes triangles.
 *
 * Requires cubic cells.
 */
void build_cell_fractions (
    Box const& bx, Geometry const& geom, MCFab const& mc_fab,
    FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& vfrac, FArrayBox& vcent, FArrayBox& barea,
    FArrayBox& bcent, FArrayBox& bnorm, int* counters);

/**
 * \brief Set cell types and connectivity in \p cellflag from the final
 *        volume and face fractions.
 *
 * \return Number of inconsistent cell classifications.
 */
int build_cell_topology (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    EBCellFlagFab& cellflag, FArrayBox const& vfrac,
    FArrayBox const& apx, FArrayBox const& apy, FArrayBox const& apz);

/**
 * \brief Mark faces with more than one fluid aperture.
 *
 * Marks are OR-ed into the arrays, so faces already marked by
 * build_face_fractions stay marked.
 */
void mark_faces_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    IArrayBox& rejected_x, IArrayBox& rejected_y,
    IArrayBox& rejected_z, int* counters);

/**
 * \brief Mark cells that cannot be represented as single-valued cut cells.
 *
 * A cell is rejected if its surface is invalid, if its fluid corners do not
 * form a single face-connected region, or if its volume fraction is below
 * \p small_volfrac.  \p rejected is set to a RejectionReason for rejected
 * cells and 0 otherwise.
 */
void mark_cells_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox const& vfrac, Real small_volfrac, IArrayBox& rejected,
    int* counters);

/**
 * \brief Set fluid nodes of \p node_box that touch a rejected cell or face
 *        to zero.
 *
 * Marks outside the rejection arrays are ignored, so ghost nodes can be
 * repaired from the marks held locally.
 */
void zero_nodes_for_cleanup (
    Box const& node_box, IArrayBox const& rejected_cells,
    IArrayBox const& rejected_x, IArrayBox const& rejected_y,
    IArrayBox const& rejected_z, FArrayBox& sdf, int* counters);

/**
 * \brief Extend the nodal level set outward from non-periodic domain faces.
 *
 * Only nodes outside the domain in exactly one direction are set; nodes
 * already zeroed by the repair are left at zero.
 */
void extend_domain_face_levelset (
    Box const& node_box, Box const& domain, GpuArray<int,3> const& is_periodic,
    FArrayBox& sdf, int* counters);

//! Extend exact crossings on edges transverse to non-periodic domain faces.
void extend_domain_face_edge_intersections (
    Box const& node_box, Box const& domain, GpuArray<int,3> const& is_periodic,
    Array<Array4<Real>,AMREX_SPACEDIM> const& edges);

//! Write the triangles in the valid boxes of \p mc_fabs to an STL file.
void write_stl (std::string const& filename, LayoutData<MCFab> const& mc_fabs);

}

#endif
