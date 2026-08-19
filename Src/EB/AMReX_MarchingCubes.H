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
 * Geometry-source contract
 * ------------------------
 * The marching-cubes EB builder is independent of how the geometry is
 * described (STL triangulation, analytic implicit function, parser
 * expression, ...).  A geometry source feeds it through exactly two inputs:
 *
 * 1. A nodal scalar field (the "MC level set").  It uses the marching-cubes
 *    sign convention: values greater than zero are fluid, values less than
 *    zero are covered, and an exact zero marks a node lying on the boundary
 *    (it takes the covered sign for the two-state MC33 lookup, and edge
 *    crossings are placed exactly at that node).  Only the sign and local
 *    continuity are required; magnitudes are used for the linear crossing
 *    fallback and the MC33 asymptotic decider.  Note that this is the
 *    opposite of the EB2 implicit-function convention (negative in fluid),
 *    which is why the builder negates the field on ingest and again when it
 *    exports the public level set.
 *
 * 2. Optional exact crossings on Cartesian edges, stored in
 *    MCFab::m_edge_intersections as the normalized position in [0,1] measured
 *    from the low node.  Edges without an exact crossing keep
 *    invalid_edge_intersection and fall back to linear interpolation of the
 *    nodal field.
 *
 * See AMReX_EB2_Level_MC.H for the adapters that produce these inputs from
 * EB2::GeometryShop implicit functions and from STLtools.
 */

namespace amrex::MC {

//! Sentinel used when an original-STL edge crossing is unavailable.
inline constexpr Real invalid_edge_intersection = -1.0_rt;

enum RejectionReason : int {
    invalid_topology = 1,
    small_volume = 2
};

/**
 * Components stored in MCFab::m_cell_data.  Besides compact triangle and
 * interior-vertex ranges, these retain the MC33 face decisions consumed by
 * downstream face construction.
 */
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
 * \brief Slots of the per-FAB counter block that the entry points below
 *        accumulate into.
 *
 * The caller allocates num_fab_counters ints per FAB in device-accessible
 * memory (Gpu::Buffer or Gpu::DeviceVector), zeroes them once per repair
 * pass, passes each FAB's block to the entry points, and copies all blocks to
 * the host once per pass.  The entry points themselves never synchronize with
 * the host to report their counts.
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
    //! Per-cell triangle/interior ranges and resolved MC33 face decisions.
    BaseFab<int> m_cell_data;
    //! Exact normalized STL crossings on Cartesian edges, when available.
    Array<BaseFab<Real>,AMREX_SPACEDIM> m_edge_intersections;
    Triangle m_triangles;
    Vertex m_vertices;

    //! Allocate edge-intersection caches for the usable nodal SDF region.
    void defineEdgeIntersections (Box const& sdf_box);
};

/**
 * \brief Run marching cubes on signed-distance field \p sdf_fab.
 *
 * \p counters is the FAB's counter block (see FabCounter); a nonzero
 * counter_invalid_triangles afterwards indicates a lookup-table invariant
 * violation and the caller must abort.
 */
void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab, int* counters);

/**
 * \brief Build Cartesian-face aperture fractions and centroids from the MC
 *        signed-distance field.
 *
 * \p bx is cell-centered; all faces surrounding it are updated.  Faces whose
 * fluid polygon is degenerate (for example a four-crossing face whose
 * crossings all sit on repaired zero nodes) receive zero aperture and are
 * marked with 1 in \p rejected_x / \p rejected_y / \p rejected_z where those
 * arrays contain the face, so the nodal repair pass covers them.
 *
 * Adds the number of faces whose two cells resolved the MC33 ambiguity
 * differently (an invariant violation) to counter_face_decision_errors and
 * the number of degenerate faces marked for repair to counter_degenerate_faces
 * of \p counters.
 */
void build_face_fractions (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& fcx, FArrayBox& fcy, FArrayBox& fcz,
    IArrayBox& rejected_x, IArrayBox& rejected_y, IArrayBox& rejected_z,
    int* counters);

/**
 * \brief Build EB edge-centroid data from the MC signed-distance field.
 *
 * Regular and covered edges receive the standard EB2 sentinels 1 and -1.
 * Cut-edge values are the centroid of the fluid part of the unit edge.
 */
void build_edge_centroids (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox& ecx, FArrayBox& ecy, FArrayBox& ecz);

/**
 * \brief Replace single-valued cut-cell geometry with moments computed from
 *        the marching-cubes triangles.
 *
 * Coordinate-face apertures previously built from the MC signed-distance
 * field and the cell-owned embedded-boundary triangles form the closed
 * boundary used to compute volume and first moments.
 *
 * Single-valued cells whose triangulation produced invalid geometry are
 * counted in counter_closure_errors, counter_volume_errors,
 * counter_centroid_errors and counter_area_vector_errors of \p counters.
 */
void build_cell_fractions (
    Box const& bx, Geometry const& geom, MCFab const& mc_fab,
    FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& vfrac, FArrayBox& vcent, FArrayBox& barea,
    FArrayBox& bcent, FArrayBox& bnorm, int* counters);

/**
 * \brief Classify cells and build all EBCellFlag connectivity from MC SDF and
 *        face apertures.
 *
 * Cell types are derived from the converged volume and face apertures,
 * including coincident nodal surfaces. MC ownership is checked against the
 * repaired SDF before coordinate connectivity is populated.
 *
 * \return Number of inconsistent cell classifications.
 */
int build_cell_topology (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    EBCellFlagFab& cellflag, FArrayBox const& vfrac,
    FArrayBox const& apx, FArrayBox const& apy, FArrayBox const& apz);

/**
 * \brief Mark Cartesian faces with more than one resolved fluid aperture.
 *
 * Marks are OR-ed into the arrays, so faces already rejected by
 * build_face_fractions are kept.  The number of newly rejected faces across
 * all three coordinate directions is added to counter_face_rejections of
 * \p counters.
 */
void mark_faces_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    IArrayBox& rejected_x, IArrayBox& rejected_y,
    IArrayBox& rejected_z, int* counters);

/**
 * \brief Mark cells that the single-VoF MC builder must remove.
 *
 * A cell is rejected if its MC boundary is invalid, if its fluid corners form
 * more than one group, or if its volume fraction is below \p small_volfrac.
 * Corner groups are joined along cube edges and across ambiguous faces whose
 * MC33 decision is "fluid connected" (the same decision the surface and the
 * face apertures use), but never through the cell interior: MC33 tunnel
 * tilings (4.1.2, 10.1.2, ...) and the split tilings 7.3/10.2/12.2/13.x are
 * always rejected, so an exported single-valued cell holds exactly one
 * face-connected fluid region.  Face-aperture rejection is handled separately
 * by mark_faces_for_cleanup and is consistent with this rule (a face whose
 * decision is "not connected" is a rejected face and also splits the corner
 * groups of both cells).
 *
 * Topology and small-cell rejections are added to counter_topology_rejections
 * and counter_small_cell_rejections of \p counters.  \p rejected records the
 * reason of each rejected cell as one of the RejectionReason values (topology
 * takes precedence) and 0 otherwise.
 */
void mark_cells_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox const& vfrac, Real small_volfrac, IArrayBox& rejected,
    int* counters);

/**
 * \brief Apply the legacy EB repair transform to the MC working level set.
 *
 * Positive (fluid) nodes incident to a rejected face or cell are moved to
 * exact zero. The operation is monotone and leaves all other samples intact.
 * The number of nodes changed from positive to zero is added to
 * counter_changed_nodes of \p counters.
 */
void zero_nodes_for_cleanup (
    Box const& node_box, IArrayBox const& rejected_cells,
    IArrayBox const& rejected_x, IArrayBox const& rejected_y,
    IArrayBox const& rejected_z, FArrayBox& sdf, int* counters);

/**
 * \brief Extrude the nodal level set straight outward from physical domain faces.
 *
 * An exterior node is assigned the value of the corresponding physical-face
 * node only when it is outside in exactly one non-periodic direction. This
 * preserves the feature's transverse footprint and does not extend diagonally
 * from domain edges or corners.  Exterior nodes that the nodal repair has
 * moved to exact zero are left at zero, so the repair stays monotone.  The
 * number of exterior nodal values changed is added to counter_extended_nodes
 * of \p counters.
 */
void extend_domain_face_levelset (
    Box const& node_box, Box const& domain, GpuArray<int,3> const& is_periodic,
    FArrayBox& sdf, int* counters);

/**
 * \brief Extrude exact transverse-edge crossings with the domain-face level set.
 *
 * Crossings on edges transverse to the extension normal are copied straight
 * outward from the physical face. Edges parallel to the normal are unchanged.
 */
void extend_domain_face_edge_intersections (
    Box const& domain, GpuArray<int,3> const& is_periodic, MCFab& mc_fab);

//! Write the converged MC triangles owned by the valid boxes of \p mc_fabs to an STL file.
void write_stl (std::string const& filename, LayoutData<MCFab> const& mc_fabs);

}

#endif
