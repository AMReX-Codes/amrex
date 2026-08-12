#ifndef AMREX_MARCHING_CUBES_H_
#define AMREX_MARCHING_CUBES_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>

#include <map>
#include <memory>
#include <string>

/**
 * \file AMReX_MarchingCubes.H
 *
 * GPU-enabled marching-cubes routines used to extract EB surfaces.
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
class MCFab
{
public:
    //! Per-cell triangle/interior ranges and resolved MC33 face decisions.
    BaseFab<int> m_cell_data;
    //! Exact normalized STL crossings on Cartesian edges, when available.
    Array<BaseFab<Real>,AMREX_SPACEDIM> m_edge_intersections;
    Triangle m_triangles;
    Vertex m_vertices;

    //! Allocate edge-intersection caches for the usable nodal SDF region.
    void defineEdgeIntersections (Box const& sdf_box);
};

//! Run marching cubes on signed-distance field \p sdf_fab.
void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab);

/**
 * \brief Build Cartesian-face aperture fractions and centroids from the MC
 *        signed-distance field.
 *
 * \p bx is cell-centered; all faces surrounding it are updated.
 *
 * \return Number of faces with unsupported contour topology.
 */
int build_face_fractions (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& fcx, FArrayBox& fcy, FArrayBox& fcz);

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
 * \return Number of single-valued cells whose triangulation produced invalid
 *         geometry.
 */
int build_cell_fractions (
    Box const& bx, Geometry const& geom, MCFab const& mc_fab,
    FArrayBox const& sdf,
    FArrayBox& apx, FArrayBox& apy, FArrayBox& apz,
    FArrayBox& vfrac, FArrayBox& vcent, FArrayBox& barea,
    FArrayBox& bcent, FArrayBox& bnorm);

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
 * \return Number of rejected faces across all three coordinate directions.
 */
int mark_faces_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    IArrayBox& rejected_x, IArrayBox& rejected_y,
    IArrayBox& rejected_z);

/**
 * \brief Mark cells that the single-VoF MC builder must remove.
 *
 * A cell is rejected if its MC boundary is invalid, its resolved fluid
 * topology has more than one fluid component, or its volume fraction is below
 * \p small_volfrac. Face-aperture rejection is handled separately.
 *
 * \return Counts of topology and small-cell rejections, respectively. The
 *         corresponding bits in \p rejected preserve those reasons.
 */
GpuArray<int,2> mark_cells_for_cleanup (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf,
    FArrayBox const& vfrac, Real small_volfrac, IArrayBox& rejected);

/**
 * \brief Apply the legacy EB repair transform to the MC working level set.
 *
 * Positive (fluid) nodes incident to a rejected face or cell are moved to
 * exact zero. The operation is monotone and leaves all other samples intact.
 *
 * \return Number of nodes changed from positive to zero.
 */
Long zero_nodes_for_cleanup (
    Box const& node_box, IArrayBox const& rejected_cells,
    IArrayBox const& rejected_x, IArrayBox const& rejected_y,
    IArrayBox const& rejected_z, FArrayBox& sdf);

/**
 * \brief Extrude the nodal level set straight outward from physical domain faces.
 *
 * An exterior node is assigned the value of the corresponding physical-face
 * node only when it is outside in exactly one non-periodic direction. This
 * preserves the feature's transverse footprint and does not extend diagonally
 * from domain edges or corners.
 *
 * \return Number of exterior nodal values changed.
 */
Long extend_domain_face_levelset (
    Box const& node_box, Box const& domain, GpuArray<int,3> const& is_periodic,
    FArrayBox& sdf);

/**
 * \brief Extrude exact transverse-edge crossings with the domain-face level set.
 *
 * Crossings on edges transverse to the extension normal are copied straight
 * outward from the physical face. Edges parallel to the normal are unchanged.
 */
void extend_domain_face_edge_intersections (
    Box const& domain, GpuArray<int,3> const& is_periodic, MCFab& mc_fab);

//! Write the converged MC triangles owned by \p grids to an STL file.
void write_stl (
    std::string const& filename,
    std::map<int,std::unique_ptr<MCFab>> const& mc_fabs,
    BoxArray const& grids);

}

#endif
