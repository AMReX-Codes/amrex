#ifndef AMREX_EB_STL_UTILS_H_
#define AMREX_EB_STL_UTILS_H_

#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Dim3.H>
#include <AMReX_EB2_Graph.H>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>

/**
 * \file AMReX_EB_STL_utils.H
 *
 * Routines for sampling STL triangulations into EB level sets.
 */

namespace amrex
{

namespace MC {
struct MCFab;
}

//! Utility class that loads STL meshes and samples them onto AMReX grids.
class STLtools
{
public:

    /// \cond DOXYGEN_IGNORE
    struct Triangle {
        XDim3 v1, v2, v3;

        [[nodiscard]] Real cent (int d) const
        {
            static_assert(sizeof(XDim3) == sizeof(Real)*3);
            return Real(1./3.)*((&v1.x)[d] + (&v2.x)[d] + (&v3.x)[d]); // NOLINT(clang-analyzer-security.ArrayBound)
        }

        [[nodiscard]] std::pair<Real,Real> minmax (int d) const
        {
            static_assert(sizeof(XDim3) == sizeof(Real)*3);
            return std::minmax({(&v1.x)[d], (&v2.x)[d], (&v3.x)[d]}); // NOLINT(clang-analyzer-security.ArrayBound)
        }
    };
    /// \endcond

    /// \cond DOXYGEN_IGNORE
    template <int M, int N>
    struct BVHNodeT
    {
        RealBox boundingbox{AMREX_D_DECL(std::numeric_limits<Real>::max(),
                                         std::numeric_limits<Real>::max(),
                                         std::numeric_limits<Real>::max()),
                            AMREX_D_DECL(std::numeric_limits<Real>::lowest(),
                                         std::numeric_limits<Real>::lowest(),
                                         std::numeric_limits<Real>::lowest())};
        STLtools::Triangle triangles[M];
        XDim3 trinorm[M];
        int children[N];
        std::int8_t ntriangles = 0;
        std::int8_t nchildren = 0;
    };
    /// \endcond

    static constexpr int m_bvh_max_size = 4; // max # of triangles in a leaf node
    static constexpr int m_bvh_max_splits = 4; // max # of children
    static constexpr int m_bvh_max_stack_size = 12; // max depth of the tree

    using Node = BVHNodeT<m_bvh_max_size,m_bvh_max_splits>;

    static constexpr int allregular = -1;
    static constexpr int mixedcells = 0;
    static constexpr int allcovered = 1;

    //! Enable/disable BVH acceleration when sampling.
    void setBVHOptimization (bool flag) { m_bvh_optimization = flag; }

    //! Select the strict STL contract required by the marching-cubes builder.
    void setUseMarchingCubes (bool flag) noexcept { m_use_marching_cubes = flag; }
    [[nodiscard]] bool usesMarchingCubes () const noexcept {
        return m_use_marching_cubes;
    }

    /**
     * \brief Read an STL file, apply scaling/translation, and store triangles.
     *
     * \param fname           Path to the STL mesh.
     * \param scale           Uniform scaling factor applied to vertices.
     * \param center          Translation applied to recenter the mesh.
     * \param reverse_normal  Non-zero to flip normals.
     */
    void read_stl_file (std::string const& fname, Real scale, Array<Real,3> const& center,
                        int reverse_normal);

    /**
     * \brief Fill a level-set MultiFab sampled from the STL mesh.
     */
    void fill (MultiFab& mf, IntVect const& nghost, Geometry const& geom,
               Real outside_value = -1._rt, Real inside_value = 1._rt) const;

    //! Return whether the EB is regular/covered/mixed inside \p box.
    [[nodiscard]] int getBoxType (Box const& box, Geometry const& geom, RunOn) const;

    static constexpr bool isGPUable () noexcept { return true; }

    //! Fill a BaseFab with level-set data (host/device aware).
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn,
                  Box const& bounding_box) const;

    //! Compute intercepts and face types for EB geometry based on the STL mesh.
    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const& lst, Geometry const& geom, RunOn,
                       Box const& bounding_box) const;

    //! Update intercepts.
    static void updateIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                          Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                          Array4<Real const> const& lst, Geometry const& geom) ;

    //! Fill \p mf with signed-distance values sampled from the STL mesh.
    void fillSignedDistance (MultiFab& mf, IntVect const& nghost, Geometry const& geom) const;

    /**
     * \brief Fill the nodal scalar field used by the marching-cubes builder.
     *
     * Exact STL distances are evaluated only in the two-node band needed by
     * cut-cell interpolation and vertex-normal stencils.  Away from that
     * band, the field retains unit-magnitude inside/outside values.
     */
    void fillMarchingCubesLevelSet (MultiFab& mf, IntVect const& nghost,
                                    Geometry const& geom) const;

    /**
     * \brief Cache exact normalized Cartesian-edge crossings for MC geometry.
     *
     * The cache is evaluated against the original triangles rather than the
     * sampled signed-distance field. It remains valid through the monotone
     * nodal cleanup used by marching cubes.
     */
    void fillMarchingCubesEdgeIntersections (
        FArrayBox const& levelset, Geometry const& geom,
        MC::MCFab& mc_fab) const;

    //! Upload triangle data to device memory (exposed for CUDA workflows).
    void prepare (Gpu::PinnedVector<Triangle> a_tri_pts);  // public for cuda

private:

    bool m_bvh_optimization = true;
    bool m_use_marching_cubes = false;

    Gpu::DeviceVector<Triangle> m_tri_pts_d;
    Gpu::DeviceVector<XDim3> m_tri_normals_d;
    Gpu::DeviceVector<Node> m_bvh_nodes;

    int m_num_tri=0;

    XDim3 m_ptmin;  // All triangles are inside the bounding box defined by
    XDim3 m_ptmax;  //     m_ptmin and m_ptmax.
    XDim3 m_ptref;  // The reference point is slightly outside the bounding box.
    bool m_boundry_is_outside; // Is the bounding box boundary outside or inside the object?

    void read_ascii_stl_file (std::string const& fname, Real scale,
                              Array<Real,3> const& center, int reverse_normal,
                              Gpu::PinnedVector<Triangle>& a_tri_pts);
    void read_binary_stl_file (std::string const& fname, Real scale,
                               Array<Real,3> const& center, int reverse_normal,
                               Gpu::PinnedVector<Triangle>& a_tri_pts);

    static void build_bvh (Triangle* begin, Triangle * end, Gpu::PinnedVector<Node>& bvh_nodes);
    static void bvh_size (int ntri, std::size_t& nnodes);
};

}
#endif
