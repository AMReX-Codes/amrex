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

namespace amrex
{

class STLtools
{
public:
    struct Triangle {
        XDim3 v1, v2, v3;

        [[nodiscard]] Real cent (int d) const
        {
            static_assert(sizeof(XDim3) == sizeof(Real)*3);
            return Real(1./3.)*((&v1.x)[d] + (&v2.x)[d] + (&v3.x)[d]);
        }

        [[nodiscard]] std::pair<Real,Real> minmax (int d) const
        {
            static_assert(sizeof(XDim3) == sizeof(Real)*3);
            return std::minmax({(&v1.x)[d], (&v2.x)[d], (&v3.x)[d]});
        }
    };

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

    static constexpr int m_bvh_max_size = 4; // max # of triangles in a leaf node
    static constexpr int m_bvh_max_splits = 4; // max # of children
    static constexpr int m_bvh_max_stack_size = 12; // max depth of the tree

    using Node = BVHNodeT<m_bvh_max_size,m_bvh_max_splits>;

    static constexpr int allregular = -1;
    static constexpr int mixedcells = 0;
    static constexpr int allcovered = 1;

    void setBVHOptimization (bool flag) { m_bvh_optimization = flag; }

    void read_stl_file (std::string const& fname, Real scale, Array<Real,3> const& center,
                        int reverse_normal);

    void fill (MultiFab& mf, IntVect const& nghost, Geometry const& geom,
               Real outside_value = -1._rt, Real inside_value = 1._rt) const;

    [[nodiscard]] int getBoxType (Box const& box, Geometry const& geom, RunOn) const;

    static constexpr bool isGPUable () noexcept { return true; }

    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn,
                  Box const& bounding_box) const;

    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const& lst, Geometry const& geom, RunOn,
                       Box const& bounding_box) const;

    static void updateIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                          Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                          Array4<Real const> const& lst, Geometry const& geom) ;

    void prepare (Gpu::PinnedVector<Triangle> a_tri_pts);  // public for cuda

private:

    bool m_bvh_optimization = true;

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
