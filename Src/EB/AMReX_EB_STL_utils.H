#ifndef AMREX_EB_STL_UTILS_H_
#define AMREX_EB_STL_UTILS_H_

#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Dim3.H>
#include <AMReX_EB2_Graph.H>

namespace amrex
{

class STLtools
{
public:
    struct Triangle {
        XDim3 v1, v2, v3;
    };

    static constexpr int allregular = -1;
    static constexpr int mixedcells = 0;
    static constexpr int allcovered = 1;

private:

    Gpu::PinnedVector<Triangle> m_tri_pts_h;
    Gpu::DeviceVector<Triangle> m_tri_pts_d;
    Gpu::DeviceVector<XDim3> m_tri_normals_d;

    int m_num_tri=0;

    XDim3 m_ptmin;  // All triangles are inside the bounding box defined by
    XDim3 m_ptmax;  //     m_ptmin and m_ptmax.
    XDim3 m_ptref;  // The reference point is slightly outside the bounding box.
    bool m_boundry_is_outside; // Is the bounding box boundary outside or inside the object?

    void read_ascii_stl_file (std::string const& fname, Real scale,
                              Array<Real,3> const& center, int reverse_normal);
    void read_binary_stl_file (std::string const& fname, Real scale,
                               Array<Real,3> const& center, int reverse_normal);

public: // for cuda
    void prepare ();

public:

    void read_stl_file (std::string const& fname, Real scale, Array<Real,3> const& center,
                        int reverse_normal);

    void fill (MultiFab& mf, IntVect const& nghost, Geometry const& geom,
               Real outside_value = -1._rt, Real inside_value = 1._rt) const;

    int getBoxType (Box const& box, Geometry const& geom, RunOn) const;

    bool isGPUable () const noexcept { return true; }

    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn,
                  Box const& bounding_box) const;

    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const& lst, Geometry const& geom, RunOn,
                       Box const& bounding_box) const;

    static void updateIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                          Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                          Array4<Real const> const& lst, Geometry const& geom) ;

};

}
#endif
