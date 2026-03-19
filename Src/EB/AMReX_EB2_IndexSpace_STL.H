#ifndef AMREX_EB2_INDEXSPACE_STL_H_
#define AMREX_EB2_INDEXSPACE_STL_H_
#include <AMReX_Config.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level_STL.H>

#include <string>

/**
 * \file AMReX_EB2_IndexSpace_STL.H
 *
 * IndexSpace implementation backed by STL triangulated geometry.
 */

namespace amrex::EB2 {

//! IndexSpace that samples implicit geometry from STL meshes.
class IndexSpaceSTL
    : public IndexSpace
{
public:

   /**
    * \brief Build an EB hierarchy from an STL file.
    *
    * \param stl_file            Path to the STL mesh.
    * \param stl_scale           Uniform scaling applied to STL coordinates.
    * \param stl_center          Translation applied to center the mesh.
    * \param stl_reverse_normal  If non-zero, flip triangle normals.
    * \param geom                Finest-level geometry descriptor.
    * \param required_coarsening_level Number of coarse levels required by the app.
    * \param max_coarsening_level      Maximum coarse level to attempt.
    * \param ngrow               Number of EB ghost cells to build on level 0.
    * \param build_coarse_level_by_coarsening Whether to build coarse levels by coarsening from fine levels.
    * \param extend_domain_face  Extend EB description beyond the domain from the domain faces.
    * \param num_coarsen_opt     Number of coarse probes to try before rebuilding at full resolution.
    * \param bvh_optimization    Enable BVH acceleration for intersection tests.
    * \param support_mvmc        Whether to generate multi-valued/multi-cut data.
    *                            (This path is not supported yet.)
    */
   IndexSpaceSTL (const std::string& stl_file, Real stl_scale,
                  Array<Real,3> const& stl_center, int stl_reverse_normal,
                  const Geometry& geom, int required_coarsening_level,
                  int max_coarsening_level, int ngrow,
                  bool build_coarse_level_by_coarsening,
                  bool extend_domain_face, int num_coarsen_opt,
                  bool bvh_optimization, bool support_mvmc);

    IndexSpaceSTL (IndexSpaceSTL const&) = delete;
    IndexSpaceSTL (IndexSpaceSTL &&) = delete;
    void operator= (IndexSpaceSTL const&) = delete;
    void operator= (IndexSpaceSTL &&) = delete;

    ~IndexSpaceSTL () override = default;

    //! Return the EB level that matches \p geom.Domain().
    [[nodiscard]] const Level& getLevel (const Geometry& geom) const final;
    //! Return the Geometry whose domain equals \p dom.
    [[nodiscard]] const Geometry& getGeometry (const Box& dom) const final;
    //! Domain of the coarsest EB representation.
    [[nodiscard]] const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }
    //! Dynamically add a number of finer levels to the hierarchy.
    void addFineLevels (int num_new_fine_levels) final;
    //! Prepend regular coarse levels below the existing hierarchy.
   void addRegularCoarseLevels (int num_new_coarse_levels) final;

    //! Apply a periodic shift (in cells) to the entire IndexSpace.
    void setShift (int direction, int ncells) override;

private:

    Vector<STLLevel> m_stllevel;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
};

}

#endif
