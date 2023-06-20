#ifndef AMREX_EB2_INDEXSPACE_STL_H_
#define AMREX_EB2_INDEXSPACE_STL_H_
#include <AMReX_Config.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level_STL.H>

#include <string>

namespace amrex::EB2 {

class IndexSpaceSTL
    : public IndexSpace
{
public:

   IndexSpaceSTL (const std::string& stl_file, Real stl_scale,
                  Array<Real,3> const& stl_center, int stl_reverse_normal,
                  const Geometry& geom, int required_coarsening_level,
                  int max_coarsening_level, int ngrow,
                  bool build_coarse_level_by_coarsening,
                  bool extend_domain_face, int num_coarsen_opt);

    IndexSpaceSTL (IndexSpaceSTL const&) = delete;
    IndexSpaceSTL (IndexSpaceSTL &&) = delete;
    void operator= (IndexSpaceSTL const&) = delete;
    void operator= (IndexSpaceSTL &&) = delete;

    ~IndexSpaceSTL () override = default;

    [[nodiscard]] const Level& getLevel (const Geometry& geom) const final;
    [[nodiscard]] const Geometry& getGeometry (const Box& dom) const final;
    [[nodiscard]] const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }
    void addFineLevels (int num_new_fine_levels) final;

private:

    Vector<STLLevel> m_stllevel;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
};

}

#endif
