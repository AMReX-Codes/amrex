#include <AMReX_EB2_IndexSpace_STL.H>

namespace amrex::EB2 {

IndexSpaceSTL::IndexSpaceSTL (const std::string& stl_file, Real stl_scale,
                              Array<Real,3> const& stl_center, int stl_reverse_normal,
                              const Geometry& geom, int required_coarsening_level,
                              int max_coarsening_level, int ngrow,
                              bool build_coarse_level_by_coarsening,
                              bool extend_domain_face, int num_coarsen_opt)
{
    Gpu::LaunchSafeGuard lsg(true); // Always use GPU

    STLtools stl_tools;
    stl_tools.read_stl_file(stl_file, stl_scale, stl_center, stl_reverse_normal);

    // build finest level (i.e., level 0) first
    AMREX_ALWAYS_ASSERT(required_coarsening_level >= 0 && required_coarsening_level <= 30);
    max_coarsening_level = std::max(required_coarsening_level,max_coarsening_level);
    max_coarsening_level = std::min(30,max_coarsening_level);

    int ngrow_finest = std::max(ngrow,0);
    for (int i = 1; i <= required_coarsening_level; ++i) {
        ngrow_finest *= 2;
    }

    m_geom.push_back(geom);
    m_domain.push_back(geom.Domain());
    m_ngrow.push_back(ngrow_finest);
    m_stllevel.reserve(max_coarsening_level+1);
    m_stllevel.emplace_back(this, stl_tools, geom, EB2::max_grid_size, ngrow_finest,
                            extend_domain_face, num_coarsen_opt);

    for (int ilev = 1; ilev <= max_coarsening_level; ++ilev)
    {
        bool coarsenable = m_geom.back().Domain().coarsenable(2,2);
        if (!coarsenable) {
            if (ilev <= required_coarsening_level) {
                amrex::Abort("IndexSpaceImp: domain is not coarsenable at level "+std::to_string(ilev));
            } else {
                break;
            }
        }

        int ng = (ilev > required_coarsening_level) ? 0 : m_ngrow.back()/2;

        Box cdomain = amrex::coarsen(m_geom.back().Domain(),2);
        Geometry cgeom = amrex::coarsen(m_geom.back(),2);
        m_stllevel.emplace_back(this, ilev, EB2::max_grid_size, ng, cgeom, m_stllevel[ilev-1]);
        if (!m_stllevel.back().isOK()) {
            m_stllevel.pop_back();
            if (ilev <= required_coarsening_level) {
                if (build_coarse_level_by_coarsening) {
                    amrex::Abort("Failed to build required coarse EB level "+std::to_string(ilev));
                } else {
                    m_stllevel.emplace_back(this, stl_tools, cgeom, EB2::max_grid_size, ng,
                                            extend_domain_face, num_coarsen_opt-ilev);
                }
            } else {
                break;
            }
        }
        m_geom.push_back(cgeom);
        m_domain.push_back(cdomain);
        m_ngrow.push_back(ng);
    }
}

const Level&
IndexSpaceSTL::getLevel (const Geometry& geom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), geom.Domain());
    auto i = std::distance(m_domain.begin(), it);
    return m_stllevel[i];
}

const Geometry&
IndexSpaceSTL::getGeometry (const Box& dom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), dom);
    auto i = std::distance(m_domain.begin(), it);
    return m_geom[i];
}

void
IndexSpaceSTL::addFineLevels (int /*num_new_fine_levels*/)
{
    amrex::Abort("IndexSpaceSTL::addFineLevels: todo");
}

}
