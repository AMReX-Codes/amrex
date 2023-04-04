#include <AMReX_EB2_IndexSpace_chkpt_file.H>

namespace amrex::EB2 {

IndexSpaceChkptFile::IndexSpaceChkptFile (const ChkptFile& chkpt_file,
                                          const Geometry& geom, int required_coarsening_level,
                                          int max_coarsening_level, int ngrow,
                                          bool build_coarse_level_by_coarsening,
                                          bool extend_domain_face)
{
    Gpu::LaunchSafeGuard lsg(true); // Always use GPU

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
    m_chkpt_file_level.reserve(max_coarsening_level+1);
    m_chkpt_file_level.emplace_back(this, chkpt_file, geom, EB2::max_grid_size, ngrow_finest,
            extend_domain_face);

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
        m_chkpt_file_level.emplace_back(this, ilev, EB2::max_grid_size, ng, cgeom, m_chkpt_file_level[ilev-1]);
        if (!m_chkpt_file_level.back().isOK()) {
            m_chkpt_file_level.pop_back();
            if (ilev <= required_coarsening_level) {
                if (build_coarse_level_by_coarsening) {
                    amrex::Abort("Failed to build required coarse EB level "+std::to_string(ilev));
                } else {
                    amrex::Abort("Chkptfile only stored for finest level. Failed to build "+std::to_string(ilev));
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
IndexSpaceChkptFile::getLevel (const Geometry& geom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), geom.Domain());
    auto i = std::distance(m_domain.begin(), it);
    return m_chkpt_file_level[i];
}

const Geometry&
IndexSpaceChkptFile::getGeometry (const Box& dom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), dom);
    auto i = std::distance(m_domain.begin(), it);
    return m_geom[i];
}

void
IndexSpaceChkptFile::addFineLevels (int /*num_new_fine_levels*/)
{
    amrex::Abort("IndexSpaceChkptFile::addFineLevels: not supported");
}

}
