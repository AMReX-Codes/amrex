
template <typename G>
IndexSpaceImp<G>::IndexSpaceImp (const G& gshop, const Geometry& geom,
                                 int required_coarsening_level,
                                 int max_coarsening_level,
                                 int ngrow, bool build_coarse_level_by_coarsening,
                                 bool extend_domain_face, int num_coarsen_opt)
    : m_gshop(gshop),
      m_build_coarse_level_by_coarsening(build_coarse_level_by_coarsening),
      m_extend_domain_face(extend_domain_face),
      m_num_coarsen_opt(num_coarsen_opt)
{
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
    m_gslevel.reserve(max_coarsening_level+1);
    m_gslevel.emplace_back(this, gshop, geom, EB2::max_grid_size, ngrow_finest, extend_domain_face,
                           num_coarsen_opt);

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
        m_gslevel.emplace_back(this, ilev, EB2::max_grid_size, ng, cgeom, m_gslevel[ilev-1]);
        if (!m_gslevel.back().isOK()) {
            m_gslevel.pop_back();
            if (ilev <= required_coarsening_level) {
                if (build_coarse_level_by_coarsening) {
                    amrex::Abort("Failed to build required coarse EB level "+std::to_string(ilev));
                } else {
                    m_gslevel.emplace_back(this, gshop, cgeom, EB2::max_grid_size, ng, extend_domain_face,
                                           num_coarsen_opt-ilev);
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


template <typename G>
const Level&
IndexSpaceImp<G>::getLevel (const Geometry& geom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), geom.Domain());
    int i = std::distance(m_domain.begin(), it);
    return m_gslevel[i];
}

template <typename G>
const Geometry&
IndexSpaceImp<G>::getGeometry (const Box& dom) const
{
    auto it = std::find(std::begin(m_domain), std::end(m_domain), dom);
    int i = std::distance(m_domain.begin(), it);
    return m_geom[i];
}

template <typename G>
void
IndexSpaceImp<G>::addFineLevels (int num_new_fine_levels)
{
    if (num_new_fine_levels <= 0) { return; }

    if (m_num_coarsen_opt > 0) {
        m_num_coarsen_opt += num_new_fine_levels;
    }

    IndexSpaceImp<G> fine_isp(m_gshop, amrex::refine(m_geom[0], 1<<num_new_fine_levels),
                              num_new_fine_levels-1, num_new_fine_levels-1,
                              m_ngrow[0], m_build_coarse_level_by_coarsening,
                              m_extend_domain_face, m_num_coarsen_opt);

    fine_isp.m_gslevel.reserve(m_domain.size()+num_new_fine_levels);
    for (int i = 0; i < m_domain.size(); ++i) {
        fine_isp.m_gslevel.emplace_back(std::move(m_gslevel[i]));
    }
    std::swap(fine_isp.m_gslevel, m_gslevel);

    m_geom.insert(m_geom.begin(), fine_isp.m_geom.begin(), fine_isp.m_geom.end());
    m_domain.insert(m_domain.begin(), fine_isp.m_domain.begin(), fine_isp.m_domain.end());
    m_ngrow.insert(m_ngrow.begin(), fine_isp.m_ngrow.begin(), fine_isp.m_ngrow.end());
}
