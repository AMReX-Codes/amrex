
template <typename G>
IndexSpaceImp<G>::IndexSpaceImp (const G& gshop, const Geometry& geom,
                                 int required_coarsening_level,
                                 int max_coarsening_level,
                                 int ngrow, bool build_coarse_level_by_coarsening,
                                 bool extend_domain_face)
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
    m_gslevel.emplace_back(this, gshop, geom, EB2::max_grid_size, ngrow_finest, extend_domain_face);

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
                    m_gslevel.emplace_back(this, gshop, cgeom, EB2::max_grid_size, ng, extend_domain_face);
                }
            } else {
                break;
            }
        }
        m_geom.push_back(cgeom);
        m_domain.push_back(cdomain);
        m_ngrow.push_back(ng);
    }

    m_impfunc.reset(new F(gshop.GetImpFunc()));
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
