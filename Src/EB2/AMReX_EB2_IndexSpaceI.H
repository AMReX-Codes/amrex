
template <typename G>
IndexSpaceImp<G>::IndexSpaceImp (const G& gshop, const Geometry& geom,
                                 int required_coarsening_level,
                                 int max_coarsening_level)
{
    // build finest level (i.e., level 0) first
    m_geom.push_back(geom);
    m_domain.push_back(geom.Domain());
    m_gslevel.reserve(std::max(max_coarsening_level+1, 1));
    m_gslevel.emplace_back(gshop, geom, EB2::max_grid_size);

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

        Box cdomain = amrex::coarsen(m_geom.back().Domain(),2);
        Geometry cgeom(cdomain);
        m_gslevel.emplace_back(ilev, EB2::max_grid_size, cgeom, m_gslevel[ilev-1]);
        if (!m_gslevel.back().isOK()) {
            m_gslevel.pop_back();
            if (ilev <= required_coarsening_level) {
                amrex::Abort("Failed to build required coarse EB level "+std::to_string(ilev));
            } else {
                break;
            }
        }
        m_geom.push_back(cgeom);
        m_domain.push_back(cdomain);
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

