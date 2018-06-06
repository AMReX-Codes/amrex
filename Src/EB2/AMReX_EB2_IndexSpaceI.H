
template <typename G>
IndexSpaceImp<G>::IndexSpaceImp (const G& gshop, const Geometry& geom, const Info& info)
{
    // build finest level (i.e., level 0) first
    m_geom.push_back(geom);
    m_domain.push_back(geom.Domain());
    m_gslevel.reserve(std::max(info.max_coarsening_level+1, 1));
    m_gslevel.emplace_back(gshop, geom, info.max_grid_size);

    for (int ilev = 1; ilev <= info.max_coarsening_level; ++ilev)
    {
        Box cdomain = amrex::coarsen(m_geom.back().Domain(),2);
        Geometry cgeom(cdomain);
        m_gslevel.emplace_back(ilev, info.max_grid_size, cgeom, m_gslevel[ilev-1]);
//        if (m_gslevel[ilev].isValid()) {
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

