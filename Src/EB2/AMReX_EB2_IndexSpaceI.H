
template <typename G>
IndexSpaceImp<G>::IndexSpaceImp (G&& gshop, const Geometry& geom, const Info& info)
{
    // build finest level (i.e., level 0) first
    m_gslevel.emplace_back(gshop, geom, info.max_grid_size);

//    for (int 

    m_impfunc.reset(new F(gshop.GetImpFunc()));
}
