
template <int NR, int NI, int NA>
ParIter<NR,NI,NA>::ParIter (ContainerType& pc, int level)
    : 
    MFIter(*pc.m_dummy_mf[level],
           pc.do_tiling ? pc.tile_size : IntVect::TheZeroVector()),
    m_level(level),
    m_pariter_index(0)
{
    auto& particles = pc.GetParticles(level);
    
    for (int i = beginIndex; i < endIndex; ++i)
    {
        int grid = (*index_map)[i];
        int tile = local_tile_index_map ? (*local_tile_index_map)[i] : 0;
        auto key = std::make_pair(grid,tile);
        auto f = particles.find(key);
        if (f != particles.end() && f->second.numParticles() > 0)
        {
            m_valid_index.push_back(i);
            m_particle_tiles.push_back(&(f->second));
        }
    }

    if (m_valid_index.empty())
    {
        endIndex = beginIndex;
    }
    else
    {
        currentIndex = beginIndex = m_valid_index.front();
        m_valid_index.push_back(endIndex);
    }
}

template <int NR, int NI, int NA>
void
ParIter<NR,NI,NA>::SetPosition (D_DECL(const Array<Real>& x,
                                       const Array<Real>& y,
                                       const Array<Real>& z))
{
    auto& aos = GetArrayOfStructs();
    BL_ASSERT(D_TERM(x.size() == aos.size(), && x.size() == y.size(), && x.size() == z.size()));
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    amrex_particle_set_position(p, shape.first, shape.second,
                                D_DECL(x.data(), y.data(), z.data()));
}

template <int NR, int NI, int NA>
void
ParIter<NR,NI,NA>::GetPosition (D_DECL(Array<Real>& x,
                                       Array<Real>& y,
                                       Array<Real>& z)) const
{
    const auto& aos = GetArrayOfStructs();
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    D_TERM(x.resize(shape.second);,
           y.resize(shape.second);,
           z.resize(shape.second););
    amrex_particle_get_position(p, shape.first, shape.second,
                                D_DECL(x.data(), y.data(), z.data()));    
}
