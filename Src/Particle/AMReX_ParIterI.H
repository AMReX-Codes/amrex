
template <int NR, int NI, int NA>
ParIter<NR,NI,NA>::ParIter (ContainerType& pc, int level)
    : 
    MFIter(*pc.m_dummy_mf[level],
           pc.do_tiling ? pc.tile_size : IntVect::TheZeroVector()),
    m_level(level),
    m_pariter_index(0)
{
    auto& aos_map = pc.m_particles[level];
    auto& soa_map = pc.m_arraydata[level];
    
    for (int i = beginIndex; i < endIndex; ++i)
    {
        int grid = (*index_map)[i];
        int tile = local_tile_index_map ? (*local_tile_index_map)[i] : 0;
        auto key = std::make_pair(grid,tile);
        auto f_aos = aos_map.find(key);
        if (f_aos != aos_map.end() && f_aos->second.numParticles() > 0)
        {
            auto f_soa = soa_map.find(key);
            BL_ASSERT(f_soa != soa_map.end());
            m_valid_index.push_back(i);
            m_aos.push_back(&(f_aos->second));
            m_soa.push_back(&(f_soa->second));
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

