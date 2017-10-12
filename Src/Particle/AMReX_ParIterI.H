
template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt>::ParIterBase 
  (ContainerRef pc, int level, MFItInfo& info)
    : 
      MFIter(*pc.m_dummy_mf[level], info.EnableTiling(pc.tile_size)),
      m_level(level),
      m_pariter_index(0)
{
    auto& particles = pc.GetParticles(level);

    int start = dynamic ? 0 : beginIndex;
    for (int i = start; i < endIndex; ++i)
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
        if (dynamic) {
#ifdef _OPENMP            
            int ind = omp_get_thread_num();
            m_pariter_index += ind;
            if (ind < m_valid_index.size()) {
                currentIndex = beginIndex = m_valid_index[ind];
            } else {
                currentIndex = endIndex;
            }
            for (int i = 0; i < omp_get_num_threads(); ++i) {
                m_valid_index.push_back(endIndex);
            }
#endif
        }
        m_valid_index.push_back(endIndex);
    }
}

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt>::ParIterBase 
  (ContainerRef pc, int level)
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

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt>::GetPosition
  (AMREX_D_DECL(Vector<Real>& x, Vector<Real>& y, Vector<Real>& z)) const
{
    const auto& aos = GetArrayOfStructs();
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    AMREX_D_TERM(x.resize(shape.second);,
           y.resize(shape.second);,
           z.resize(shape.second););
    amrex_particle_get_position(p, shape.first, shape.second,
                                AMREX_D_DECL(x.data(), y.data(), z.data()));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetPosition
  (AMREX_D_DECL(const Vector<Real>& x, const Vector<Real>& y, const Vector<Real>& z)) const
{
    auto& aos = this->GetArrayOfStructs();
    BL_ASSERT(AMREX_D_TERM(x.size() == aos.size(), && x.size() == y.size(), && x.size() == z.size()));
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    amrex_particle_set_position(p, shape.first, shape.second,
                                AMREX_D_DECL(x.data(), y.data(), z.data()));
}

