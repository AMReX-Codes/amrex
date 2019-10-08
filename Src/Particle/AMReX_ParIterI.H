
template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt>::ParIterBase 
  (ContainerRef pc, int level, MFItInfo& info)
    : 
      MFIter(*pc.m_dummy_mf[level], pc.do_tiling ? info.EnableTiling(pc.tile_size) : info),
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
template <typename Container>
void
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt>::GetPosition
(AMREX_D_DECL(Container& x, Container& y, Container& z)) const
{
    const auto& aos = GetArrayOfStructs();
    const auto np = aos.numParticles();

    AMREX_D_TERM(x.resize(np);, y.resize(np);, z.resize(np););
    
    const auto pstruct_ptr = aos().data();

    AMREX_D_TERM(auto x_ptr = x.data();,
                 auto y_ptr = y.data();,
                 auto z_ptr = z.data();)
    
    AMREX_FOR_1D( np, i,
    {
        AMREX_D_TERM(x_ptr[i] = pstruct_ptr[i].pos(0);,
                     y_ptr[i] = pstruct_ptr[i].pos(1);,
                     z_ptr[i] = pstruct_ptr[i].pos(2);)
    });

    Gpu::streamSynchronize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <typename Container>
void
ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetPosition
(AMREX_D_DECL(const Container& x, const Container& y, const Container& z)) const
{
    auto& aos = this->GetArrayOfStructs();
    const auto np = aos.numParticles();

    auto pstruct_ptr = aos().data();

    AMREX_D_TERM(const auto x_ptr = x.data();,
                 const auto y_ptr = y.data();,
                 const auto z_ptr = z.data();)
    
    AMREX_ASSERT(AMREX_D_TERM(x.size() == np, && y.size() == np, && z.size() == np));
    
    AMREX_FOR_1D( np, i,
    {
        AMREX_D_TERM(pstruct_ptr[i].pos(0) = x_ptr[i];,
                     pstruct_ptr[i].pos(1) = y_ptr[i];,
                     pstruct_ptr[i].pos(2) = z_ptr[i];)
    });

    Gpu::streamSynchronize();
}
