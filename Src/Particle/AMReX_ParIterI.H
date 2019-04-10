
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
    AMREX_D_TERM(x.resize(aos.size());,
                 y.resize(aos.size());,
                 z.resize(aos.size()););
    
#ifdef AMREX_USE_CUDA
    thrust::transform(
        thrust::device,
        thrust::device_ptr<const ParticleType>(static_cast<const ParticleType*>(aos().data())),
        thrust::device_ptr<const ParticleType>(static_cast<const ParticleType*>(aos().data()) + aos.size()), 
        thrust::make_zip_iterator(thrust::make_tuple(AMREX_D_DECL(thrust::device_ptr<Real>(x.data()), 
                                                                  thrust::device_ptr<Real>(y.data()), 
                                                                  thrust::device_ptr<Real>(z.data())))),
        functors::getPosition() );
#else
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    amrex_particle_get_position(p, shape.first, shape.second,
                                AMREX_D_DECL(x.data(), y.data(), z.data()));
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <typename Container>
void
ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetPosition
(AMREX_D_DECL(const Container& x, const Container& y, const Container& z)) const
{
    auto& aos = this->GetArrayOfStructs();
    BL_ASSERT(AMREX_D_TERM(x.size() == aos.size(), && x.size() == y.size(), && x.size() == z.size()));
    auto N = aos.size();

#ifdef AMREX_USE_CUDA
    thrust::transform(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(AMREX_D_DECL(thrust::device_ptr<const Real>(x.data()), 
                                                                  thrust::device_ptr<const Real>(y.data()), 
                                                                  thrust::device_ptr<const Real>(z.data())))),
        thrust::make_zip_iterator(thrust::make_tuple(AMREX_D_DECL(thrust::device_ptr<const Real>(x.data() + N), 
                                                                  thrust::device_ptr<const Real>(y.data() + N), 
                                                                  thrust::device_ptr<const Real>(z.data() + N)))),
        thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())),
        functors::setPosition<ParticleType>() );
#else
    const auto  p     = aos.data();
    const auto& shape = aos.dataShape();
    amrex_particle_set_position(p, shape.first, shape.second,
                                AMREX_D_DECL(x.data(), y.data(), z.data()));
#endif
}
