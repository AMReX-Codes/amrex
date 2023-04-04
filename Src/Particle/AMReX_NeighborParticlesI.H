
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::use_mask = false;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::enable_inverse = false;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::NeighborParticleContainer (ParGDBBase* gdb, int ncells)
    : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> (gdb),
    m_num_neighbor_cells(ncells)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::NeighborParticleContainer (const Geometry            & geom,
                             const DistributionMapping & dmap,
                             const BoxArray            & ba,
                             int                         nneighbor)
    : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> (geom, dmap, ba),
    m_num_neighbor_cells(nneighbor)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::NeighborParticleContainer (const Vector<Geometry>            & geom,
                             const Vector<DistributionMapping> & dmap,
                             const Vector<BoxArray>            & ba,
                             const Vector<int>                 & rr,
                             int                               nneighbor)
    : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> (geom, dmap, ba, rr),
    m_num_neighbor_cells(nneighbor)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::initializeCommComps () {
    for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal + this->NumRealComps(); ++ii)
        ghost_real_comp.push_back(1);
    for (int ii = 0; ii < 2 + NStructInt + this->NumIntComps(); ++ii)
        ghost_int_comp.push_back(1);
    calcCommSize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::setRealCommComp (int i, bool value) {
    ghost_real_comp[i] = value;
    calcCommSize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::setIntCommComp (int i, bool value) {
    ghost_int_comp[i] = value;
    calcCommSize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::calcCommSize () {
    size_t comm_size = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal + this->NumRealComps(); ++ii) {
        if (ghost_real_comp[ii]) {
            comm_size += sizeof(typename ParticleType::RealType);
        }
    }
    for (int ii = 0; ii < 2 + NStructInt + this->NumIntComps(); ++ii) {
        if (ghost_int_comp[ii]) {
            comm_size += sizeof(int);
        }
    }
    if ( enableInverse() ) comm_size += 4*sizeof(int);
    cdata_size = comm_size;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Regrid (const DistributionMapping &dmap, const BoxArray &ba ) {
    const int lev = 0;
    AMREX_ASSERT(this->finestLevel() == 0);
    this->SetParticleBoxArray(lev, ba);
    this->SetParticleDistributionMap(lev, dmap);
    this->Redistribute();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Regrid (const DistributionMapping &dmap, const BoxArray &ba, int lev) {
    AMREX_ASSERT(lev <= this->finestLevel());
    this->SetParticleBoxArray(lev, ba);
    this->SetParticleDistributionMap(lev, dmap);
    this->Redistribute();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Regrid (const Vector<DistributionMapping>& dmap, const Vector<BoxArray>& ba) {
    AMREX_ASSERT(ba.size() == this->finestLevel()+1);
    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        this->SetParticleBoxArray(lev, ba[lev]);
        this->SetParticleDistributionMap(lev, dmap[lev]);
    }
    this->Redistribute();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::areMasksValid () {

    BL_PROFILE("NeighborParticleContainer::areMasksValid");

    resizeContainers(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        BoxArray ba = this->ParticleBoxArray(lev);
        const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

        if (mask_ptr[lev] == nullptr ||
            ! BoxArray::SameRefs(mask_ptr[lev]->boxArray(), ba) ||
            ! DistributionMapping::SameRefs(mask_ptr[lev]->DistributionMap(), dmap))
        {
            return false;
        }
    }
    return true;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::BuildMasks () {

    BL_PROFILE("NeighborParticleContainer::BuildMasks");

    if (this->numLevels() == 1) use_mask = true;
    else use_mask = false;

    resizeContainers(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        BoxArray ba = this->ParticleBoxArray(lev);
        const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

        const Geometry& geom = this->Geom(lev);

        mask_ptr[lev] = std::make_unique<iMultiFab>(ba, dmap, int(num_mask_comps), m_num_neighbor_cells);
        mask_ptr[lev]->setVal(-1, m_num_neighbor_cells);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*mask_ptr[lev],this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
             mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            (*mask_ptr[lev])[mfi].template setVal<RunOn::Host>(grid_id, box, MaskComps::grid,  1);
            (*mask_ptr[lev])[mfi].template setVal<RunOn::Host>(tile_id, box, MaskComps::tile,  1);
            (*mask_ptr[lev])[mfi].template setVal<RunOn::Host>(lev    , box, MaskComps::level, 1);
        }

        mask_ptr[lev]->FillBoundary(geom.periodicity());
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::GetNeighborCommTags ()
{
    BL_PROFILE("NeighborParticleContainer::GetNeighborCommTags");

    local_neighbors.clear();
    neighbor_procs.clear();

    if (use_mask)
    {
        AMREX_ASSERT(this->finestLevel() == 0);
        const int lev = 0;
        for (MFIter mfi(*mask_ptr[lev],this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
             mfi.isValid(); ++mfi) {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) {
                const int grid = (*mask_ptr[lev])[mfi](iv, MaskComps::grid);
                if (grid >= 0) {
                    const int tile = (*mask_ptr[lev])[mfi](iv, MaskComps::tile);
                    const int level = (*mask_ptr[lev])[mfi](iv, MaskComps::level);
                    const int global_proc = this->ParticleDistributionMap(level)[grid];
                    const int proc = ParallelContext::global_to_local_rank(global_proc);
                    NeighborCommTag comm_tag(proc, level, grid, tile);
                    local_neighbors.push_back(comm_tag);
                    if (proc != ParallelContext::MyProcSub())
                        neighbor_procs.push_back(proc);
                }
            }
        }
    }
    else
    {
        for (int lev = 0; lev < this->numLevels(); ++lev)
        {
            for (MFIter mfi(*mask_ptr[lev],this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
                 mfi.isValid(); ++mfi) {
                const Box& box = mfi.validbox();
                Vector<NeighborCommTag> comm_tags;
                GetCommTagsBox(comm_tags, lev, box);
                for (auto const& tag : comm_tags) {
                    local_neighbors.push_back(tag);
                    if (tag.proc_id != ParallelContext::MyProcSub())
                        neighbor_procs.push_back(tag.proc_id);
                }
            }
        }
    }

    RemoveDuplicates(local_neighbors);
    RemoveDuplicates(neighbor_procs);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
IntVect
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::computeRefFac (int src_lev, int lev)
{
    IntVect ref_fac = IntVect(AMREX_D_DECL(1,1,1));
    if (src_lev < lev) {
        for (int l = src_lev; l < lev; ++l) {
            ref_fac *= this->GetParGDB()->refRatio(l);
        }
    } else if (src_lev > lev) {
        for (int l = src_lev; l > lev; --l) {
            ref_fac *= this->GetParGDB()->refRatio(l-1);
        }
        ref_fac *= -1;
    }
    return ref_fac;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::GetCommTagsBox (Vector<NeighborCommTag>& tags, int src_lev, const Box& in_box)
{
    std::vector< std::pair<int, Box> > isects;
    Box tbx;

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        Box box = in_box;
        const IntVect& ref_fac = computeRefFac(src_lev, lev);
        if (ref_fac < IntVect::TheZeroVector())
        {
            box.coarsen(-1*ref_fac);
        }
        else if (ref_fac > IntVect::TheZeroVector())
        {
            box.refine(ref_fac);
        }
        box.grow(computeRefFac(0, lev)*m_num_neighbor_cells);
        const Periodicity& periodicity = this->Geom(lev).periodicity();
        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
        const BoxArray& ba = this->ParticleBoxArray(lev);

        for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
        {
            const Box& pbox = box + (*pit);
            bool first_only = false;
            ba.intersections(pbox, isects, first_only, 0);
            for (const auto& isec : isects) {
                const int grid = isec.first;
                const int global_proc = this->ParticleDistributionMap(lev)[grid];
                const int proc = ParallelContext::global_to_local_rank(global_proc);
                for (IntVect iv = pbox.smallEnd(); iv <= pbox.bigEnd(); pbox.next(iv))
                {
                    if (ba[grid].contains(iv))
                    {
                        int tile = getTileIndex(iv, ba[grid],
                                                this->do_tiling, this->tile_size, tbx);
                        tags.push_back(NeighborCommTag(proc, lev, grid, tile));
                    }
                }
            }
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::cacheNeighborInfo () {

    BL_PROFILE("NeighborParticleContainer::cacheNeighborInfo");

    AMREX_ASSERT(this->OK());

    resizeContainers(this->numLevels());

    clearNeighbors();

    AMREX_ASSERT(hasNeighbors() == false);

    const int MyProc = ParallelContext::MyProcSub();

    amrex::Vector<std::map<PairIndex,       Vector<NeighborIndexMap> > > local_map;
    std::map<NeighborCommTag, Vector<NeighborIndexMap> > remote_map;

    // tmp data structures used for OMP reduction
    amrex::Vector<std::map<PairIndex,       Vector<Vector<NeighborIndexMap> > > > tmp_local_map;
    std::map<NeighborCommTag, Vector<Vector<NeighborIndexMap> > > tmp_remote_map;

    local_map.resize(this->numLevels());
    tmp_local_map.resize(this->numLevels());

    int num_threads = OpenMP::get_max_threads();

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        // resize our temporaries in serial
        for (int i = 0; i < static_cast<int>(local_neighbors.size()); ++i) {
            const NeighborCommTag& comm_tag = local_neighbors[i];
            tmp_remote_map[comm_tag].resize(num_threads);
            remote_map[comm_tag];
            PairIndex index(comm_tag.grid_id, comm_tag.tile_id);
            tmp_local_map[lev][index].resize(num_threads);
            local_map[lev][index];
            buffer_tag_cache[lev][index].resize(num_threads);
        }
    }

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        // First pass - each thread collects the NeighborIndexMaps it owes to other
        // grids / tiles / procs
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            Vector<NeighborCopyTag> tags;
            tags.reserve(AMREX_D_TERM(3, *3, *3));
            for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
                int thread_num = OpenMP::get_thread_num();
                const int& grid = pti.index();
                const int& tile = pti.LocalTileIndex();
                PairIndex src_index(grid, tile);

                NeighborCopyTag src_tag(lev, grid, tile);

                auto& cache = buffer_tag_cache[lev][src_index][thread_num];

                auto& particles = pti.GetArrayOfStructs();
                for (int i = 0; i < pti.numParticles(); ++i) {
                    const ParticleType& p = particles[i];

                    getNeighborTags(tags, p, m_num_neighbor_cells, src_tag, pti);

                    // Add neighbors to buffers
                    for (int j = 0; j < static_cast<int>(tags.size()); ++j) {
                        NeighborCopyTag& tag = tags[j];
                        PairIndex dst_index(tag.grid, tag.tile);
                        if (tag.grid < 0) continue;

                        tag.src_index = i;
                        const int cache_index = cache.size();
                        cache.push_back(tag);

                        const int global_who = this->ParticleDistributionMap(tag.level)[tag.grid];
                        const int who = ParallelContext::global_to_local_rank(global_who);
                        NeighborIndexMap nim(tag.level, dst_index.first, dst_index.second, -1,
                                             lev, src_index.first, src_index.second,
                                             cache_index, thread_num);
                        if (who == MyProc) {
                            auto& tmp = tmp_local_map[tag.level][dst_index];
                            Vector<NeighborIndexMap>& buffer = tmp[thread_num];
                            buffer.push_back(nim);
                        } else {
                            NeighborCommTag comm_tag(who, tag.level, tag.grid, tag.tile);
                            Vector<NeighborIndexMap>& buffer = tmp_remote_map[comm_tag][thread_num];
                            buffer.push_back(nim);
                        }
                    }
                    tags.clear();
                }
            }
        }
    }

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        // second pass - for each tile, collect the neighbors owed from all threads
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const int grid = mfi.index();
            const int tile = mfi.LocalTileIndex();
            PairIndex index(grid, tile);
            for (int i = 0; i < num_threads; ++i) {
                local_map[lev][index].insert(local_map[lev][index].end(),
                                             tmp_local_map[lev][index][i].begin(),
                                             tmp_local_map[lev][index][i].end());
                tmp_local_map[lev][index][i].erase(tmp_local_map[lev][index][i].begin(),
                                                   tmp_local_map[lev][index][i].end());
            }
        }
    }

    // do the same for the remote neighbors
    typename std::map<NeighborCommTag, Vector<Vector<NeighborIndexMap> > >::iterator it;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#pragma omp single nowait
#endif
    for (it=tmp_remote_map.begin(); it != tmp_remote_map.end(); it++) {
#ifdef AMREX_USE_OMP
#pragma omp task firstprivate(it)
#endif
        {
            const NeighborCommTag& tag = it->first;
            Vector<Vector<NeighborIndexMap> >& tmp = it->second;
            for (int i = 0; i < num_threads; ++i) {
                remote_map[tag].insert(remote_map[tag].end(), tmp[i].begin(), tmp[i].end());
                tmp[i].erase(tmp[i].begin(), tmp[i].end());
            }
        }
    }

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        // now for the local neighbors, allocate buffers and cache
        for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const int grid = mfi.index();
            const int tile = mfi.LocalTileIndex();
            PairIndex dst_index(grid, tile);
            const Vector<NeighborIndexMap>& map = local_map[lev][dst_index];
            const int num_ghosts = map.size();
            neighbors[lev][dst_index].define(this->NumRuntimeRealComps(),
                                             this->NumRuntimeIntComps());
            neighbors[lev][dst_index].resize(num_ghosts);
            local_neighbor_sizes[lev][dst_index] = neighbors[lev][dst_index].size();
        }
    }

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const int grid = mfi.index();
            const int tile = mfi.LocalTileIndex();
            PairIndex dst_index(grid, tile);
            const Vector<NeighborIndexMap>& map = local_map[lev][dst_index];
            const int num_ghosts = map.size();
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < num_ghosts; ++i) {
                const NeighborIndexMap& nim = map[i];
                PairIndex src_index(nim.src_grid, nim.src_tile);
                Vector<NeighborCopyTag>& tags = buffer_tag_cache[nim.src_level][src_index][nim.thread_num];
                AMREX_ASSERT(nim.src_index < tags.size());
                tags[nim.src_index].dst_index = i;
                AMREX_ASSERT(size_t(tags[nim.src_index].dst_index) < neighbors[nim.dst_level][dst_index].size());
            }
        }
    }

    // now we allocate the send buffers and cache the remotes
    std::map<int, int> tile_counts;
    for (const auto& kv: remote_map) {
        tile_counts[kv.first.proc_id] += 1;
    }

    for (const auto& kv: remote_map) {
        if (kv.first.proc_id == MyProc) continue;
        Vector<char>& buffer = send_data[kv.first.proc_id];
        buffer.resize(sizeof(int));
        std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
    }

    for (auto& kv : remote_map) {
        if (kv.first.proc_id == MyProc) continue;
        int np = kv.second.size();
        int data_size = np * cdata_size;
        Vector<char>& buffer = send_data[kv.first.proc_id];
        size_t old_size = buffer.size();
        size_t new_size = buffer.size() + 4*sizeof(int) + data_size;
        buffer.resize(new_size);
        char* dst = &buffer[old_size];
        std::memcpy(dst, &(kv.first.level_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.grid_id ), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.tile_id ), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &data_size,           sizeof(int)); dst += sizeof(int);
        size_t buffer_offset = old_size + 4*sizeof(int);
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < np; ++i) {
            const NeighborIndexMap& nim = kv.second[i];
            PairIndex src_index(nim.src_grid, nim.src_tile);
            Vector<NeighborCopyTag>& tags = buffer_tag_cache[nim.src_level][src_index][nim.thread_num];
            tags[nim.src_index].dst_index = buffer_offset + i*cdata_size;
        }
    }

    if ( enableInverse() )
    {
        for (int lev = 0; lev < this->numLevels(); ++lev)
        {
            for (const auto& kv : neighbors[lev])
            {
                inverse_tags[lev][kv.first].resize(kv.second.size());
            }
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
getNeighborTags (Vector<NeighborCopyTag>& tags, const ParticleType& p,
                 int nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti)
{
    getNeighborTags(tags, p, IntVect(AMREX_D_DECL(nGrow, nGrow, nGrow)), src_tag, pti);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
getNeighborTags (Vector<NeighborCopyTag>& tags, const ParticleType& p,
                 const IntVect& nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti)
{
    Box shrink_box = pti.tilebox();
    shrink_box.grow(-nGrow);

    if (use_mask) {
        const BaseFab<int>& mask = (*mask_ptr[src_tag.level])[src_tag.grid];
        AMREX_ASSERT(this->finestLevel() == 0);
        AMREX_ASSERT(src_tag.level == 0);

        const int lev = 0;
        const IntVect& iv = this->Index(p, lev);
        if (shrink_box.contains(iv)) return;

        const Periodicity& periodicity = this->Geom(lev).periodicity();
        const Box& domain = this->Geom(lev).Domain();
        const IntVect& lo = domain.smallEnd();
        const IntVect& hi = domain.bigEnd();

        // Figure out all our neighbors, removing duplicates
        AMREX_D_TERM(
                 for (int ii = -nGrow[0]; ii < nGrow[0] + 1; ii += nGrow[0]) {,
                     for (int jj = -nGrow[1]; jj < nGrow[1] + 1; jj += nGrow[1]) {,
                         for (int kk = -nGrow[2]; kk < nGrow[2] + 1; kk += nGrow[2]) {)
                             if (AMREX_D_TERM((ii == 0), && (jj == 0), && (kk == 0))) continue;
                             IntVect shift(AMREX_D_DECL(ii, jj, kk));
                             IntVect neighbor_cell = iv + shift;

                             NeighborCopyTag tag;
                             tag.grid  = mask(neighbor_cell, MaskComps::grid);
                             tag.tile  = mask(neighbor_cell, MaskComps::tile);
                             tag.level = mask(neighbor_cell, MaskComps::level);
                             if (periodicity.isAnyPeriodic()) {
                                 for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                     if (! periodicity.isPeriodic(dir)) continue;
                                     if (neighbor_cell[dir] < lo[dir])
                                         tag.periodic_shift[dir] = -1;
                                     else if (neighbor_cell[dir] > hi[dir])
                                         tag.periodic_shift[dir] =  1;
                                 }
                             }

                             if (tag != src_tag) tags.push_back(tag);

                             AMREX_D_TERM(
                                          },
                                 },
                         })

        RemoveDuplicates(tags);
        return;
    }
    else
    {
        std::vector< std::pair<int, Box> > isects;
        Box tbx;
        for (int lev = 0; lev < this->numLevels(); ++lev)
        {
            IntVect ref_fac = computeRefFac(0, lev);
            const Periodicity& periodicity = this->Geom(lev).periodicity();
            const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
            const BoxArray& ba = this->ParticleBoxArray(lev);
            const IntVect& iv = this->Index(p, lev);
            for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
            {
                Box pbox = amrex::grow(Box(iv, iv), ref_fac*nGrow) + (*pit);
                bool first_only = false;
                ba.intersections(pbox, isects, first_only, 0);
                for (const auto& isec : isects)
                {
                    const Box& grid_box = ba[isec.first];
                    for (IntVect cell = pbox.smallEnd(); cell <= pbox.bigEnd(); pbox.next(cell)) {
                        if ( !grid_box.contains(cell) ) continue;
                        int tile = getTileIndex(cell, grid_box,
                                                this->do_tiling, this->tile_size, tbx);
                        auto nbor = NeighborCopyTag(lev, isec.first, tile);
                        nbor.periodic_shift = -1*(*pit);
                        if (src_tag != nbor) tags.push_back(nbor);
                    }
                }
            }
        }

        RemoveDuplicates(tags);
        return;
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::fillNeighbors () {
#ifdef AMREX_USE_GPU
    fillNeighborsGPU();
#else
    fillNeighborsCPU();
#endif
    m_has_neighbors = true;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::sumNeighbors (int real_start_comp, int real_num_comp,
                int int_start_comp,  int int_num_comp) {
#ifdef AMREX_USE_GPU
    amrex::ignore_unused(real_start_comp,real_num_comp,int_start_comp,int_num_comp);
    amrex::Abort("Not implemented.");
#else
    sumNeighborsCPU(real_start_comp, real_num_comp, int_start_comp, int_num_comp);
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::updateNeighbors (bool boundary_neighbors_only)
{

  AMREX_ASSERT(hasNeighbors());

#ifdef AMREX_USE_GPU
    updateNeighborsGPU(boundary_neighbors_only);
#else
    amrex::ignore_unused(boundary_neighbors_only);
    updateNeighborsCPU(true);
#endif
    m_has_neighbors = true;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::clearNeighbors ()
{
#ifdef AMREX_USE_GPU
    clearNeighborsGPU();
#else
    clearNeighborsCPU();
#endif
    m_has_neighbors = false;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class CheckPair>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
buildNeighborList (CheckPair&& check_pair, bool /*sort*/)
{
    AMREX_ASSERT(numParticlesOutOfRange(*this, m_num_neighbor_cells) == 0);

    BL_PROFILE("NeighborParticleContainer::buildNeighborList");

    resizeContainers(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        m_neighbor_list[lev].clear();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            m_neighbor_list[lev][index];
        }

#ifndef AMREX_USE_GPU
        neighbor_list[lev].clear();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            neighbor_list[lev][index];
        }
#endif

              auto& plev = this->GetParticles(lev);
        const auto& geom = this->Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            int gid = pti.index();
            int tid = pti.LocalTileIndex();
            auto index = std::make_pair(gid, tid);

            auto& ptile = plev[index];

            if (ptile.numParticles() == 0) continue;

            Box bx = pti.tilebox();
            int ng = computeRefFac(0, lev).max()*m_num_neighbor_cells;
            bx.grow(ng);

            Gpu::DeviceVector<int> off_bins_v;
            Gpu::DeviceVector<Dim3>      lo_v;
            Gpu::DeviceVector<Dim3>      hi_v;
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> dxi_v;
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> plo_v;

            off_bins_v.push_back(0);
            off_bins_v.push_back(bx.numPts());
            lo_v.push_back(lbound(bx));
            hi_v.push_back(ubound(bx));
            dxi_v.push_back(geom.InvCellSizeArray());
            plo_v.push_back(geom.ProbLoArray());

            m_neighbor_list[lev][index].build(ptile,
                                              std::forward<CheckPair>(check_pair),
                                              off_bins_v, dxi_v, plo_v, lo_v, hi_v, ng);

#ifndef AMREX_USE_GPU
            const auto& counts = m_neighbor_list[lev][index].GetCounts();
            const auto& list   = m_neighbor_list[lev][index].GetList();

            int li = 0;
            for (int i = 0; i < ptile.numParticles(); ++i)
            {
                auto cnt = counts[i];
                neighbor_list[lev][index].push_back(cnt);
                for (size_t j = 0; j < cnt; ++j)
                {
                    neighbor_list[lev][index].push_back(list[li++]+1);
                }
            }
#endif
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class CheckPair, class OtherPCType>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
buildNeighborList (CheckPair&& check_pair, OtherPCType& other,
                   Vector<std::map<std::pair<int, int>, amrex::NeighborList<typename OtherPCType::ParticleType> > >& neighbor_lists,
                   bool /*sort*/)
{
    BL_PROFILE("NeighborParticleContainer::buildNeighborList");

    AMREX_ASSERT(numParticlesOutOfRange(*this, m_num_neighbor_cells) == 0);
    AMREX_ASSERT(numParticlesOutOfRange(other, m_num_neighbor_cells) == 0);
    AMREX_ASSERT(SameIteratorsOK(*this, other));

    EnsureThreadSafeTiles(*this);
    EnsureThreadSafeTiles(other);

    resizeContainers(this->numLevels());
    neighbor_lists.resize(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        neighbor_lists[lev].clear();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            neighbor_lists[lev][index];
        }

              auto& plev = this->GetParticles(lev);
        const auto& geom = this->Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            int gid = pti.index();
            int tid = pti.LocalTileIndex();
            auto index = std::make_pair(gid, tid);

            const auto& ptile = plev[index];
            auto& other_ptile = other.ParticlesAt(lev, pti);
            if (ptile.numParticles() == 0) continue;

            Box bx = pti.tilebox();
            int ng = computeRefFac(0, lev).max()*m_num_neighbor_cells;
            bx.grow(ng);

            Gpu::DeviceVector<int> off_bins_v;
            Gpu::DeviceVector<Dim3>      lo_v;
            Gpu::DeviceVector<Dim3>      hi_v;
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> dxi_v;
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> plo_v;

            off_bins_v.push_back(0);
            off_bins_v.push_back(bx.numPts());
            lo_v.push_back(lbound(bx));
            hi_v.push_back(ubound(bx));
            dxi_v.push_back(geom.InvCellSizeArray());
            plo_v.push_back(geom.ProbLoArray());

            neighbor_lists[lev][index].build(ptile, other_ptile,
                                             std::forward<CheckPair>(check_pair),
                                             off_bins_v, dxi_v, plo_v, lo_v, hi_v, ng);
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class CheckPair>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
buildNeighborList (CheckPair&& check_pair, int type_ind, int* ref_ratio,
                   int num_bin_types, bool /*sort*/)
{
    AMREX_ASSERT(numParticlesOutOfRange(*this, m_num_neighbor_cells) == 0);

    if (num_bin_types == 1) AMREX_ASSERT(ref_ratio[0] == 1);

    BL_PROFILE("NeighborParticleContainer::buildNeighborList");

    resizeContainers(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        m_neighbor_list[lev].clear();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            m_neighbor_list[lev][index];
        }

#ifndef AMREX_USE_GPU
        neighbor_list[lev].clear();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            neighbor_list[lev][index];
        }
#endif

              auto& plev = this->GetParticles(lev);
        const auto& geom = this->Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            int gid = pti.index();
            int tid = pti.LocalTileIndex();
            auto index  = std::make_pair(gid, tid);
            auto& ptile = plev[index];

            if (ptile.numParticles() == 0) continue;

            Box bx = pti.tilebox();
            int ng = 1;

            auto& soa     = pti.GetStructOfArrays();
            auto  TypeVec = soa.GetIntData(type_ind);
            int*  bin_type_array = TypeVec.data();

            Gpu::DeviceVector<int> off_bins_v(num_bin_types+1,0);
            Gpu::DeviceVector<int>    nbins_v(num_bin_types+1,0);
            Gpu::DeviceVector<Dim3>      lo_v(num_bin_types);
            Gpu::DeviceVector<Dim3>      hi_v(num_bin_types);
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> dxi_v(num_bin_types);
            Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> plo_v(num_bin_types);

            for (int type(0); type<num_bin_types; ++type) {
                // Domain, RB, Coord, Per
                Box dom = geom.Domain();
                const Real* plo = geom.ProbLo();
                const Real* phi = geom.ProbHi();
                auto lcoord = geom.Coord();
                Array<int,AMREX_SPACEDIM> lper = geom.isPeriodic();

                // Refined tile box and domain
                Box  lbx(  bx.smallEnd(), bx.bigEnd(), bx.ixType() );
                Box ldom( dom.smallEnd(),dom.bigEnd(),dom.ixType() );
                 lbx.refine( ref_ratio[type] );
                ldom.refine( ref_ratio[type] );

                // Local copy of RB for refined geom
                RealBox lrb(plo,phi);

                // New geometry with refined domain
                Geometry lgeom(ldom,lrb,lcoord,lper);

                // Grow for ghost cells
                int NGhost = ref_ratio[type]*m_num_neighbor_cells;
                lbx.grow(NGhost);

                // Store for memcpy
                int nbins = lbx.numPts();
                Dim3 lo = lbound( lbx );
                Dim3 hi = ubound( lbx );

                auto dxInv = lgeom.InvCellSizeArray();
                auto ploa  = lgeom.ProbLoArray();

#ifdef AMREX_USE_GPU
                Gpu::htod_memcpy_async( dxi_v.data()   + type, dxInv.data(), sizeof(dxInv) );
                Gpu::htod_memcpy_async( plo_v.data()   + type, ploa.data() , sizeof(ploa) );
                Gpu::htod_memcpy_async( lo_v.data()    + type, &lo         , sizeof(lo)    );
                Gpu::htod_memcpy_async( hi_v.data()    + type, &hi         , sizeof(hi)    );
                Gpu::htod_memcpy_async( nbins_v.data() + type, &nbins      , sizeof(nbins) );
#else
                std::memcpy( dxi_v.data()   + type, dxInv.data(), sizeof(dxInv) );
                std::memcpy( plo_v.data()   + type, ploa.data() , sizeof(ploa)  );
                std::memcpy( lo_v.data()    + type, &lo         , sizeof(lo)    );
                std::memcpy( hi_v.data()    + type, &hi         , sizeof(hi)    );
                std::memcpy( nbins_v.data() + type, &nbins      , sizeof(nbins) );
#endif
            }

            Gpu::exclusive_scan(nbins_v.begin(), nbins_v.end(), off_bins_v.begin());

            m_neighbor_list[lev][index].build(ptile,
                                              std::forward<CheckPair>(check_pair),
                                              off_bins_v, dxi_v, plo_v, lo_v, hi_v,
                                              ng, num_bin_types, bin_type_array);

#ifndef AMREX_USE_GPU
              BL_PROFILE_VAR("CPU_CopyNeighborList()",CPUCNL);

              const auto& counts = m_neighbor_list[lev][index].GetCounts();
              const auto& list   = m_neighbor_list[lev][index].GetList();

              int li = 0;
              for (int i = 0; i < ptile.numParticles(); ++i) {
                  auto cnt = counts[i];
                  neighbor_list[lev][index].push_back(cnt);
                  for (size_t j = 0; j < cnt; ++j) {
                      neighbor_list[lev][index].push_back(list[li++]+1);
                  }
              }

              BL_PROFILE_VAR_STOP(CPUCNL);
#endif
        } //ParIter
    } //Lev
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class CheckPair>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
selectActualNeighbors (CheckPair&& check_pair, int num_cells)
{
    BL_PROFILE("NeighborParticleContainer::selectActualNeighbors");
    const auto& geom_fine = this->Geom(0);
    const auto& ba_fine   = this->ParticleBoxArray(0);
    if (ba_fine.size() == 1 && !geom_fine.isAnyPeriodic()) {
        return;
    }

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        // clear previous neighbor particle ids
        if (!m_boundary_particle_ids.empty()) {
          for (auto& keyval: m_boundary_particle_ids[lev]) {
            keyval.second.clear();
          }
        }

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());

            // id of actual particles that need to be sent
            m_boundary_particle_ids[lev][index];
            m_boundary_particle_ids[lev][index].resize(pti.numNeighborParticles());
            auto* p_boundary_particle_ids = m_boundary_particle_ids[lev][index].dataPtr();

            const auto& aos     = pti.GetArrayOfStructs();
            const auto* pstruct = aos().dataPtr();
            const auto ptile_data = this->ParticlesAt(lev, pti).getConstParticleTileData();

            Box box       = pti.validbox();
            Box grownBox  = pti.tilebox();
            grownBox.grow(computeRefFac(0, lev).max()*m_num_neighbor_cells);
            const auto lo = lbound(grownBox);
            const auto hi = ubound(grownBox);

            const auto& geom    = this->Geom(lev);
            const auto  domain  = geom.Domain();
            const auto  dxi     = geom.InvCellSizeArray();
            const auto  plo     = geom.ProbLoArray();

            const size_t  np_real  = pti.numRealParticles();
            const size_t  np_total = aos().size();

            DenseBins<ParticleType> bins;
            bins.build(np_total, pstruct, grownBox,
                      [=] AMREX_GPU_DEVICE (const ParticleType& p) noexcept -> IntVect
                      {
                          AMREX_D_TERM(int i = static_cast<int>(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0]) - lo.x);,
                                       int j = static_cast<int>(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1]) - lo.y);,
                                       int k = static_cast<int>(amrex::Math::floor((p.pos(2)-plo[2])*dxi[2]) - lo.z));
                          AMREX_D_TERM(AMREX_ASSERT(i >= 0);, AMREX_ASSERT(j >= 0);, AMREX_ASSERT(k >= 0));

                          return IntVect(AMREX_D_DECL(i, j, k));
                        });

            auto pperm   = bins.permutationPtr();
            auto poffset = bins.offsetsPtr();

            Gpu::Buffer<unsigned int> np_boundary({0});
            unsigned int* p_np_boundary = np_boundary.data();

            AMREX_FOR_1D ( np_real, i,
            {
                IntVect iv(AMREX_D_DECL(
                    static_cast<int>(amrex::Math::floor((pstruct[i].pos(0)-plo[0])*dxi[0])) - lo.x,
                    static_cast<int>(amrex::Math::floor((pstruct[i].pos(1)-plo[1])*dxi[1])) - lo.y,
                    static_cast<int>(amrex::Math::floor((pstruct[i].pos(2)-plo[2])*dxi[2])) - lo.z));
                auto iv3 = iv.dim3();

                int ix = iv3.x;
                int iy = iv3.y;
                int iz = iv3.z;

                int nx = hi.x-lo.x+1;
                int ny = hi.y-lo.y+1;
                int nz = hi.z-lo.z+1;

                bool isActualNeighbor = false;
                for (int ii = amrex::max(ix-num_cells, 0); ii <= amrex::min(ix+num_cells, nx); ++ii) {
                    for (int jj = amrex::max(iy-num_cells, 0); jj <= amrex::min(iy+num_cells, ny); ++jj) {
                        for (int kk = amrex::max(iz-num_cells, 0); kk <= amrex::min(iz+num_cells, nz); ++kk) {
                            if (isActualNeighbor) break;
                            int nbr_cell_id = (ii * ny + jj) * nz + kk;
                            for (auto p = poffset[nbr_cell_id]; p < poffset[nbr_cell_id+1]; ++p) {
                                if (pperm[p] == i) continue;
                                if (call_check_pair(check_pair, ptile_data, ptile_data, i, pperm[p])) {
                                    IntVect cell_ijk = getParticleCell(pstruct[pperm[p]], plo, dxi, domain);
                                    if (!box.contains(cell_ijk)) {
                                        unsigned int loc = Gpu::Atomic::Add(p_np_boundary, 1u);
                                        p_boundary_particle_ids[loc] = i;
                                        isActualNeighbor = true;
                                        break;
                                    }
                                }// end if check_pair
                            }
                        }
                    }
                }
            });// end amrex_for_1d

            unsigned int* p_np_boundary_h = np_boundary.copyToHost();
            m_boundary_particle_ids[lev][index].resize(*p_np_boundary_h);

        }// end mypariter
    }// end lev
}
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
printNeighborList ()
{
    BL_PROFILE("NeighborParticleContainer::printNeighborList");

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto index = std::make_pair(gid, tid);
            m_neighbor_list[lev][index].print();
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
resizeContainers (int num_levels)
{
    this->reserveData();
    this->resizeData();
    if ( static_cast<int>(neighbors.size()) <= num_levels )
    {
        neighbors.resize(num_levels);
        m_neighbor_list.resize(num_levels);
        neighbor_list.resize(num_levels);
        mask_ptr.resize(num_levels);
        buffer_tag_cache.resize(num_levels);
        local_neighbor_sizes.resize(num_levels);
        if ( enableInverse() ) inverse_tags.resize(num_levels);
    }

    AMREX_ASSERT((neighbors.size() == m_neighbor_list.size()) &&
                 (neighbors.size() == mask_ptr.size()     )    );
}
