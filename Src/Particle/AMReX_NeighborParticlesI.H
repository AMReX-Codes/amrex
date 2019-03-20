
template <int NStructReal, int NStructInt>
bool NeighborParticleContainer<NStructReal, NStructInt>::use_mask = false;

template <int NStructReal, int NStructInt>
NeighborParticleContainer<NStructReal, NStructInt>
::NeighborParticleContainer(ParGDBBase* gdb, int ncells)
    : ParticleContainer<NStructReal, NStructInt, 0, 0> (gdb),
    m_num_neighbor_cells(ncells)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt>
NeighborParticleContainer<NStructReal, NStructInt>
::NeighborParticleContainer(const Geometry            & geom,
                            const DistributionMapping & dmap,
                            const BoxArray            & ba,
                            int                         ncells)
    : ParticleContainer<NStructReal, NStructInt, 0, 0> (geom, dmap, ba),
    m_num_neighbor_cells(ncells)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt>
NeighborParticleContainer<NStructReal, NStructInt>
::NeighborParticleContainer(const Vector<Geometry>            & geom,
                            const Vector<DistributionMapping> & dmap,
                            const Vector<BoxArray>            & ba,
                            const Vector<int>                 & rr,
                            int                               ncells)
    : ParticleContainer<NStructReal, NStructInt, 0, 0> (geom, dmap, ba, rr),
    m_num_neighbor_cells(ncells)
{
    initializeCommComps();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::initializeCommComps() {
    for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal; ++ii)
        rc[ii] = true;
    for (int ii = 0; ii < 2 + NStructInt; ++ii)
        ic[ii] = true;
    calcCommSize();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::setRealCommComp(int i, bool value) {
    rc[i] = value;
    calcCommSize();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::setIntCommComp(int i, bool value) {
    ic[i] = value;
    calcCommSize();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::calcCommSize() {
    size_t comm_size = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal; ++ii) {
        if (rc[ii]) {
            comm_size += sizeof(typename ParticleType::RealType);
        }
    }
    for (int ii = 0; ii < 2 + NStructInt; ++ii) {
        if (ic[ii]) {
            comm_size += sizeof(int);
        }
    }
    cdata_size = comm_size;
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::Regrid(const DistributionMapping &dmap, const BoxArray &ba ) {
    const int lev = 0;
    AMREX_ASSERT(this->finestLevel() == 0);
    this->SetParticleBoxArray(lev, ba);
    this->SetParticleDistributionMap(lev, dmap);
    this->Redistribute();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::Regrid(const DistributionMapping &dmap, const BoxArray &ba, const int lev) {
    AMREX_ASSERT(lev <= this->finestLevel());
    this->SetParticleBoxArray(lev, ba);
    this->SetParticleDistributionMap(lev, dmap);
    this->Redistribute();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::Regrid(const Vector<DistributionMapping>& dmap, const Vector<BoxArray>& ba) {
    AMREX_ASSERT(ba.size() == this->finestLevel()+1);
    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        this->SetParticleBoxArray(lev, ba[lev]);
        this->SetParticleDistributionMap(lev, dmap[lev]);
    }
    this->Redistribute();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::BuildMasks() {

    BL_PROFILE("NeighborParticleContainer::BuildMasks");

    if (this->numLevels() == 1) use_mask = true;
    else use_mask = false;

    resizeContainers(this->numLevels());

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        BoxArray ba = this->ParticleBoxArray(lev);
        const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

        if (mask_ptr[lev] == nullptr ||
            ! BoxArray::SameRefs(mask_ptr[lev]->boxArray(), ba) ||
            ! DistributionMapping::SameRefs(mask_ptr[lev]->DistributionMap(), dmap))
        {
            const Geometry& geom = this->Geom(lev);

            mask_ptr[lev].reset(new iMultiFab(ba, dmap, num_mask_comps, m_num_neighbor_cells));
            mask_ptr[lev]->setVal(-1, m_num_neighbor_cells);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*mask_ptr[lev],this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
                 mfi.isValid(); ++mfi) {
                const Box& box = mfi.tilebox();
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                mask_ptr[lev]->setVal(grid_id, box, MaskComps::grid,  1);
                mask_ptr[lev]->setVal(tile_id, box, MaskComps::tile,  1);
                mask_ptr[lev]->setVal(lev    , box, MaskComps::level, 1);
            }

            mask_ptr[lev]->FillBoundary(geom.periodicity());
        }
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::GetNeighborCommTags()
{
    BL_PROFILE("NeighborParticleContainer::GetNeighborCommTags");

    local_neighbors.clear();
    neighbor_procs.clear();

    if (use_mask)
    {
        BL_ASSERT(this->finestLevel() == 0);
        const int lev = 0;
        for (MFIter mfi(*mask_ptr[lev],this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
             mfi.isValid(); ++mfi) {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) {
                const int grid = (*mask_ptr[lev])[mfi](iv, MaskComps::grid);
                if (grid >= 0) {
                    const int tile = (*mask_ptr[lev])[mfi](iv, MaskComps::tile);
                    const int level = (*mask_ptr[lev])[mfi](iv, MaskComps::level);
                    const int proc = this->ParticleDistributionMap(level)[grid];
                    NeighborCommTag comm_tag(proc, level, grid, tile);
                    local_neighbors.push_back(comm_tag);
                    if (proc != ParallelDescriptor::MyProc())
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
                    if (tag.proc_id != ParallelDescriptor::MyProc())
                        neighbor_procs.push_back(tag.proc_id);
                }
            }
        }
    }

    RemoveDuplicates(local_neighbors);
    RemoveDuplicates(neighbor_procs);
}

template <int NStructReal, int NStructInt>
IntVect
NeighborParticleContainer<NStructReal, NStructInt>
::computeRefFac(const int src_lev, const int lev)
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

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::GetCommTagsBox(Vector<NeighborCommTag>& tags, const int src_lev, const Box& in_box)
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
        box.grow(computeRefFac(0, src_lev)*m_num_neighbor_cells);
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
                const int proc = this->ParticleDistributionMap(lev)[grid];
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

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::cacheNeighborInfo() {

    BL_PROFILE("NeighborParticleContainer::cacheNeighborInfo");

    AMREX_ASSERT(this->OK());

    resizeContainers(this->numLevels());

    clearNeighbors();

    const int MyProc = ParallelDescriptor::MyProc();

    amrex::Vector<std::map<PairIndex,       Vector<NeighborIndexMap> > > local_map;
    std::map<NeighborCommTag, Vector<NeighborIndexMap> > remote_map;

    // tmp data structures used for OMP reduction
    amrex::Vector<std::map<PairIndex,       Vector<Vector<NeighborIndexMap> > > > tmp_local_map;
    std::map<NeighborCommTag, Vector<Vector<NeighborIndexMap> > > tmp_remote_map;

    local_map.resize(this->numLevels());
    tmp_local_map.resize(this->numLevels());

    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    num_threads = omp_get_num_threads();
#endif

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
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Vector<NeighborCopyTag> tags;
            tags.reserve(AMREX_D_TERM(3, *3, *3));
            for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
#ifdef _OPENMP
                int thread_num = omp_get_thread_num();
#else
                int thread_num = 0;
#endif
                const int& grid = pti.index();
                const int& tile = pti.LocalTileIndex();
                PairIndex src_index(grid, tile);

                NeighborCopyTag src_tag(lev, grid, tile);

                auto& cache = buffer_tag_cache[lev][src_index][thread_num];

                auto& particles = pti.GetArrayOfStructs();
                for (unsigned i = 0; i < pti.numParticles(); ++i) {
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

                        const int who = this->ParticleDistributionMap(tag.level)[tag.grid];
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
#ifdef _OPENMP
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
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
    for (it=tmp_remote_map.begin(); it != tmp_remote_map.end(); it++) {
#ifdef _OPENMP
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < num_ghosts; ++i) {
                const NeighborIndexMap& nim = map[i];
                PairIndex src_index(nim.src_grid, nim.src_tile);
                Vector<NeighborCopyTag>& tags = buffer_tag_cache[nim.src_level][src_index][nim.thread_num];
                AMREX_ASSERT(nim.src_index < tags.size());
                tags[nim.src_index].dst_index = i;
                AMREX_ASSERT(tags[nim.src_index].dst_index < neighbors[nim.dst_level][dst_index].size());
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < np; ++i) {
            const NeighborIndexMap& nim = kv.second[i];
            PairIndex src_index(nim.src_grid, nim.src_tile);
            Vector<NeighborCopyTag>& tags = buffer_tag_cache[nim.src_level][src_index][nim.thread_num];
            tags[nim.src_index].dst_index = buffer_offset + i*cdata_size;
        }
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
getNeighborTags(Vector<NeighborCopyTag>& tags, const ParticleType& p,
                const int nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti)
{
    getNeighborTags(tags, p, IntVect(AMREX_D_DECL(nGrow, nGrow, nGrow)), src_tag, pti);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
getNeighborTags(Vector<NeighborCopyTag>& tags, const ParticleType& p,
                const IntVect& nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti)
{
    Box shrink_box = pti.tilebox();
    shrink_box.grow(-nGrow);

    if (use_mask) {
        const BaseFab<int>& mask = (*mask_ptr[src_tag.level])[src_tag.grid];
        BL_ASSERT(this->finestLevel() == 0);
        BL_ASSERT(src_tag.level == 0);

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
                             if (AMREX_D_TERM((ii == 0), and (jj == 0), and (kk == 0))) continue;
                             IntVect shift(AMREX_D_DECL(ii, jj, kk));
                             IntVect neighbor_cell = iv + shift;

                             NeighborCopyTag tag;
                             tag.grid  = mask(neighbor_cell, MaskComps::grid);
                             tag.tile  = mask(neighbor_cell, MaskComps::tile);
                             tag.level = mask(neighbor_cell, MaskComps::level);
                             if (periodicity.isAnyPeriodic()) {
                                 for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                     if (not periodicity.isPeriodic(dir)) continue;
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

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::fillNeighbors() {
#ifdef AMREX_USE_CUDA
    fillNeighborsGPU();
#else
    fillNeighborsCPU();
#endif
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::updateNeighbors()
{
#ifdef AMREX_USE_CUDA
    updateNeighborsGPU();
#else
    updateNeighborsCPU(true);
#endif
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::clearNeighbors()
{
#ifdef AMREX_USE_CUDA
    clearNeighborsGPU();
#else
    clearNeighborsCPU();
#endif
}

template <int NStructReal, int NStructInt>
template <class CheckPair>
void
NeighborParticleContainer<NStructReal, NStructInt>::
buildNeighborList(CheckPair check_pair, bool sort) 
{
#ifdef AMREX_USE_CUDA
    buildNeighborListGPU(check_pair);
#else
    buildNeighborListCPU(check_pair, sort);
#endif
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
printNeighborList()
{
#ifdef AMREX_USE_CUDA
    printNeighborListGPU();
#else
    printNeighborListCPU("nl");
#endif
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
resizeContainers(const int num_levels)
{
    if ( static_cast<int>(neighbors.size()) <= num_levels )
    {
        neighbors.resize(num_levels);
        neighbor_list.resize(num_levels);
        mask_ptr.resize(num_levels);
        buffer_tag_cache.resize(num_levels);
        local_neighbor_sizes.resize(num_levels);
    }

    AMREX_ASSERT((neighbors.size() == neighbor_list.size()) and
                 (neighbors.size() == mask_ptr.size()     )    );
}
