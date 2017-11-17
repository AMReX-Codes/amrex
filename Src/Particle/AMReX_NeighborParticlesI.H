#include "AMReX_NeighborParticles.H"

#include "AMReX_Particles_F.H"

using namespace amrex;

template <int NStructReal, int NStructInt, int NNeighborReal>
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::NeighborParticleContainer(ParGDBBase* gdb, int ncells)
    : ParticleContainer<NStructReal, NStructInt, 0, 0> (gdb),
      num_neighbor_cells(ncells)
{
}


template <int NStructReal, int NStructInt, int NNeighborReal>
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::NeighborParticleContainer(const Geometry            & geom,
                            const DistributionMapping & dmap,
                            const BoxArray            & ba,
                            int                         ncells)
    : ParticleContainer<NStructReal, NStructInt, 0, 0> (geom, dmap, ba),
      num_neighbor_cells(ncells)
{
    BuildLevelMask(0, geom, dmap, ba);
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::Regrid(const DistributionMapping &dmap, const BoxArray &ba ) {
    const int lev = 0;
    this->SetParticleBoxArray(lev, ba);
    this->SetParticleDistributionMap(lev, dmap);
    BuildLevelMask(lev, this->Geom(lev), dmap, ba);
    this->Redistribute();
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::BuildLevelMask(int lev,
                 const Geometry &geom, 
                 const DistributionMapping &dmap, 
                 const BoxArray &ba ) {
    
    BL_PROFILE("NeighborParticleContainer::BuildLevelMask");
    BL_ASSERT(lev == 0);
    BL_ASSERT(dmap == this->ParticleDistributionMap(lev));
    BL_ASSERT(ba == this->ParticleBoxArray(lev));

    if (mask_ptr == nullptr ||
        ! BoxArray::SameRefs(mask_ptr->boxArray(), this->ParticleBoxArray(lev)) ||
        ! DistributionMapping::SameRefs(mask_ptr->DistributionMap(), this->ParticleDistributionMap(lev))) {        
        mask_ptr.reset(new FabArray<BaseFab<int> > (ba, dmap, 2, num_neighbor_cells));
        mask_ptr->setVal(-1, num_neighbor_cells);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*mask_ptr, this->tile_size); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            mask_ptr->setVal(grid_id, box, 0, 1);
            mask_ptr->setVal(tile_id, box, 1, 1);
        }
        
        mask_ptr->FillBoundary(geom.periodicity());

        local_neighbors.clear();
        neighbor_procs.clear();
        for (MFIter mfi(*mask_ptr, this->tile_size); mfi.isValid(); ++mfi) {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) {
                const int grid = (*mask_ptr)[mfi](iv, 0);
                const int tile = (*mask_ptr)[mfi](iv, 1);
                const int proc = this->ParticleDistributionMap(lev)[grid];
                NeighborCommTag comm_tag(proc, grid, tile);
                local_neighbors.push_back(comm_tag);
                if (proc != ParallelDescriptor::MyProc() and grid >= 0) {
                    neighbor_procs.push_back(proc);
                }
            }
        }

        RemoveDuplicates(local_neighbors);
        RemoveDuplicates(neighbor_procs);
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::cacheNeighborInfo(int lev) {
    
    BL_PROFILE("NeighborParticleContainer::cacheNeighborInfo");
    BL_ASSERT(lev == 0);

    clearNeighbors(lev);
    
    const int MyProc = ParallelDescriptor::MyProc();
    const DistributionMapping& dmap = this->ParticleDistributionMap(lev);
    const Periodicity& periodicity = this->Geom(lev).periodicity();
    const RealBox& prob_domain = this->Geom(lev).ProbDomain();
    const Box& domain = this->Geom(lev).Domain();
    const IntVect& lo = domain.smallEnd();
    const IntVect& hi = domain.bigEnd();
    const int nc = num_neighbor_cells;
    
    std::map<PairIndex,       Vector<NeighborIndexMap> > local_map;
    std::map<NeighborCommTag, Vector<NeighborIndexMap> > remote_map;
    
    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    num_threads = omp_get_num_threads();
#endif

    // tmp data structures used for OMP reduction
    std::map<PairIndex,       Vector<Vector<NeighborIndexMap> > > tmp_local_map;
    std::map<NeighborCommTag, Vector<Vector<NeighborIndexMap> > > tmp_remote_map;

    // resize our temporaries in serial
    for (int i = 0; i < static_cast<int>(local_neighbors.size()); ++i) {
        const NeighborCommTag& comm_tag = local_neighbors[i];
        tmp_remote_map[comm_tag].resize(num_threads);
        remote_map[comm_tag];
        PairIndex index(comm_tag.grid_id, comm_tag.tile_id);
        tmp_local_map[index].resize(num_threads);
        local_map[index];
        buffer_tag_cache[index].resize(num_threads);
    }

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
            const BaseFab<int>& mask = (*mask_ptr)[grid];

            auto& cache = buffer_tag_cache[src_index][thread_num];

            Box shrink_box = pti.tilebox();
            shrink_box.grow(-num_neighbor_cells);
            
            auto& particles = pti.GetArrayOfStructs();
            for (unsigned i = 0; i < pti.numParticles(); ++i) {
                const ParticleType& p = particles[i];
                const IntVect& iv = this->Index(p, lev);
                
                // if the particle is more than one cell away from 
                // the tile boundary, its not anybody's neighbor
                if (shrink_box.contains(iv)) continue;
                
                // Figure out all our neighbors, removing duplicates
                AMREX_D_TERM(
                for (int ii = -nc; ii < nc + 1; ii += nc) {,
                    for (int jj = -nc; jj < nc + 1; jj += nc) {,
                        for (int kk = -nc; kk < nc + 1; kk += nc) {)
                            if (AMREX_D_TERM((ii == 0), and (jj == 0), and (kk == 0))) continue;
                            IntVect shift(AMREX_D_DECL(ii, jj, kk));
                            IntVect neighbor_cell = iv + shift;

                            NeighborCopyTag tag;
                            tag.grid = mask(neighbor_cell, 0);
                            tag.tile = mask(neighbor_cell, 1);
                            for (int dim = 0; dim < 3; ++dim)
                                tag.periodic_shift[dim] = 0;
                            if (periodicity.isAnyPeriodic()) {
                                for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
                                    if (not periodicity.isPeriodic(dir)) continue;
                                    if (neighbor_cell[dir] < lo[dir]) 
                                        tag.periodic_shift[dir] = -1;
                                    else if (neighbor_cell[dir] > hi[dir]) 
                                        tag.periodic_shift[dir] =  1;
                                }
                            }

                            if ((grid == tag.grid)           and
                                (tile == tag.tile)           and
                                (tag.periodic_shift[0] == 0) and
                                (tag.periodic_shift[1] == 0) and
                                (tag.periodic_shift[2] == 0)) continue;

                            tags.push_back(tag);

                    AMREX_D_TERM(
                        },
                    },
                })
                    
                RemoveDuplicates(tags);
                
                // Add neighbors to buffers
                for (int j = 0; j < static_cast<int>(tags.size()); ++j) {
                    NeighborCopyTag& tag = tags[j];
                    PairIndex dst_index(tag.grid, tag.tile);
                    if (tag.grid < 0) continue;

                    tag.src_index = i;
                    const int cache_index = cache.size();
                    cache.push_back(tag);
                    
                    const int who = dmap[tag.grid];
                    NeighborIndexMap nim(dst_index.first, dst_index.second, -1,
                                         src_index.first, src_index.second, 
                                         cache_index, thread_num);
                    if (who == MyProc) {
                        auto& tmp = tmp_local_map[dst_index];
                        Vector<NeighborIndexMap>& buffer = tmp[thread_num];
                        buffer.push_back(nim);
                    } else {
                        NeighborCommTag comm_tag(who, tag.grid, tag.tile);
                        Vector<NeighborIndexMap>& buffer = tmp_remote_map[comm_tag][thread_num];
                        buffer.push_back(nim);
                    }
                }
                tags.clear();
            }
        }
    }

    // second pass - for each tile, collect the neighbors owed from all threads
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        const int tile = mfi.LocalTileIndex();
        PairIndex index(grid, tile);
        for (int i = 0; i < num_threads; ++i) {
            local_map[index].insert(local_map[index].end(),
                                    tmp_local_map[index][i].begin(),
                                    tmp_local_map[index][i].end());
            tmp_local_map[index][i].erase(tmp_local_map[index][i].begin(),
                                          tmp_local_map[index][i].end());
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

    // remove any empty map entries from the remote_map
    for (auto pmap_it = remote_map.begin(); pmap_it != remote_map.end(); /* no ++ */) {
        if (pmap_it->second.empty()) {
            remote_map.erase(pmap_it++);
        }
        else {
            ++pmap_it;
        }
    }

    // now for the local neighbors, allocate buffers and cache
    for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        const int tile = mfi.LocalTileIndex();
        PairIndex dst_index(grid, tile);
        const Vector<NeighborIndexMap>& map = local_map[dst_index];
        const int num_ghosts = map.size();
        neighbors[dst_index].resize(num_ghosts * pdata_size);
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int i = 0; i < num_ghosts; ++i) {
            const NeighborIndexMap& nim = map[i];
            PairIndex src_index(nim.src_grid, nim.src_tile);
            Vector<NeighborCopyTag>& tags = buffer_tag_cache[src_index][nim.thread_num];
            BL_ASSERT(nim.src_index < tags.size());
            tags[nim.src_index].dst_offset = i * pdata_size;
            BL_ASSERT(tags[nim.src_index].dst_offset < neighbors[dst_index].size());
        }
    }
    
    for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        const int tile = mfi.LocalTileIndex();
        PairIndex dst_index(grid, tile);
        local_neighbor_sizes[dst_index] = neighbors[dst_index].size();
    }

    // now we allocate the send buffers and cache the remotes
    std::map<int, int> tile_counts;
    for (const auto& kv: remote_map) {
        tile_counts[kv.first.proc_id] += 1;
    }
    
    for (const auto& kv: remote_map) {
        Vector<char>& buffer = send_data[kv.first.proc_id];
        buffer.resize(sizeof(int));
        std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
    }
    
    for (auto& kv : remote_map) {
        int np = kv.second.size();
        int data_size = np * pdata_size;
        Vector<char>& buffer = send_data[kv.first.proc_id];
        size_t old_size = buffer.size();
        size_t new_size = buffer.size() + 2*sizeof(int) + sizeof(int) + data_size;
        buffer.resize(new_size);
        char* dst = &buffer[old_size];
        std::memcpy(dst, &(kv.first.grid_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.tile_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &data_size,          sizeof(int)); dst += sizeof(int);        
        size_t buffer_offset = old_size + 2*sizeof(int) + sizeof(int);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < np; ++i) {
            const NeighborIndexMap& nim = kv.second[i];
            PairIndex src_index(nim.src_grid, nim.src_tile);
            Vector<NeighborCopyTag>& tags = buffer_tag_cache[src_index][nim.thread_num];
            tags[nim.src_index].dst_offset = buffer_offset + i*pdata_size;
        }
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::fillNeighbors(int lev) {
    BL_PROFILE("NeighborParticleContainer::fillNeighbors");
    cacheNeighborInfo(lev);    
    updateNeighbors(lev, false);
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::updateNeighbors(int lev, bool reuse_rcv_counts) {
    
    BL_PROFILE("NeighborParticleContainer::updateNeighbors");
    BL_ASSERT(lev == 0);

    const int MyProc = ParallelDescriptor::MyProc();
    const Periodicity& periodicity = this->Geom(lev).periodicity();
    const RealBox& prob_domain = this->Geom(lev).ProbDomain();
    const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    num_threads = omp_get_num_threads();
#endif
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        PairIndex src_index(pti.index(), pti.LocalTileIndex());
        auto& particles = pti.GetArrayOfStructs();
        for (int j = 0; j < num_threads; ++j) {
            auto& tags = buffer_tag_cache[src_index][j];
            int num_tags = tags.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned i = 0; i < num_tags; ++i) {
                const NeighborCopyTag& tag = tags[i];
                const int who = dmap[tag.grid];
                ParticleType p = particles[tag.src_index];
                if (periodicity.isAnyPeriodic()) {
                    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
                        if (not periodicity.isPeriodic(dir)) continue;
                        if (tag.periodic_shift[dir] == -1) 
                            p.pos(dir) += prob_domain.length(dir);
                        else if (tag.periodic_shift[dir] ==  1) 
                            p.pos(dir) -= prob_domain.length(dir);
                    }
                }
                if (who == MyProc) {
                    PairIndex dst_index(tag.grid, tag.tile);
                    Vector<char>& buffer = neighbors[dst_index];
                    BL_ASSERT(tag.dst_offset < buffer.size());
                    std::memcpy(&buffer[tag.dst_offset], &p, pdata_size);
                } else {
                    std::memcpy(&send_data[who][tag.dst_offset], &p, pdata_size);
                }
            }
        }
    }

    for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        const int tile = mfi.LocalTileIndex();
        PairIndex dst_index(grid, tile);
        neighbors[dst_index].resize(local_neighbor_sizes[dst_index]);
    }
    
    fillNeighborsMPI(reuse_rcv_counts);
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::clearNeighbors(int lev) 
{

    BL_PROFILE("NeighborParticleContainer::clearNeighbors");
    BL_ASSERT(lev == 0);

    neighbors.clear();
    buffer_tag_cache.clear();
    send_data.clear();
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
getRcvCountsMPI() {

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();

    BL_ASSERT(send_data.size() == neighbor_procs.size());

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Vector<long> snds(NProcs, 0);
    rcvs.resize(NProcs, 0);

    num_snds = 0;
    for (const auto& kv : send_data) {
        num_snds      += kv.second.size();
        snds[kv.first] = kv.second.size();
    }
    ParallelDescriptor::ReduceLongMax(num_snds);
    if (num_snds == 0) return;

    const int num_rcvs = neighbor_procs.size();
    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Post receives
    for (int i = 0; i < num_rcvs; ++i) {
        const int Who = neighbor_procs[i];
        const long Cnt = 1;

        BL_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&rcvs[Who], Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (int i = 0; i < num_rcvs; ++i) {
        const int Who = neighbor_procs[i];
        const long Cnt = 1;

        BL_ASSERT(Who >= 0 && Who < NProcs);

        ParallelDescriptor::Send(&snds[Who], Cnt, Who, SeqNum);        
    }

    if (num_rcvs > 0) {
        BL_MPI_REQUIRE( MPI_Waitall(num_rcvs, rreqs.data(), stats.data()) );
    }

#endif // BL_USE_MPI
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
fillNeighborsMPI(bool reuse_rcv_counts) {
    
    BL_PROFILE("NeighborParticleContainer::fillNeighborsMPI");
    
#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    
    // each proc figures out how many bytes it will send, and how
    // many it will receive
    if (!reuse_rcv_counts) getRcvCountsMPI();
    if (num_snds == 0) return;
    
    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);
    
    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.
    Vector<char> recvdata(TotRcvBytes);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }
    
    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0) {
        BL_MPI_REQUIRE( MPI_Waitall(nrcvs, rreqs.data(), stats.data()) );
        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, gid, tid, size;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);
                
                if (size == 0) continue;
                
                PairIndex dst_index(gid, tid);
                size_t old_size = neighbors[dst_index].size();
                size_t new_size = neighbors[dst_index].size() + size;
                neighbors[dst_index].resize(new_size);
                std::memcpy(&neighbors[dst_index][old_size], buffer, size); buffer += size;
            }
        }
    }
#endif
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
buildNeighborList(int lev, bool sort) {
    
    BL_PROFILE("NeighborParticleContainer::buildNeighborList");
    BL_ASSERT(lev == 0);

    neighbor_list.clear();
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        neighbor_list[index] = Vector<int>();
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

    Vector<IntVect> cells;
    Vector<ParticleType> tmp_particles;
    BaseFab<int> head;
    Vector<int>  list;

    for (MyParIter pti(*this, lev, MFItInfo().SetDynamic(true)); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        Vector<int>& nl = neighbor_list[index];
        AoS& particles = pti.GetArrayOfStructs();

        int Np = particles.size();
        int Nn = neighbors[index].size() / pdata_size;
        int N = Np + Nn;

        cells.resize(N);
        tmp_particles.resize(N);
        std::memcpy(&tmp_particles[0], particles.data(), Np*sizeof(ParticleType));
        if (Nn > 0)
           std::memcpy(&tmp_particles[Np], neighbors[index].dataPtr(), Nn*pdata_size); 

        // For each cell on this tile, we build linked lists storing the
        // indices of the particles belonging to it.
        Box box = pti.tilebox();
        box.grow(num_neighbor_cells);
        head.resize(box);
        head.setVal(-1);
        list.resize(N, -1);

        for (int i = 0; i < N; ++i) {
            const ParticleType& p = tmp_particles[i];
            const IntVect& cell = this->Index(p, lev);
            cells[i] = cell;
            list[i] = head(cell);
            head(cell) = i;
        }
        
        // using these lists, we build a neighbor list containing both
        // kinds of particles.
        int p_start_index = 0;
        for (unsigned i = 0; i < N; ++i) {
            const ParticleType& p = tmp_particles[i];
            
            int num_neighbors = 0;
            nl.push_back(0);
            
            const IntVect& cell = cells[i];
            Box bx(cell, cell);
            bx.grow(num_neighbor_cells);
            bx &= head.box();
            
            for (IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
                int j = head(iv);
                while (j >= 0) {
                    if (i == j) {
                        j = list[j];
                        continue;
                    }
                    if ( check_pair(p, tmp_particles[j]) ) {
                        nl.push_back(j+1);
                        num_neighbors += 1;
                    }
                    j = list[j];
                }
            }

            nl[p_start_index] = num_neighbors;
            p_start_index += num_neighbors + 1;
        }

        if (sort) {
            for (unsigned i = 0; i < nl.size(); i += nl[i] +1) {
                std::sort(nl.begin() + i + 1,
                          nl.begin() + nl[i] + i + 1);
            }
        }
    }
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
buildNeighborListFort(int lev, bool sort) {
    
    BL_PROFILE("NeighborParticleContainer::buildNeighborListFort");
    BL_ASSERT(lev == 0);

    neighbor_list.clear();

    const Geometry& gm  = this->Geom(lev);
    const Real*     plo = gm.ProbLo();
    const Real*     dx  = gm.CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        
        PairIndex index(pti.index(), pti.LocalTileIndex());
        Vector<int>& nl = neighbor_list[index];
        AoS& particles = pti.GetArrayOfStructs();

        int Np = particles.size();
        int Nn = neighbors[index].size() / pdata_size;
        int Ns = particles.dataShape().first;
        int N = Np + Nn;

        Vector<ParticleType> tmp_particles(N);
        std::memcpy(&tmp_particles[0], particles.data(), Np*sizeof(ParticleType));

        for (int i = 0; i < Nn; ++i) {
            std::memcpy(&tmp_particles[i + Np],
                        neighbors[index].dataPtr() + i*pdata_size,
                        pdata_size);
        }

        int size_nl = 0;
        nl.resize(100000);

        Box box = pti.tilebox();
        box.grow(num_neighbor_cells);        
        const int* lo = box.loVect();
        const int* hi = box.hiVect();

        amrex_build_nl(tmp_particles.data(), Ns, N, nl.data(), &size_nl, 
                       lo, hi, plo, dx);

        nl.resize(size_nl);

        if (sort) {
            for (unsigned i = 0; i < nl.size(); i += nl[i] +1) {
                std::sort(nl.begin() + i + 1,
                          nl.begin() + nl[i] + i + 1);
            }
        }
    }
}
