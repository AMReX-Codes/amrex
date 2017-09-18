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
        std::cout << "actually rebuilding " << std::endl;
        mask_ptr.reset(new FabArray<BaseFab<int> > (ba, dmap, 2, num_neighbor_cells));
        mask_ptr->setVal(-1, num_neighbor_cells);
        for (MFIter mfi(*mask_ptr, this->tile_size); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            mask_ptr->setVal(grid_id, box, 0, 1);
            mask_ptr->setVal(tile_id, box, 1, 1);
        }
        
        mask_ptr->FillBoundary(geom.periodicity());
    };
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::fillNeighbors(int lev) {

    BL_PROFILE("NeighborParticleContainer::fillNeighbors");
    BL_ASSERT(lev == 0);
    NeighborCommMap neighbors_to_comm;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int& grid = pti.index();
        const int& tile = pti.LocalTileIndex(); 
        PairIndex src_index(grid, tile);

        Box shrink_box = pti.tilebox();
        shrink_box.grow(-num_neighbor_cells);
        
        auto& particles = pti.GetArrayOfStructs();
        for (unsigned i = 0; i < pti.numParticles(); ++i) {
            const ParticleType& p = particles[i];
            const IntVect& iv = this->Index(p, lev);
            const int nc = num_neighbor_cells;
            
            // if the particle is more than one cell away from 
            // the tile boundary, it's not anybody's neighbor
            if (shrink_box.contains(iv)) continue;
            
            // Figure out all our neighbors, removing duplicates
            std::vector<NeighborTag> tags;
            for (int ii = -nc; ii < nc + 1; ii += nc) {
                for (int jj = -nc; jj < nc + 1; jj += nc) {
                    for (int kk = -nc; kk < nc + 1; kk += nc) {
                        if ((ii == 0) and (jj == 0) and (kk == 0)) continue;
                        IntVect shift(AMREX_D_DECL(ii, jj, kk));
                        IntVect neighbor_cell = iv + shift;
                        const NeighborTag& tag = getNeighborTag(lev, neighbor_cell, (*mask_ptr)[pti]);
                        if ((grid == tag.grid)           and 
                            (tile == tag.tile)           and 
                            (tag.periodic_shift[0] == 0) and 
                            (tag.periodic_shift[1] == 0) and
                            (tag.periodic_shift[2] == 0)) continue;
                        tags.push_back(tag);
                    }
                }
            }
            std::sort(tags.begin(), tags.end());
            auto it = std::unique(tags.begin(), tags.end());
            tags.erase(it, tags.end());

            // Add neighbors to buffers
            for (int j = 0; j < static_cast<int>(tags.size()); ++j) {
                cacheNeighborParticle(lev, src_index, i, tags[j]);
                packNeighborParticle(lev, p, tags[j], neighbors_to_comm);
            }
        }
    }
    
    fillNeighborsMPI(neighbors_to_comm);
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::updateNeighbors(int lev) {

    BL_PROFILE("NeighborParticleContainer::updateNeighbors");
    BL_ASSERT(lev == 0);

    neighbors.clear();
    
    NeighborCommMap neighbors_to_comm;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        PairIndex src_index(pti.index(), pti.LocalTileIndex());
        auto& particles = pti.GetArrayOfStructs();
        for (unsigned i = 0; i < buffer_id_cache[src_index].size(); ++i) {
            int particle_index = buffer_id_cache[src_index][i];
            const NeighborTag& tag = buffer_tag_cache[src_index][i];
            const ParticleType& p = particles[particle_index];
            packNeighborParticle(lev, p, tag, neighbors_to_comm);
        }
    }
    
    bool reuse_rcv_counts = true;
    fillNeighborsMPI(neighbors_to_comm, reuse_rcv_counts);
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>
::clearNeighbors(int lev) 
{

    BL_PROFILE("NeighborParticleContainer::clearNeighbors");
    BL_ASSERT(lev == 0);

    neighbors.clear();
    buffer_id_cache.clear();
    buffer_tag_cache.clear();
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
applyPeriodicShift(int lev, ParticleType& p, const NeighborTag& tag) {
    
    BL_ASSERT(lev == 0);
    
    const Periodicity& periodicity = this->Geom(lev).periodicity();
    if (not periodicity.isAnyPeriodic()) return;
    
    const RealBox& prob_domain = this->Geom(lev).ProbDomain();
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        if (not periodicity.isPeriodic(dir)) continue;
        if      (tag.periodic_shift[dir] == -1) p.pos(dir) += prob_domain.length(dir);
        else if (tag.periodic_shift[dir] ==  1) p.pos(dir) -= prob_domain.length(dir);
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
typename NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::NeighborTag 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
getNeighborTag(int lev,
               const IntVect& neighbor_cell,
               const BaseFab<int>& mask) {

    BL_ASSERT(lev == 0);
    
    NeighborTag tag;
    tag.grid = mask(neighbor_cell, 0);
    tag.tile = mask(neighbor_cell, 1);
    for (int dim = 0; dim < 3; ++dim)
        tag.periodic_shift[dim] = 0;

    const Periodicity& periodicity = this->Geom(lev).periodicity();
    if (not periodicity.isAnyPeriodic()) return tag;

    const Box& domain = this->Geom(lev).Domain();
    const IntVect& lo = domain.smallEnd();
    const IntVect& hi = domain.bigEnd();
    
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        if (not periodicity.isPeriodic(dir)) continue;
        if      (neighbor_cell[dir] < lo[dir]) tag.periodic_shift[dir] = -1;
        else if (neighbor_cell[dir] > hi[dir]) tag.periodic_shift[dir] =  1;
    }

    return tag;
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
packNeighborParticle(int lev,
                     const ParticleType& p,
                     const NeighborTag& tag,
                     NeighborCommMap& neighbors_to_comm) {
    
    BL_ASSERT(lev == 0);
    
    if (tag.grid >= 0) {
        const int who = this->ParticleDistributionMap(lev)[tag.grid];
        const int MyProc = ParallelDescriptor::MyProc();
        PairIndex dst_index(tag.grid, tag.tile);
        ParticleType particle = p;
        applyPeriodicShift(lev, particle, tag);
        if (who == MyProc) {
            size_t old_size = neighbors[dst_index].size();
            size_t new_size = neighbors[dst_index].size() + pdata_size;
            neighbors[dst_index].resize(new_size);
            std::memcpy(&neighbors[dst_index][old_size], &particle, pdata_size);
        } else {
            NeighborCommTag comm_tag(who, tag.grid, tag.tile);
            Array<char>& buffer = neighbors_to_comm[comm_tag];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + pdata_size;
            buffer.resize(new_size);
            std::memcpy(&buffer[old_size], &particle, pdata_size);
        }
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void 
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
cacheNeighborParticle(int lev, const PairIndex& index,
                      int particle_index, const NeighborTag& tag) {
    
    BL_ASSERT(lev == 0);
    
    if (tag.grid >= 0) {
        PairIndex dst_index(tag.grid, tag.tile);
        buffer_id_cache[index].push_back(particle_index);
        buffer_tag_cache[index].push_back(tag);
    }
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
getRcvCountsMPI(const std::map<int, Array<char> >& send_data) {

#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Array<long> snds(NProcs, 0);
    rcvs.resize(NProcs, 0);

    num_snds = 0;
    for (const auto& kv : send_data) {
        num_snds      += kv.second.size();
        snds[kv.first] = kv.second.size();
    }
    ParallelDescriptor::ReduceLongMax(num_snds);
    if (num_snds == 0) return;
    
    // communicate that information
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(rcvs[MyProc] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

#endif // BL_USE_MPI
}

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
fillNeighborsMPI(NeighborCommMap& neighbors_to_comm, bool reuse_rcv_counts) {
    
    BL_PROFILE("NeighborParticleContainer::fillNeighborsMPI");
    
#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    
    // count the number of tiles to be sent to each proc
    std::map<int, int> tile_counts;
    for (const auto& kv: neighbors_to_comm) {
        tile_counts[kv.first.proc_id] += 1;
    }
    
    // flatten all the data for each proc into a single buffer
    // once this is done, each dst proc will have an Array<char>
    // the buffer will be packed like:
    // ntiles, gid1, tid1, size1, data1....  gid2, tid2, size2, data2... etc. 
    std::map<int, Array<char> > send_data;
    for (const auto& kv: neighbors_to_comm) {
        Array<char>& buffer = send_data[kv.first.proc_id];
        buffer.resize(sizeof(int));
        std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
    }
    
    for (auto& kv : neighbors_to_comm) {
        int data_size = kv.second.size();
        Array<char>& buffer = send_data[kv.first.proc_id];
        size_t old_size = buffer.size();
        size_t new_size = buffer.size() + 2*sizeof(int) + sizeof(int) + data_size;
        buffer.resize(new_size);
        char* dst = &buffer[old_size];
        std::memcpy(dst, &(kv.first.grid_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.tile_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &data_size,          sizeof(int)); dst += sizeof(int);
        if (data_size == 0) continue;
        std::memcpy(dst, &kv.second[0], data_size);
        Array<char>().swap(kv.second);
    }
    
    // each proc figures out how many bytes it will send, and how
    // many it will receive
    if (!reuse_rcv_counts) getRcvCountsMPI(send_data);
    if (num_snds == 0) return;
    
    Array<int> RcvProc;
    Array<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Array<MPI_Status>  stats(nrcvs);
    Array<MPI_Request> rreqs(nrcvs);
    
    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.
    Array<char> recvdata(TotRcvBytes);
    
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
    
    BL_PROFILE("NeighborListParticleContainer::buildNeighborList");
    BL_ASSERT(lev == 0);

    neighbor_list.clear();
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        neighbor_list[index] = Array<int>();
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        
        PairIndex index(pti.index(), pti.LocalTileIndex());
        Array<int>& nl = neighbor_list[index];
        AoS& particles = pti.GetArrayOfStructs();

        int Np = particles.size();
        int Nn = neighbors[index].size() / pdata_size;
        int N = Np + Nn;

        Array<IntVect> cells(N);
        Array<ParticleType> tmp_particles(N);
        std::memcpy(&tmp_particles[0], particles.data(), Np*sizeof(ParticleType));

        for (int i = 0; i < Nn; ++i) {
            std::memcpy(&tmp_particles[i + Np],
                        neighbors[index].dataPtr() + i*pdata_size,
                        pdata_size);
        }

        // For each cell on this tile, we build linked lists storing the
        // indices of the particles belonging to it.
        Box box = pti.tilebox();
        box.grow(num_neighbor_cells);
        BaseFab<int> head(box);
        head.setVal(-1);
        Array<int> list(N, -1);

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
                        if ((p.pos(0) - tmp_particles[j].pos(0))*(p.pos(0) - tmp_particles[j].pos(0)) + 
                            (p.pos(1) - tmp_particles[j].pos(1))*(p.pos(1) - tmp_particles[j].pos(1)) +
                            (p.pos(2) - tmp_particles[j].pos(2))*(p.pos(2) - tmp_particles[j].pos(2)) == 0.0) {
                            std::cout << p << " " << tmp_particles[j] << " " << i << " " << j << std::endl; 
                        }
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

template <int NStructReal, int NStructInt, int NNeighborReal>
void
NeighborParticleContainer<NStructReal, NStructInt, NNeighborReal>::
buildNeighborListFort(int lev, bool sort) {
    
    BL_PROFILE("NeighborListParticleContainer::buildNeighborListFort");
    BL_ASSERT(lev == 0);

    neighbor_list.clear();

    const Geometry& gm  = this->Geom(lev);
    const Real*     plo = gm.ProbLo();
    const Real*     dx  = gm.CellSize();
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        
        PairIndex index(pti.index(), pti.LocalTileIndex());
        Array<int>& nl = neighbor_list[index];
        AoS& particles = pti.GetArrayOfStructs();

        int Np = particles.size();
        int Nn = neighbors[index].size() / pdata_size;
        int Ns = particles.dataShape().first;
        int N = Np + Nn;

        Array<ParticleType> tmp_particles(N);
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
