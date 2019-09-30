#ifndef AMREX_NEIGHBORPARTICLESCPUIMPL_H_
#define AMREX_NEIGHBORPARTICLESCPUIMPL_H_

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::fillNeighborsCPU () {
    BL_PROFILE("NeighborParticleContainer::fillNeighborsCPU");
    BuildMasks();
    GetNeighborCommTags();
    cacheNeighborInfo();
    updateNeighborsCPU(false);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::sumNeighborsCPU (int real_start_comp, int real_num_comp,
                   int int_start_comp,  int int_num_comp)
{
    BL_PROFILE("NeighborParticleContainer::sumNeighborsCPU");
    
    if ( not enableInverse() )
    {
        amrex::Abort("Need to enable inverse to true to use sumNeighbors. \n");
    }
    
    const int MyProc = ParallelDescriptor::MyProc();

    std::map<int, Vector<char> > isend_data;
    
    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            PairIndex src_index(pti.index(), pti.LocalTileIndex());
            const auto& tags = inverse_tags[lev][src_index];
            const auto& neighbs = neighbors[lev][src_index];            
            AMREX_ASSERT(tags.size() == neighbs.size());
            
            const int num_neighbs = neighbs.size();
            for (int i = 0; i < num_neighbs; ++i)
            {
                const auto& neighb = neighbs[i];
                const auto& tag = tags[i];
                const int dst_grid = tag.src_grid;
                const int dst_proc = this->ParticleDistributionMap(lev)[dst_grid];
                const int dst_tile = tag.src_tile;
                const int dst_index = tag.src_index;
                const int dst_level = tag.src_level;
                
                if (dst_proc == MyProc)
                {
                    auto pair = std::make_pair(dst_grid, dst_tile);
                    auto& dst_ptile = this->GetParticles(dst_level)[pair];
                    auto& dst_parts = dst_ptile.GetArrayOfStructs();
                    auto& p = dst_parts[dst_index];
                    
                    for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
                    {
                        p.rdata(comp) += neighb.rdata(comp);                        
                    }

                    for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
                    {
                        p.idata(comp) += neighb.idata(comp);                        
                    }
                }
                
                else
                {
                    auto& sdata = isend_data[dst_proc];
                    auto old_size = sdata.size();
                    auto new_size = old_size + real_num_comp*sizeof(Real) + int_num_comp*sizeof(int) + 4*sizeof(int);
                    sdata.resize(new_size);
                    char* dst = &sdata[old_size];
                    std::memcpy(dst, &dst_grid, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_tile, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_index, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_level, sizeof(int)); dst += sizeof(int);
                    for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
                    {
                        Real data = neighb.rdata(comp);
                        std::memcpy(dst, &data, sizeof(Real));
                        dst += sizeof(Real);
                    }
                    for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
                    {
                        int data = neighb.idata(comp);
                        std::memcpy(dst, &data, sizeof(int));
                        dst += sizeof(int);
                    }
                }
            }                    
        }
    }

    sumNeighborsMPI(isend_data, real_start_comp, real_num_comp, int_start_comp, int_num_comp);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
sumNeighborsMPI (std::map<int, Vector<char> >& not_ours,
                 int real_start_comp, int real_num_comp,
                 int int_start_comp, int int_num_comp) 
{
    BL_PROFILE("NeighborParticleContainer::sumNeighborsMPI");
    
#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    
    Vector<long> isnds(NProcs, 0);
    Vector<long> ircvs(NProcs, 0);
    for (int i = 0; i < NProcs; ++i)
        ircvs[i] = 0;
    
    {
        // each proc figures out how many bytes it will send, and how
        // many it will receive
        
        long num_isnds = 0;
        for (const auto& kv : not_ours)
        {
            num_isnds      += kv.second.size();
            isnds[kv.first] = kv.second.size();
        }
        ParallelDescriptor::ReduceLongMax(num_isnds);
        
        if (num_isnds == 0) return;
        
        const int num_ircvs = neighbor_procs.size();
        Vector<MPI_Status>  stats(num_ircvs);
        Vector<MPI_Request> rreqs(num_ircvs);
        
        const int SeqNum = ParallelDescriptor::SeqNum();
        
        // Post receives
        for (int i = 0; i < num_ircvs; ++i)
        {
            const int Who = neighbor_procs[i];
            const long Cnt = 1;
            
            AMREX_ASSERT(Who >= 0 && Who < NProcs);
            
            rreqs[i] = ParallelDescriptor::Arecv(&ircvs[Who], Cnt, Who, SeqNum).req();
        }
        
        // Send.
        for (int i = 0; i < num_ircvs; ++i) {
        const int Who = neighbor_procs[i];
        const long Cnt = 1;
        
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        
        ParallelDescriptor::Send(&isnds[Who], Cnt, Who, SeqNum);
        }
        
        if (num_ircvs > 0) ParallelDescriptor::Waitall(rreqs, stats);        
    }
    
    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (ircvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += ircvs[i];
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
        const auto Cnt    = ircvs[Who];

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }

    // Send.
    for (const auto& kv : not_ours) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0)
    {
        ParallelDescriptor::Waitall(rreqs, stats);

        const size_t data_size = real_num_comp*sizeof(Real) + int_num_comp*sizeof(int) + 4 * sizeof(int);
        
        if (recvdata.size() % data_size != 0) {
            amrex::Print() << recvdata.size() << " " << data_size << "\n";
            if (this->m_verbose) {
                amrex::AllPrint() << "NeighborParticles::sumNeighbors: sizes = "
                                  << recvdata.size() << ", " << data_size << "\n";
            }
            amrex::Abort("NeighborParticles::sumNeighbors: How did this happen?");
        }

        int npart = recvdata.size() / data_size;
        
        char* buffer = recvdata.data();
        for (int j = 0; j < npart; ++j)
        {
            int grid, tile, index, lev;
            std::memcpy(&grid,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&tile,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&index,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&lev, buffer, sizeof(int)); buffer += sizeof(int);

            auto pair = std::make_pair(grid, tile);
            auto& ptile = this->GetParticles(lev)[pair];
            auto& parts = ptile.GetArrayOfStructs();
            auto& p = parts[index];

            for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
            {
                Real data;
                std::memcpy(&data, buffer, sizeof(Real));
                p.rdata(comp) += data;
                buffer += sizeof(Real);
            }

            for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
            {
                int data;
                std::memcpy(&data, buffer, sizeof(int));
                p.idata(comp) += data;
                buffer += sizeof(int);
            }
        }
    }
#endif
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::updateNeighborsCPU (bool reuse_rcv_counts) {

    BL_PROFILE_VAR("NeighborParticleContainer::updateNeighborsCPU", update);

    const int MyProc = ParallelDescriptor::MyProc();

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        const Periodicity& periodicity = this->Geom(lev).periodicity();
        const RealBox& prob_domain = this->Geom(lev).ProbDomain();

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
                auto& tags = buffer_tag_cache[lev][src_index][j];
                int num_tags = tags.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (unsigned i = 0; i < num_tags; ++i) {
                    const NeighborCopyTag& tag = tags[i];
                    const int who = this->ParticleDistributionMap(tag.level)[tag.grid];
                    ParticleType p = particles[tag.src_index];  // copy
                    if (periodicity.isAnyPeriodic()) {
                        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                            if (not periodicity.isPeriodic(dir)) continue;
                            if (tag.periodic_shift[dir] < 0)
                                p.pos(dir) += prob_domain.length(dir);
                            else if (tag.periodic_shift[dir] > 0)
                                p.pos(dir) -= prob_domain.length(dir);
                        }
                    }
                    if (who == MyProc) {
                        PairIndex dst_index(tag.grid, tag.tile);
                        ParticleVector& buffer = neighbors[tag.level][dst_index];
                        AMREX_ASSERT(tag.dst_index < buffer.size());
                        std::memcpy(&buffer[tag.dst_index], &p, pdata_size);
                        if ( enableInverse() )
                        {
                            auto& itags = inverse_tags[tag.level][dst_index];
                            AMREX_ASSERT(tag.dst_index < itags.size());
                            itags[tag.dst_index].src_grid = src_index.first;
                            itags[tag.dst_index].src_tile = src_index.second;
                            itags[tag.dst_index].src_index = tag.src_index;
                            itags[tag.dst_index].src_level = lev;
                        }
                    } else {
                        char* dst = &send_data[who][tag.dst_index];
                        char* src = (char *) &p;
                        for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal; ++ii) {
                            if (rc[ii]) {
                                std::memcpy(dst, src, sizeof(typename ParticleType::RealType));
                                dst += sizeof(typename ParticleType::RealType);
                            }
                            src += sizeof(typename ParticleType::RealType);
                        }
                        for (int ii = 0; ii < 2 + NStructInt; ++ii) {
                            if (ic[ii]) {
                                std::memcpy(dst, src, sizeof(int));
                                dst += sizeof(int);
                            }
                            src += sizeof(int);
                        }
                        if ( enableInverse() )
                        {
                            std::memcpy(dst,&(src_index.first),sizeof(int)); dst += sizeof(int);
                            std::memcpy(dst,&(src_index.second),sizeof(int)); dst += sizeof(int);
                            std::memcpy(dst,&(tag.src_index),sizeof(int)); dst += sizeof(int);
                            std::memcpy(dst,&(lev),sizeof(int)); dst += sizeof(int);
                        }
                    }
                }
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const int grid = mfi.index();
            const int tile = mfi.LocalTileIndex();
            PairIndex dst_index(grid, tile);
            neighbors[lev][dst_index].resize(local_neighbor_sizes[lev][dst_index]);
            if ( enableInverse() )
               inverse_tags[lev][dst_index].resize(local_neighbor_sizes[lev][dst_index]);
        }
    }
    BL_PROFILE_VAR_STOP(update);

    fillNeighborsMPI(reuse_rcv_counts);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::clearNeighborsCPU ()
{
    BL_PROFILE("NeighborParticleContainer::clearNeighborsCPU");

    resizeContainers(this->numLevels());
    for (int lev = 0; lev < this->numLevels(); ++lev) {
        neighbors[lev].clear();
        if ( enableInverse() )
           inverse_tags[lev].clear();
        buffer_tag_cache[lev].clear();
    }

    send_data.clear();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
getRcvCountsMPI () {

    BL_PROFILE("NeighborParticleContainer::getRcvCountsMPI");

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Vector<long> snds(NProcs, 0);
    rcvs.resize(NProcs);
    for (int i = 0; i < NProcs; ++i)
        rcvs[i] = 0;

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

        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&rcvs[Who], Cnt, Who, SeqNum).req();
    }

    // Send.
    for (int i = 0; i < num_rcvs; ++i) {
        const int Who = neighbor_procs[i];
        const long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        ParallelDescriptor::Send(&snds[Who], Cnt, Who, SeqNum);
    }

    if (num_rcvs > 0) ParallelDescriptor::Waitall(rreqs, stats);

#endif // BL_USE_MPI
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
fillNeighborsMPI (bool reuse_rcv_counts) {

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

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }

    // Send.
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);
        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, lev, gid, tid, size, np;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&lev,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);

                if (size == 0) continue;

                np = size / cdata_size;

                AMREX_ASSERT(size % cdata_size == 0);

                PairIndex dst_index(gid, tid);
                size_t old_size = neighbors[lev][dst_index].size();
                size_t new_size = neighbors[lev][dst_index].size() + np;
                if ( enableInverse() )
                {
                    AMREX_ASSERT(neighbors[lev][dst_index].size() ==
                                 inverse_tags[lev][dst_index].size());
                    inverse_tags[lev][dst_index].resize(new_size);
                }
                neighbors[lev][dst_index].resize(new_size);
                
                char* dst = (char*) &neighbors[lev][dst_index][old_size];
                char* src = buffer;

                for (int n = 0; n < np; ++n) {
                    for (int ii = 0; ii < AMREX_SPACEDIM + NStructReal; ++ii) {
                        if (rc[ii]) {
                            std::memcpy(dst, src, sizeof(typename ParticleType::RealType));
                            src += sizeof(typename ParticleType::RealType);
                        }
                        dst += sizeof(typename ParticleType::RealType);
                    }
                    for (int ii = 0; ii < 2 + NStructInt; ++ii) {
                        if (ic[ii]) {
                            std::memcpy(dst, src, sizeof(int));
                            src += sizeof(int);
                        }
                        dst += sizeof(int);
                    }
                    dst = ((char*) &neighbors[lev][dst_index][old_size+n]) + pdata_size;

                    if ( enableInverse() )
                    {
                        auto& tag = inverse_tags[lev][dst_index][old_size+n];
                        std::memcpy(&(tag.src_grid),src,sizeof(int));
                        src += sizeof(int);
                        
                        std::memcpy(&(tag.src_tile),src,sizeof(int));
                        src += sizeof(int);
                        
                        std::memcpy(&(tag.src_index),src,sizeof(int));
                        src += sizeof(int);

                        std::memcpy(&(tag.src_level),src,sizeof(int));
                        src += sizeof(int);
                    }
                }
                buffer += size;
            }
        }
    }
#endif
}

template <int NStructReal, int NStructInt>
template <class CheckPair>
void
NeighborParticleContainer<NStructReal, NStructInt>::
buildNeighborListCPU (CheckPair check_pair, bool sort) {

    BL_PROFILE("NeighborParticleContainer::buildNeighborList");
    AMREX_ASSERT(this->OK());

    for (int lev = 0; lev < this->numLevels(); ++lev) {

        neighbor_list[lev].clear();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            neighbor_list[lev][index];
        }

        IntVect ref_fac = computeRefFac(0, lev);

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
#ifdef AMREX_USE_CUDA
            Cuda::HostVector<int> nl;
#else
            IntVector& nl = neighbor_list[lev][index];
#endif
            AoS& particles = pti.GetArrayOfStructs();

            int Np = particles.size();
            int Nn = neighbors[lev][index].size();
            int N = Np + Nn;

            cells.resize(N);
            tmp_particles.resize(N);
            std::memcpy(&tmp_particles[0], particles.data(), Np*sizeof(ParticleType));
            if (Nn > 0)
                std::memcpy(&tmp_particles[Np], neighbors[lev][index].dataPtr(), Nn*pdata_size);

            // For each cell on this tile, we build linked lists storing the
            // indices of the particles belonging to it.
            Box box = pti.tilebox();
            box.coarsen(ref_fac);
            box.grow(m_num_neighbor_cells+1); // need an extra cell to account for roundoff errors.
            head.resize(box);
            head.setVal(-1);
            list.resize(N, -1);

            for (int i = 0; i < N; ++i) {
                const ParticleType& p = tmp_particles[i];
                const IntVect& cell = this->Index(p, 0);  // we always bin on level 0
                cells[i] = cell;
                list[i] = head(cell);
                head(cell) = i;
            }

            // using these lists, we build a neighbor list containing both
            // kinds of particles.
            int p_start_index = 0;
            for (unsigned i = 0; i < Np; ++i) {
                const ParticleType& p = tmp_particles[i];

                int num_neighbors = 0;
                nl.push_back(0);

                const IntVect& cell = cells[i];
                Box bx(cell, cell);
                bx.grow(m_num_neighbor_cells);

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
#ifdef AMREX_USE_CUDA
                    thrust::sort(nl.begin() + i + 1,
                                 nl.begin() + nl[i] + i + 1);
#else
                    std::sort(nl.begin() + i + 1,
                              nl.begin() + nl[i] + i + 1);
#endif
                }
            }
#ifdef AMREX_USE_CUDA
            neighbor_list[lev][index].resize(nl.size());
            thrust::copy(nl.begin(), nl.end(), neighbor_list[lev][index].begin());
#endif
        }
        }
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
printNeighborListCPU (const std::string& prefix)
{
    BL_PROFILE("NeighborParticleContainer::printNeighborListCPU");

    for (int lev = 0; lev < this->numLevels(); ++lev) {

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            for (MyParIter pti(*this, lev, MFItInfo().SetDynamic(true)); pti.isValid(); ++pti) {
                PairIndex index(pti.index(), pti.LocalTileIndex());
#ifdef AMREX_USE_CUDA
                Cuda::HostVector<int> nl;
#else
                IntVector& nl = neighbor_list[lev][index];
#endif
                if (nl.size() == 0) continue;
                std::stringstream ss;
                ss << prefix << "_level_" << lev;
                int ind = 0;
                while (ind < nl.size()) {
                    int num_partners = nl[ind++];
                    amrex::AllPrintToFile(ss.str()) << num_partners << ": \n ";
                    amrex::AllPrintToFile(ss.str()) << "\t";
                    for (int i = ind; i < ind + num_partners; ++i) {
                        amrex::AllPrintToFile(ss.str()) << nl[i] << " ";
                    }
                    amrex::AllPrintToFile(ss.str()) << "\n";
                    ind += num_partners;
                }
            }
        }
    }
}

#endif
