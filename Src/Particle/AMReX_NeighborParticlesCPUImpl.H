#ifndef AMREX_NEIGHBORPARTICLESCPUIMPL_H_
#define AMREX_NEIGHBORPARTICLESCPUIMPL_H_

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::fillNeighborsCPU() {
    BL_PROFILE("NeighborParticleContainer::fillNeighborsCPU");
    BuildMasks();
    GetNeighborCommTags();
    cacheNeighborInfo();
    updateNeighborsCPU(false);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::updateNeighborsCPU(bool reuse_rcv_counts) {

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
                    ParticleType p = particles[tag.src_index];
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
        }
    }
    BL_PROFILE_VAR_STOP(update);

    fillNeighborsMPI(reuse_rcv_counts);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>
::clearNeighborsCPU()
{
    BL_PROFILE("NeighborParticleContainer::clearNeighborsCPU");

    resizeContainers(this->numLevels());
    for (int lev = 0; lev < this->numLevels(); ++lev) {
        neighbors[lev].clear();
        buffer_tag_cache[lev].clear();
    }

    send_data.clear();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
getRcvCountsMPI() {

    BL_PROFILE("NeighborParticleContainer::getRcvCountsMPI");

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();

    AMREX_ASSERT(send_data.size() == neighbor_procs.size());

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

                BL_ASSERT(size % cdata_size == 0);

                PairIndex dst_index(gid, tid);
                size_t old_size = neighbors[lev][dst_index].size();
                size_t new_size = neighbors[lev][dst_index].size() + np;
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
buildNeighborListCPU(CheckPair check_pair, bool sort) {

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
printNeighborListCPU(const std::string& prefix)
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
