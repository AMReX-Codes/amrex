#ifdef AMREX_USE_CUDA

#include <thrust/device_vector.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_ptr.h>
#include <thrust/tuple.h>

#include <AMReX_RedistributeStrategy.H>

namespace {
    
    struct DeviceBox
    {
        int lo[AMREX_SPACEDIM];
        int hi[AMREX_SPACEDIM];
        
        DeviceBox() 
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) 
            {
                lo[i] = 0;
                hi[i] = 0;
            }
        }
        
        DeviceBox(const Box& a_box)
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            {
                lo[i] = a_box.smallEnd(i);
                hi[i] = a_box.bigEnd(i);
            }
        }
    };

    struct DeviceDomain
    {
        Real left_edge[AMREX_SPACEDIM];
        Real inverse_dx[AMREX_SPACEDIM];
        
        DeviceDomain(const Geometry& geom) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                left_edge[i] = geom.ProbLo(i);            
                inverse_dx[i] = geom.InvCellSize(i);

            }
        }
    };

    struct assignParticle
    {        
        DeviceDomain domain;
        DeviceBox box;
        const Real* mask_ptr;
        
        __host__ __device__
        assignParticle(const DeviceDomain& a_domain,
                       const DeviceBox&    a_box,
                       const Real*         a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Tuple>
        __host__ __device__
        int operator()(Tuple tup) const {
            int i = floor((thrust::get<0>(tup) - domain.left_edge[0]) * domain.inverse_dx[0]);
            int j = floor((thrust::get<1>(tup) - domain.left_edge[1]) * domain.inverse_dx[1]);
            int k = floor((thrust::get<2>(tup) - domain.left_edge[2]) * domain.inverse_dx[2]);
            
            long offset = i - box.lo[0];
#if   AMREX_SPACEDIM==2
            offset += (box.hi[0] - box.lo[0] + 1) * (j - box.lo[1]);
#elif AMREX_SPACEDIM==3
            offset += (box.hi[0] - box.lo[0] + 1) * (j - box.lo[1]
                                                     + (k - box.lo[2]) * (box.hi[1] - box.lo[1] + 1));
#endif
            return mask_ptr[offset];
        }
    };
    
    struct grid_is
    {
        int grid_id;
        
        __host__ __device__
        grid_is(int a_grid_id) : grid_id(a_grid_id) {}
        
        __host__ __device__
        bool operator()(int grid) const
        {
            return grid_id == grid;
        }
    };
    
    struct grid_is_not
    {
        int grid_id;
        
        __host__ __device__
        grid_is_not(int a_grid_id) : grid_id(a_grid_id) {}
    
        __host__ __device__
        bool operator()(int grid) const
        {
            return grid_id != grid;
        }
    };
}

void RedistributeStrategyGPU::Redistribute (std::map<std::pair<int, int>, Particles>&         a_particles,
                                            const amrex::BoxArray&            a_ba,
                                            const amrex::DistributionMapping& a_dm,
                                            const amrex::Geometry&            a_geom,
                                            const amrex::MultiFab*            a_mask_ptr)
{
    BL_PROFILE("RedistributeStrategyGPU::Redistribute");
    
    amrex::StructOfArrays<PIdx::nattribs, 0> particles_to_redistribute;
    
    thrust::device_vector<int> grids;
    thrust::device_vector<int> grids_copy;    
    thrust::device_vector<int> grids_to_redistribute;
    
    for(MFIter mfi(*a_mask_ptr); mfi.isValid(); ++mfi) {
        int src_grid = mfi.index();
        const int tid = 0;
        
        auto& attribs = a_particles[std::make_pair(src_grid, tid)].attribs;
        const size_t old_num_particles = attribs.numParticles();
        
        grids.resize(old_num_particles);
        grids_copy.resize(old_num_particles);
        
        thrust::transform(thrust::device, attribs.begin(), attribs.end(),
                          grids.begin(),
                          assignParticle(DeviceDomain(a_geom),
                                         DeviceBox((*a_mask_ptr)[src_grid].box()),
                                         (*a_mask_ptr)[src_grid].dataPtr()));
                
        thrust::copy(grids.begin(), grids.end(), grids_copy.begin());

        auto begin = thrust::make_zip_iterator(thrust::make_tuple(attribs.begin(), grids_copy.begin()));
        auto end   = thrust::make_zip_iterator(thrust::make_tuple(attribs.end(),   grids_copy.end()));
        auto mid   = thrust::partition(begin, end, grids.begin(), grid_is(src_grid));
                
        const size_t num_to_redistribute = thrust::distance(mid, end);
        const size_t new_num_particles   = old_num_particles - num_to_redistribute;
        
        const size_t old_size = particles_to_redistribute.numParticles();
        const size_t new_size = old_size + num_to_redistribute;
        
        particles_to_redistribute.resize(new_size);
        grids_to_redistribute.resize(new_size);
        
        auto pos = thrust::make_zip_iterator(thrust::make_tuple(particles_to_redistribute.begin() + old_size, 
                                                                grids_to_redistribute.begin()     + old_size));
        thrust::copy(mid, end, pos);
        
        attribs.resize(new_num_particles);
    }

    const int num_grids = a_ba.size();
    thrust::device_vector<int> grid_begin(num_grids);
    thrust::device_vector<int> grid_end(num_grids);

    thrust::sort_by_key(grids_to_redistribute.begin(),
                        grids_to_redistribute.end(),
                        particles_to_redistribute.begin());
    
    thrust::counting_iterator<unsigned> search_begin(0);
    thrust::lower_bound(grids_to_redistribute.begin(),
                        grids_to_redistribute.end(),
                        search_begin,
                        search_begin + num_grids,
                        grid_begin.begin());
    
    thrust::upper_bound(grids_to_redistribute.begin(),
                        grids_to_redistribute.end(),
                        search_begin,
                        search_begin + num_grids,
                        grid_end.begin());

    thrust::host_vector<int> start(grid_begin);
    thrust::host_vector<int> stop(grid_end);

    thrust::host_vector<int> proc_map(num_grids);
    for (int i = 0; i < num_grids; ++i) proc_map[i] = m_dmap[i];
 
    std::map<int, thrust::device_vector<char> > not_ours;

    std::map<int, size_t> grid_counts;
    for (int i = 0; i < num_grids; ++i)
    {
        const int dest_proc = proc_map[i];
        if (dest_proc != ParallelDescriptor::MyProc())
        {
            grid_counts[dest_proc] += 1;
        }
    }

    for (int i = 0; i < num_grids; ++i)
    {
        auto begin = particles_to_redistribute.begin();
        thrust::advance(begin, start[i]);
        
        auto end = particles_to_redistribute.begin();
        thrust::advance(end, stop[i]);
        
        const size_t num_to_add = stop[i] - start[i];
        const int dest_proc = proc_map[i];

        if (dest_proc == ParallelDescriptor::MyProc())
        {
            const size_t old_size = m_particles[i].attribs.numParticles();
            const size_t new_size = old_size + num_to_add;
            m_particles[i].attribs.resize(new_size);
            thrust::copy(begin, end, m_particles[i].attribs.begin() + old_size);
            m_particles[i].temp.resize(m_particles[i].attribs.numParticles());
        }
        else
        {           
            char* dst;
            const size_t old_size = not_ours[dest_proc].size();
            const size_t new_size 
                = old_size + num_to_add*m_superparticle_size + sizeof(size_t) + 2*sizeof(int);
            if (old_size == 0) {
                not_ours[dest_proc].resize(new_size + sizeof(size_t));
                cudaMemcpy(thrust::raw_pointer_cast(not_ours[dest_proc].data()), 
                           &grid_counts[dest_proc], sizeof(size_t), cudaMemcpyHostToDevice);
                dst = thrust::raw_pointer_cast(not_ours[dest_proc].data() + old_size + sizeof(size_t));
            } else {
                not_ours[dest_proc].resize(new_size);
                dst = thrust::raw_pointer_cast(not_ours[dest_proc].data() + old_size);
            } 

            cudaMemcpy(thrust::raw_pointer_cast(dst), &num_to_add, sizeof(size_t), cudaMemcpyHostToDevice);
            dst += sizeof(size_t);

            cudaMemcpy(thrust::raw_pointer_cast(dst), &i, sizeof(int), cudaMemcpyHostToDevice);
            dst += sizeof(int);

            cudaMemcpy(thrust::raw_pointer_cast(dst), &dest_proc, sizeof(int), cudaMemcpyHostToDevice);
            dst += sizeof(int);

            for (int j = 0; j < PIdx::nattribs; ++j)
            {
                auto& attrib = particles_to_redistribute.GetRealData(j);
                cudaMemcpy(thrust::raw_pointer_cast(dst), attrib.data() + start[i],
                           num_to_add*sizeof(Real), cudaMemcpyDeviceToDevice);
                dst += num_to_add*sizeof(Real);
            }
        }
    }
    
    RedistributeMPI(not_ours);
}

void RedistributeStrategyGPU::OK (std::map<std::pair<int, int>, Particles>&         a_particles,
                                  const amrex::BoxArray&            a_ba,
                                  const amrex::DistributionMapping& a_dm,
                                  const amrex::Geometry&            a_geom,
                                  const amrex::MultiFab*            a_mask_ptr)
{
    BL_PROFILE("RedistributeStrategyGPU::OK");

    thrust::device_vector<int> grid_indices;

    long total_np = 0;
    for(MFIter mfi(*a_mask_ptr); mfi.isValid(); ++mfi) {
        int i = mfi.index();
        const int tid = 0;
        
        auto& particles = a_particles[std::make_pair(i, tid)];
        const int np = particles.numParticles();
        auto& soa = particles.attribs;

        total_np += np;

        if (np == 0) continue;
        
        grid_indices.resize(np);
        
        thrust::transform(soa.begin(),
                          soa.end(),
                          grid_indices.begin(),
                          assignParticle(DeviceDomain(a_geom),
                                         DeviceBox((*a_mask_ptr)[i].box()),
                                         (*a_mask_ptr)[i].dataPtr()));

        int count = thrust::count_if(grid_indices.begin(),
                                     grid_indices.end(),
                                     grid_is_not(i));

        AMREX_ALWAYS_ASSERT(count == 0);
    }

    ParallelDescriptor::ReduceLongMax(total_np);
    amrex::Print() << "I have " << total_np << " particles." << std::endl;
}

void RedistributeStrategyGPU::RedistributeMPI (std::map<int, thrust::device_vector<char> >& not_ours,
                                               std::map<std::pair<int, int>, Particles>& a_particles)
{
    BL_PROFILE("RedistributeStrategyGPU::RedistributeMPI()");
    
#if BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;

    for (const auto& kv : not_ours)
    {
        const int np = kv.second.size(); 
        Snds[kv.first] = np;
        NumSnds += np;
    }
    
    ParallelDescriptor::ReduceLongMax(NumSnds);

    if (NumSnds == 0) return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    
    BL_ASSERT(Rcvs[ParallelDescriptor::MyProc()] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.
    thrust::device_vector<char> recvdata(TotRcvBytes);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(thrust::raw_pointer_cast(recvdata.data() + offset),
                                             Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : not_ours) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(thrust::raw_pointer_cast(kv.second.data()),
                                 Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = thrust::raw_pointer_cast(recvdata.data() + offset);
            size_t num_grids, num_particles;
            int gid, pid;
            cudaMemcpy(&num_grids, buffer, sizeof(size_t), cudaMemcpyDeviceToHost); buffer += sizeof(size_t);
            for (int g = 0; g < num_grids; ++g) {
                cudaMemcpy(&num_particles, buffer, sizeof(size_t), cudaMemcpyDeviceToHost); buffer += sizeof(size_t);
                cudaMemcpy(&gid, buffer, sizeof(int), cudaMemcpyDeviceToHost); buffer += sizeof(int);
                cudaMemcpy(&pid, buffer, sizeof(int), cudaMemcpyDeviceToHost); buffer += sizeof(int);
                
                if (num_particles == 0) continue;

                amrex::StructOfArrays<PIdx::nattribs, 0> redistributed_particles;

                AMREX_ALWAYS_ASSERT(pid == ParallelDescriptor::MyProc());
                {
                    const size_t old_size = redistributed_particles.numParticles();
                    const size_t new_size = old_size + num_particles;        
                    redistributed_particles.resize(new_size);
                    
                    for (int j = 0; j < PIdx::nattribs; ++j) {
                        auto& attrib = redistributed_particles.GetRealData(j);
                        cudaMemcpy(attrib.data() + old_size, buffer, num_particles*sizeof(Real), cudaMemcpyHostToDevice);
                        buffer += num_particles*sizeof(Real);
                    }
                }
            
                {
                    const int tid = 0;
                    auto pair_index = std::make_pair(gid, tid);
                    const size_t old_size = a_particles[pair_index].attribs.numParticles();
                    const size_t new_size = old_size + num_particles;
                    a_particles[pair_index].attribs.resize(new_size);
                    thrust::copy(redistributed_particles.begin(),
                                 redistributed_particles.end(),
                                 a_particles[pair_index].attribs.begin() + old_size);
                    a_particles[pair_index].temp.resize(a_particles[pair_index].attribs.numParticles());
                }
            }
        }
    }
        
#endif // MPI    
}
#endif // CUDA
