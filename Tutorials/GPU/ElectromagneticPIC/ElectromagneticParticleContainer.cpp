#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_CudaLaunch.H>

#include <thrust/device_vector.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_ptr.h>
#include <thrust/tuple.h>

#include "ElectromagneticParticleContainer.H"
#include "Constants.H"

#include "em_pic_F.H"

using namespace amrex;

namespace {

    void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
        int ny = nppc[1];
        int nz = nppc[2];
        
        int ix_part = i_part/(ny * nz);
        int iy_part = (i_part % (ny * nz)) % ny;
        int iz_part = (i_part % (ny * nz)) / ny;
        
        r[0] = (0.5+ix_part)/nx;
        r[1] = (0.5+iy_part)/ny;
        r[2] = (0.5+iz_part)/nz;    
    }
    
    void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
        Real ux_th = amrex::RandomNormal(0.0, u_std);
        Real uy_th = amrex::RandomNormal(0.0, u_std);
        Real uz_th = amrex::RandomNormal(0.0, u_std);
        
        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }


    struct assignParticle
    {        
        GeometryData domain;
        Box box;
        BaseFab<int>* mask_ptr;
        
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE 
        assignParticle(const GeometryData& a_domain,
                       const Box&          a_box,
                       BaseFab<int>*       a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Tuple>
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE 
        int operator()(Tuple tup) const {

            IntVect offset = IntVect(floor((thrust::get<0>(tup) - domain.ProbLo()[0]) / domain.CellSize()[0]),
                                     floor((thrust::get<1>(tup) - domain.ProbLo()[1]) / domain.CellSize()[1]),
                                     floor((thrust::get<2>(tup) - domain.ProbLo()[2]) / domain.CellSize()[2]));

            return (*mask_ptr)(offset);
        }
    };
    
    struct grid_is
    {
        int grid_id;
        
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        grid_is(int a_grid_id) : grid_id(a_grid_id) {}
        
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        bool operator()(int grid) const
        {
            return grid_id == grid;
        }
    };
    
    struct grid_is_not
    {
        int grid_id;
        
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        grid_is_not(int a_grid_id) : grid_id(a_grid_id) {}
    
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        bool operator()(int grid) const
        {
            return grid_id != grid;
        }
    };
}

ElectromagneticParticleContainer::
ElectromagneticParticleContainer(const Geometry            & a_geom,
                                 const DistributionMapping & a_dmap,
                                 const BoxArray            & a_ba,
                                 const int                   a_species_id,
                                 const Real                  a_charge,
                                 const Real                  a_mass)
    : m_ba(a_ba), m_geom(a_geom), m_dmap(a_dmap),
      m_species_id(a_species_id), m_charge(a_charge), m_mass(a_mass)
{
    BL_PROFILE("ElectromagneticParticleContainer::ElectromagneticParticleContainer");
    
    const int ng = 1;
    m_mask_ptr.reset(new iMultiFab(m_ba, m_dmap, 1, ng));
    m_mask_ptr->setVal(-1);
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const int grid_id = mfi.index();
        m_mask_ptr->setVal(grid_id, box, 0, 1);
    }
    m_mask_ptr->FillBoundary(m_geom.periodicity());

    superparticle_size = PIdx::nattribs * sizeof(Real);
}

void
ElectromagneticParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std, 
              const Real     a_thermal_momentum_mean,
              const Real     a_density, 
              const RealBox& a_bounds, 
              const int      a_problem)
{
    BL_PROFILE("ElectromagneticParticleContainer::InitParticles");
    
    const Real* dx = m_geom.CellSize();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0], 
                                     *a_num_particles_per_cell[1], 
                                     *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
    
    for(MFIter mfi(*m_mask_ptr); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const Real* plo = m_geom.ProbLo();
        const int grid_id = mfi.index();
        auto& particles = m_particles[grid_id];
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real u[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                if (a_problem == 0) {
                    get_gaussian_random_momentum(u, a_thermal_momentum_mean, a_thermal_momentum_std);
                }
                else if (a_problem == 1 ) {
                    u[0] = 0.01;
                    u[1] = 0.0;
                    u[2] = 0.0;
                } else {
                    amrex::Abort("problem type not valid");
                }

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;

                particles.x().push_back(x);
                particles.y().push_back(y);
                particles.z().push_back(z);
                
                particles.ux().push_back(u[0] * PhysConst::c);
                particles.uy().push_back(u[1] * PhysConst::c);
                particles.uz().push_back(u[2] * PhysConst::c);
                
                particles.w().push_back(a_density * scale_fac);
                
                particles.ex().push_back(0);
                particles.ey().push_back(0);
                particles.ez().push_back(0);
                particles.bx().push_back(0);
                particles.by().push_back(0);
                particles.bz().push_back(0);
                particles.ginv().push_back(0);
            }
        }
    }
}

void
ElectromagneticParticleContainer::
PushAndDeposeParticles(const amrex::MultiFab& Ex,
                       const amrex::MultiFab& Ey,
                       const amrex::MultiFab& Ez,
                       const amrex::MultiFab& Bx,
                       const amrex::MultiFab& By,
                       const amrex::MultiFab& Bz,
                             amrex::MultiFab& jx, 
                             amrex::MultiFab& jy, 
                             amrex::MultiFab& jz,
                             amrex::Real      dt)
{   
    BL_PROFILE("ElectromagneticParticleContainer::PushAndDeposeParticles");
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        const int np    = m_particles[mfi.index()].size();
        if (np == 0) continue;

        const FArrayBox* elecX = &(Ex[mfi]);
        const FArrayBox* elecY = &(Ey[mfi]);
        const FArrayBox* elecZ = &(Ez[mfi]);
        const FArrayBox* magX = &(Bx[mfi]);
        const FArrayBox* magY = &(By[mfi]);
        const FArrayBox* magZ = &(Bz[mfi]);
        FArrayBox* currDenX = &(jx[mfi]);
        FArrayBox* currDenY = &(jy[mfi]);
        FArrayBox* currDenZ = &(jz[mfi]);
        const GeometryData& geomData = m_geom.data();
        ParticlesData&& pData = m_particles[mfi.index()].data();

        // Create a copy of this class member to pass.
        // this->m_charge shouldn't work, since it is not available on the GPU.
        Real charge = m_charge;
        Real mass = m_mass;

        auto depose = [=] AMREX_CUDA_DEVICE () mutable
        {
            int index, threadSize;
            getThreadIndex(index, threadSize, np);

            if (threadSize == 0) return;

            gather_magnetic_field(threadSize, 
                                  &(pData.x()[index]), &(pData.y()[index]), &(pData.z()[index]),
                                  &(pData.bx()[index]), &(pData.by()[index]), &(pData.bz()[index]),
                                  BL_TO_FORTRAN_3D(*magX),
                                  BL_TO_FORTRAN_3D(*magY),
                                  BL_TO_FORTRAN_3D(*magZ),
                                  geomData.ProbLo(), geomData.CellSize());
        
            gather_electric_field(threadSize, 
                                  &(pData.x()[index]), &(pData.y()[index]), &(pData.z()[index]),
                                  &(pData.ex()[index]), &(pData.ey()[index]), &(pData.ez()[index]),
                                  BL_TO_FORTRAN_3D(*elecX),
                                  BL_TO_FORTRAN_3D(*elecY),
                                  BL_TO_FORTRAN_3D(*elecZ),
                                  geomData.ProbLo(), geomData.CellSize());
        
            push_momentum_boris(threadSize,
                                &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                                &(pData.ginv()[index]),
                                &(pData.ex()[index]), &(pData.ey()[index]), &(pData.ez()[index]),
                                &(pData.bx()[index]), &(pData.by()[index]), &(pData.bz()[index]),
                                charge, mass, dt);
         
            push_position_boris(threadSize,
                                &(pData.x()[index]),  &(pData.y()[index]),  &(pData.z()[index]),
                                &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                                &(pData.ginv()[index]), dt);
        
            deposit_current(BL_TO_FORTRAN_3D(*currDenX),
                            BL_TO_FORTRAN_3D(*currDenY),
                            BL_TO_FORTRAN_3D(*currDenZ),
                            threadSize,
                            &(pData.x()[index]),  &(pData.y()[index]),  &(pData.z()[index]),
                            &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                            &(pData.ginv()[index]), &(pData.w()[index]),
                            charge, geomData.ProbLo(), dt, geomData.CellSize());

        };

        AMREX_CUDA_LAUNCH_DEVICE(Strategy(np), depose); 
    }

}

void
ElectromagneticParticleContainer::
PushParticleMomenta(const amrex::MultiFab& Ex,
                    const amrex::MultiFab& Ey,
                    const amrex::MultiFab& Ez,
                    const amrex::MultiFab& Bx,
                    const amrex::MultiFab& By,
                    const amrex::MultiFab& Bz,
                          amrex::Real      dt)
{   
    BL_PROFILE("ElectromagneticParticleContainer::PushParticleMomenta");
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        const int np    = m_particles[mfi.index()].size();
        if (np == 0) continue;

        const FArrayBox* elecX = &(Ex[mfi]);
        const FArrayBox* elecY = &(Ey[mfi]);
        const FArrayBox* elecZ = &(Ez[mfi]);
        const FArrayBox* magX = &(Bx[mfi]);
        const FArrayBox* magY = &(By[mfi]);
        const FArrayBox* magZ = &(Bz[mfi]);
        const GeometryData& geomData = m_geom.data();
        ParticlesData&& pData = m_particles[mfi.index()].data(); 

        Real charge = m_charge;
        Real mass = m_mass;

        auto pushMom = [=] AMREX_CUDA_DEVICE () mutable
        {
            int index, threadSize;
            getThreadIndex(index, threadSize, np);

            if (threadSize == 0) return;

            gather_magnetic_field(threadSize, 
                                  &(pData.x()[index]),  &(pData.y()[index]),  &(pData.z()[index]),
                                  &(pData.bx()[index]), &(pData.by()[index]), &(pData.bz()[index]),
                                  BL_TO_FORTRAN_3D(*magX),
                                  BL_TO_FORTRAN_3D(*magY),
                                  BL_TO_FORTRAN_3D(*magZ),
                                  geomData.ProbLo(), geomData.CellSize());

            gather_electric_field(threadSize, 
                                  &(pData.x()[index]), &(pData.y()[index]), &(pData.z()[index]),
                                  &(pData.ex()[index]), &(pData.ey()[index]), &(pData.ez()[index]),
                                  BL_TO_FORTRAN_3D(*elecX),
                                  BL_TO_FORTRAN_3D(*elecY),
                                  BL_TO_FORTRAN_3D(*elecZ),
                                  geomData.ProbLo(), geomData.CellSize());
        
            push_momentum_boris(threadSize,
                                &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                                &(pData.ginv()[index]),
                                &(pData.ex()[index]), &(pData.ey()[index]), &(pData.ez()[index]),
                                &(pData.bx()[index]), &(pData.by()[index]), &(pData.bz()[index]),
                                charge, mass, dt);

        };

        AMREX_CUDA_LAUNCH_DEVICE(Strategy(np), pushMom); 
    }
}

void 
ElectromagneticParticleContainer::
PushParticlePositions(amrex::Real dt)
{
    BL_PROFILE("ElectromagneticParticleContainer::PushParticlePositions");
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        const int np    = m_particles[mfi.index()].size();
        if (np == 0) continue;

        ParticlesData&& pData = m_particles[mfi.index()].data();
        auto pushPos = [=] AMREX_CUDA_DEVICE () mutable
        {
            int index, threadSize;
            getThreadIndex(index, threadSize, np);

            if (threadSize == 0) return;

            set_gamma(threadSize, 
                      &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                      &(pData.ginv()[index]));

            push_position_boris(threadSize,
                                &(pData.x()[index]), &(pData.y()[index]), &(pData.z()[index]),
                                &(pData.ux()[index]), &(pData.uy()[index]), &(pData.uz()[index]),
                                &(pData.ginv()[index]), dt);
        };

        AMREX_CUDA_LAUNCH_DEVICE(Strategy(np), pushPos);
    }
}

void
ElectromagneticParticleContainer::
EnforcePeriodicBCs()
{
    BL_PROFILE("ElectromagneticParticleContainer::EnforcePeriodicBCs");

    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        const int np = m_particles[mfi.index()].size();
        if (np == 0) continue;

        ParticlesData&& pData = m_particles[mfi.index()].data(); 
        const GeometryData& geomData = m_geom.data();

        auto periodic = [=] AMREX_CUDA_DEVICE () mutable
        {
            int index, threadSize;
            getThreadIndex(index, threadSize, np);

            if (threadSize == 0) return;

            enforce_periodic(threadSize,
                             &(pData.x()[index]), &(pData.y()[index]), &(pData.z()[index]),
                             geomData.ProbLo(), geomData.ProbLo());
        };

        AMREX_CUDA_LAUNCH_DEVICE(Strategy(np), periodic);
    }
}

void
ElectromagneticParticleContainer::
OK ()
{
    BL_PROFILE("ElectromagneticParticleContainer::OK");

    thrust::device_vector<int> grid_indices;

    long total_np = 0;
    for(MFIter mfi(*m_mask_ptr); mfi.isValid(); ++mfi) {
        int i = mfi.index();
        
        auto& particles = m_particles[i];
        const int np = particles.size();
        auto& soa = particles.attribs;

        total_np += np;

        if (np == 0) continue;
        
        grid_indices.resize(np);
        
        thrust::transform(soa.begin(),
                          soa.end(),
                          grid_indices.begin(),
                          assignParticle( m_geom.data(),
                                          (*m_mask_ptr)[i].box(),
                                          &((*m_mask_ptr)[i]) ));

        int count = thrust::count_if(grid_indices.begin(),
                                     grid_indices.end(),
                                     grid_is_not(i));

        AMREX_ALWAYS_ASSERT(count == 0);
    }

    ParallelDescriptor::ReduceLongMax(total_np);
    amrex::Print() << "I have " << total_np << " particles." << std::endl;
}

void
ElectromagneticParticleContainer::
Redistribute()
{

    BL_PROFILE("ElectromagneticParticleContainer::Redistribute");
    
    StructOfArrays<PIdx::nattribs, 0> particles_to_redistribute;
    
    thrust::device_vector<int> grids;
    thrust::device_vector<int> grids_copy;    
    thrust::device_vector<int> grids_to_redistribute;
    
    for(MFIter mfi(*m_mask_ptr); mfi.isValid(); ++mfi) {
        int src_grid = mfi.index();
        
        auto& attribs = m_particles[src_grid].attribs;
        const size_t old_num_particles = attribs.size();
        
        grids.resize(old_num_particles);
        grids_copy.resize(old_num_particles);
        
        thrust::transform(thrust::device, attribs.begin(), attribs.end(),
                          grids.begin(),
                          assignParticle( m_geom.data(),
                                          (*m_mask_ptr)[src_grid].box(),
                                          &((*m_mask_ptr)[src_grid]) ));

                
        thrust::copy(grids.begin(), grids.end(), grids_copy.begin());

        auto begin = thrust::make_zip_iterator(thrust::make_tuple(attribs.begin(), grids_copy.begin()));
        auto end   = thrust::make_zip_iterator(thrust::make_tuple(attribs.end(),   grids_copy.end()));
        auto mid   = thrust::partition(begin, end, grids.begin(), grid_is(src_grid));
                
        const size_t num_to_redistribute = thrust::distance(mid, end);
        const size_t new_num_particles   = old_num_particles - num_to_redistribute;
        
        const size_t old_size = particles_to_redistribute.size();
        const size_t new_size = old_size + num_to_redistribute;
        
        particles_to_redistribute.resize(new_size);
        grids_to_redistribute.resize(new_size);
        
        auto pos = thrust::make_zip_iterator(thrust::make_tuple(particles_to_redistribute.begin() + old_size, 
                                                                grids_to_redistribute.begin()     + old_size));
        thrust::copy(mid, end, pos);
        
        attribs.resize(new_num_particles);
    }

    const int num_grids = m_ba.size();
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
            const size_t old_size = m_particles[i].attribs.size();
            const size_t new_size = old_size + num_to_add;
            m_particles[i].attribs.resize(new_size);
            thrust::copy(begin, end, m_particles[i].attribs.begin() + old_size);
            m_particles[i].temp.resize(m_particles[i].attribs.size());
        }
        else
        {           
            char* dst;
            const size_t old_size = not_ours[dest_proc].size();
            const size_t new_size 
                = old_size + num_to_add*superparticle_size + sizeof(size_t) + 2*sizeof(int);
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
 
void
ElectromagneticParticleContainer::
RedistributeMPI(std::map<int, thrust::device_vector<char> >& not_ours)
{


    BL_PROFILE("ParticleContainer::RedistributeMPI()");
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

                StructOfArrays<PIdx::nattribs, 0> redistributed_particles;

                AMREX_ALWAYS_ASSERT(pid == ParallelDescriptor::MyProc());
                {
                    const size_t old_size = redistributed_particles.size();
                    const size_t new_size = old_size + num_particles;        
                    redistributed_particles.resize(new_size);
                    
                    for (int j = 0; j < PIdx::nattribs; ++j) {
                        auto& attrib = redistributed_particles.GetRealData(j);
                        cudaMemcpy(attrib.data() + old_size, buffer, num_particles*sizeof(Real), cudaMemcpyHostToDevice);
                        buffer += num_particles*sizeof(Real);
                    }
                }
            
                {
                    const size_t old_size = m_particles[gid].attribs.size();
                    const size_t new_size = old_size + num_particles;
                    m_particles[gid].attribs.resize(new_size);
                    thrust::copy(redistributed_particles.begin(),
                                 redistributed_particles.end(),
                                 m_particles[gid].attribs.begin() + old_size);
                    m_particles[gid].temp.resize(m_particles[gid].attribs.size());
                }
            }
        }
    }
        
#endif // MPI
 
}
