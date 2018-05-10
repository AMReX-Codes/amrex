#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

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
        Real right_edge[AMREX_SPACEDIM];
        int  num_cells[AMREX_SPACEDIM];
        
        Real domain_width[AMREX_SPACEDIM];
        Real inverse_dx[AMREX_SPACEDIM];
        Real dx[AMREX_SPACEDIM];
        
        DeviceDomain(const Geometry& geom) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                left_edge[i] = geom.ProbLo(i);            
                right_edge[i] = geom.ProbHi(i);            
                num_cells[i] = geom.Domain().length(i);            
                domain_width[i] = right_edge[i] - left_edge[i];            
                dx[i] = domain_width[i] / num_cells[i];
                inverse_dx[i] = 1.0/dx[i];
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
    m_mask_ptr.reset(new MultiFab(m_ba, m_dmap, 1, ng));
    m_mask_ptr->setVal(-1);
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const int grid_id = mfi.index();
        m_mask_ptr->setVal(grid_id, box, 0, 1);
    }
    m_mask_ptr->FillBoundary(m_geom.periodicity());
}

void
ElectromagneticParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std, 
              const Real     a_thermal_momentum_mean,
              const Real     a_density)
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
                get_gaussian_random_momentum(u, a_thermal_momentum_mean, a_thermal_momentum_std);
                
                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                RealBox rb(tile_box, dx, plo);

                Real point[3];
                point[0] = x;
                point[1] = y;
                point[2] = z;
                AMREX_ALWAYS_ASSERT(rb.contains(point));

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
    
    const Real* dx  = m_geom.CellSize();
    const Real* plo = m_geom.ProbLo();
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.size();

        FORT_LAUNCH_PARTICLES(np,
                              gather_magnetic_field,
                              np, 
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              particles.bx().data(), particles.by().data(), particles.bz().data(),
                              BL_TO_FORTRAN_3D(Bx[mfi]),
                              BL_TO_FORTRAN_3D(By[mfi]),
                              BL_TO_FORTRAN_3D(Bz[mfi]),
                              plo, dx);
        
        FORT_LAUNCH_PARTICLES(np,
                              gather_electric_field,
                              np, 
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              particles.ex().data(), particles.ey().data(), particles.ez().data(),
                              BL_TO_FORTRAN_3D(Ex[mfi]),
                              BL_TO_FORTRAN_3D(Ey[mfi]),
                              BL_TO_FORTRAN_3D(Ez[mfi]),
                              plo, dx);
        
#ifdef AMREX_USE_CUDA           
        cudaDeviceSynchronize();
#endif
        FORT_LAUNCH_PARTICLES(np, 
                              push_momentum_boris,
                              np,
                              particles.ux().data(), particles.uy().data(), particles.uz().data(),
                              particles.ginv().data(),
                              particles.ex().data(), particles.ey().data(), particles.ez().data(),
                              particles.bx().data(), particles.by().data(), particles.bz().data(),
                              m_charge, m_mass, dt);
        
#ifdef AMREX_USE_CUDA                        
        cudaDeviceSynchronize();
#endif      

        FORT_LAUNCH_PARTICLES(np,
                              push_position_boris,
                              np,
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              particles.ux().data(), particles.uy().data(), particles.uz().data(),
                              particles.ginv().data(), dt);
        
#ifdef AMREX_USE_CUDA            
        cudaDeviceSynchronize();
#endif

        FORT_LAUNCH_PARTICLES(np, deposit_current,
                              BL_TO_FORTRAN_3D(jx[mfi]),
                              BL_TO_FORTRAN_3D(jy[mfi]),
                              BL_TO_FORTRAN_3D(jz[mfi]),
                              np, 
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              particles.ux().data(), particles.uy().data(), particles.uz().data(),
                              particles.ginv().data(), particles.w().data(),
                              m_charge, plo, dt, dx);            
    }
}

void 
ElectromagneticParticleContainer::
PushParticlesOnly(amrex::Real dt)
{
    BL_PROFILE("ElectromagneticParticleContainer::PushParticlesOnly");
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.size();
        
        FORT_LAUNCH_PARTICLES(np, set_gamma,
                              np, 
                              particles.ux().data(), particles.uy().data(), particles.uz().data(),
                              particles.ginv().data());
        
#ifdef AMREX_USE_CUDA
        cudaDeviceSynchronize();
#endif        
        
        FORT_LAUNCH_PARTICLES(np, push_position_boris,
                              np,
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              particles.ux().data(), particles.uy().data(), particles.uz().data(),
                              particles.ginv().data(), dt);
    }
}

void
ElectromagneticParticleContainer::
EnforcePeriodicBCs()
{
    BL_PROFILE("ElectromagneticParticleContainer::EnforcePeriodicBCs");

    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.size();

        const Real* plo = m_geom.ProbLo();
        const Real* phi = m_geom.ProbHi();
        
        FORT_LAUNCH_PARTICLES(np, enforce_periodic,
                              np,
                              particles.x().data(),  particles.y().data(),  particles.z().data(),
                              plo, phi);
    }
}

void
ElectromagneticParticleContainer::
OK()
{
    BL_PROFILE("ElectromagneticParticleContainer::OK");

    thrust::device_vector<int> grid_indices;

    for(MFIter mfi(*m_mask_ptr); mfi.isValid(); ++mfi) {
        int i = mfi.index();
        
        auto& particles = m_particles[i];
        const int np = particles.size();
        auto& soa = particles.attribs;
        
        grid_indices.resize(np);
        
        thrust::transform(soa.begin(),
                          soa.end(),
                          grid_indices.begin(),
                          assignParticle(DeviceDomain(m_geom),
                                         DeviceBox((*m_mask_ptr)[i].box()),
                                         (*m_mask_ptr)[i].dataPtr()));

        int count = thrust::count_if(grid_indices.begin(),
                                     grid_indices.end(),
                                     grid_is_not(i));

        AMREX_ALWAYS_ASSERT(count == 0);
    }
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
                          assignParticle(DeviceDomain(m_geom),
                                         DeviceBox((*m_mask_ptr)[src_grid].box()),
                                         (*m_mask_ptr)[src_grid].dataPtr()));
                
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

    for (int i = 0; i < num_grids; ++i)
    {
        auto begin = particles_to_redistribute.begin();
        thrust::advance(begin, start[i]);

        auto end = particles_to_redistribute.begin();
        thrust::advance(end, stop[i]);
     
        const size_t num_to_add = stop[i] - start[i];
        const size_t old_size = m_particles[i].attribs.size();
        const size_t new_size = old_size + num_to_add;

        m_particles[i].attribs.resize(new_size);

        thrust::copy(begin, end, m_particles[i].attribs.begin() + old_size);

        m_particles[i].temp.resize(m_particles[i].attribs.size());
    }
}
