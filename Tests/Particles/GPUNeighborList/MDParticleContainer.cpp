#include "MDParticleContainer.H"
#include "Constants.H"

#include <thrust/reduce.h>

#include "md_K.H"

using namespace amrex;

namespace
{    
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
}

MDParticleContainer::
MDParticleContainer(const Geometry            & a_geom,
                    const DistributionMapping & a_dmap,
                    const BoxArray            & a_ba)
    : ParticleContainer<0, 0, PIdx::nattribs, 0>(a_geom, a_dmap, a_ba)
{}

void
MDParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std,
              const Real     a_thermal_momentum_mean)
{
    BL_PROFILE("MDParticleContainer::InitParticles");

    amrex::Print() << "Generating particles... \n";

    const int lev = 0;   
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    
    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                     *a_num_particles_per_cell[1],
                                     *a_num_particles_per_cell[2]);

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);    

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        Cuda::HostVector<ParticleType> host_particles;
        std::array<Cuda::HostVector<Real>, PIdx::nattribs> host_attribs;
        
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real v[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                get_gaussian_random_momentum(v, a_thermal_momentum_mean,
                                             a_thermal_momentum_std);
                
                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();                
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;
                
                attribs[PIdx::vx] = v[0];
                attribs[PIdx::vy] = v[1];
                attribs[PIdx::vz] = v[2];

                attribs[PIdx::ax] = 0.0;
                attribs[PIdx::ay] = 0.0;
                attribs[PIdx::az] = 0.0;
                
                host_particles.push_back(p);
                for (int kk = 0; kk < PIdx::nattribs; ++kk)
                    host_attribs[kk].push_back(attribs[kk]);                
            }
        }
        
        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);
        
        Cuda::thrust_copy(host_particles.begin(),
                          host_particles.end(),
                          particle_tile.GetArrayOfStructs().begin() + old_size);
        
        for (int kk = 0; kk < PIdx::nattribs; ++kk)
        {
            Cuda::thrust_copy(host_attribs[kk].begin(),
                              host_attribs[kk].end(),
                              particle_tile.GetStructOfArrays().GetRealData(kk).begin() + old_size);
        }
    }
}

void MDParticleContainer::BuildNeighborList()
{
    BL_PROFILE("MDParticleContainer::BuildNeighborList");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dxi = Geom(lev).InvCellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();
    const BoxArray& ba   = ParticleBoxArray(lev);
    const DistributionMapping& dmap = ParticleDistributionMap(lev);
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        
        const Box& bx = mfi.tilebox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        auto& ptile = plev[std::make_pair(gid, tid)];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();

        thrust::device_vector<unsigned int> cells(np);
        unsigned int* pcell = thrust::raw_pointer_cast(cells.data());

        thrust::device_vector<unsigned int> counts(bx.numPts());
        unsigned int* pcount = thrust::raw_pointer_cast(counts.data());

        thrust::device_vector<unsigned int> offsets(bx.numPts() + 1, np);
        unsigned int* poffset = thrust::raw_pointer_cast(offsets.data());

        thrust::device_vector<unsigned int> permutation(np);
        unsigned int* pperm = thrust::raw_pointer_cast(permutation.data());

        // First we build the cell list data structure
        
        ParticleType* pstruct = &(aos[0]);
        AMREX_FOR_1D ( np, i,
        {
            int ix = (pstruct[i].pos(0)-plo[0])*dxi[0] - lo.x;
            int iy = (pstruct[i].pos(1)-plo[1])*dxi[1] - lo.y;
            int iz = (pstruct[i].pos(2)-plo[2])*dxi[2] - lo.z;
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;            
            unsigned int uix = amrex::min(nx,amrex::max(0,ix));
            unsigned int uiy = amrex::min(ny,amrex::max(0,iy));
            unsigned int uiz = amrex::min(nz,amrex::max(0,iz));
            pcell[i] = (uix * ny + uiy) * nz + uiz; 
            Cuda::Atomic::Add(&pcount[pcell[i]], 1u);
        });

        thrust::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

        thrust::copy(offsets.begin(), offsets.end()-1, counts.begin());

        AMREX_FOR_1D ( np, i,
        {
            unsigned int index = atomicInc(&pcount[pcell[i]],
                                           std::numeric_limits<unsigned int>::max());
            pperm[index] = i;
        });

        // Now count the number of neighbors for each particle

        thrust::device_vector<unsigned int> nbor_counts(np);
        unsigned int* pnbor_counts = thrust::raw_pointer_cast(nbor_counts.data());
        
        AMREX_FOR_1D ( np, i,
        {
            int ix = (pstruct[i].pos(0)-plo[0])*dxi[0] - lo.x;
            int iy = (pstruct[i].pos(1)-plo[1])*dxi[1] - lo.y;
            int iz = (pstruct[i].pos(2)-plo[2])*dxi[2] - lo.z;

            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;            

            int count = 0;

            for (int ii = amrex::max(ix-1, 0); ii <= amrex::min(ix+1, nx-1); ++ii) {
                for (int jj = amrex::max(iy-1, 0); jj <= amrex::min(iy+1, ny-1); ++jj) {
                    for (int kk = amrex::max(iz-1, 0); kk <= amrex::min(iz+1, nz-1); ++kk) {
                        int index = (ii * ny + jj) * nz + kk;
                        for (int p = poffset[index]; p < poffset[index+1]; ++p) {
                            if (pperm[p] == i) continue;
                            if (check_pair(pstruct[i], pstruct[pperm[p]]))
                                count += 1;
                        }
                    }
                }
            }

            pnbor_counts[i] = count;
        });

        // Now we can allocate and build our neighbor list

        const size_t total_nbors = thrust::reduce(nbor_counts.begin(), nbor_counts.end());
        thrust::device_vector<unsigned int> nbor_offsets(np + 1, total_nbors);
        unsigned int* pnbor_offset = thrust::raw_pointer_cast(nbor_offsets.data());

        thrust::exclusive_scan(nbor_counts.begin(), nbor_counts.end(), nbor_offsets.begin());
                
        thrust::device_vector<unsigned int> nbor_list(total_nbors);
        unsigned int* pnbor_list = thrust::raw_pointer_cast(nbor_list.data());

        AMREX_FOR_1D ( np, i,
        {
            int ix = (pstruct[i].pos(0)-plo[0])*dxi[0] - lo.x;
            int iy = (pstruct[i].pos(1)-plo[1])*dxi[1] - lo.y;
            int iz = (pstruct[i].pos(2)-plo[2])*dxi[2] - lo.z;
            
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;            

            int n = 0;            
            for (int ii = amrex::max(ix-1, 0); ii <= amrex::min(ix+1, nx-1); ++ii) {
                for (int jj = amrex::max(iy-1, 0); jj <= amrex::min(iy+1, ny-1); ++jj) {
                    for (int kk = amrex::max(iz-1, 0); kk <= amrex::min(iz+1, nz-1); ++kk) {
                        int index = (ii * ny + jj) * nz + kk;
                        for (int p = poffset[index]; p < poffset[index+1]; ++p) {
                            if (pperm[p] == i) continue;
                            if (check_pair(pstruct[i], pstruct[pperm[p]])) {
                                pnbor_list[pnbor_offset[i] + n] = pperm[p]; 
                                ++n;
                            }
                        }
                    }
                }
            }
        });

        // for (int i = 0; i < np; ++i) {
        //     amrex::Print() << "Particle " << i << " will collide with: ";
        //     for (int j = nbor_offsets[i]; j < nbor_offsets[i+1]; ++j) {
        //         amrex::Print() << nbor_list[j] << " ";
        //     }
        //     amrex::Print() << "\n";
        // }
    }
}
