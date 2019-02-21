#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_

#include <AMReX_Particles.H>
#include <AMReX_CudaContainers.H>

#include <thrust/reduce.h>

#include "Constants.H"
#include "md_K.H"

namespace amrex
{

class NeighborList
{
public:

    template <class AoS>
    void build (const AoS& aos, const amrex::Box& bx, const amrex::Geometry& geom)
    {
        BL_PROFILE("NeighborList::build");
        
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
                
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const size_t np = aos.numTotalParticles();
        
        Gpu::DeviceVector<unsigned int> cells(np);
        unsigned int* pcell = cells.dataPtr();
        
        Gpu::DeviceVector<unsigned int> counts(bx.numPts(), 0);
        unsigned int* pcount = counts.dataPtr();
        
        Gpu::DeviceVector<unsigned int> offsets(bx.numPts() + 1, np);
        unsigned int* poffset = offsets.dataPtr();
    
        Gpu::DeviceVector<unsigned int> permutation(np);
        unsigned int* pperm = permutation.dataPtr();
        
        // First we build the cell list data structure
        
        const typename AoS::ParticleType* pstruct = &(aos[0]);
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
        
        constexpr unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();
        
        AMREX_FOR_1D ( np, i,
        {
            unsigned int index = Cuda::Atomic::Inc(&pcount[pcell[i]], max_unsigned_int);
            pperm[index] = i;
        });
    
        // Now count the number of neighbors for each particle
        
        Gpu::DeviceVector<unsigned int> nbor_counts(np);
        unsigned int* pnbor_counts = nbor_counts.dataPtr();
        
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
        m_nbor_offsets.resize(np + 1, total_nbors);
        unsigned int* pnbor_offset = m_nbor_offsets.dataPtr();
        
        thrust::exclusive_scan(nbor_counts.begin(), nbor_counts.end(),
                               m_nbor_offsets.begin());
        
        m_nbor_list.resize(total_nbors);
        unsigned int* pm_nbor_list = m_nbor_list.dataPtr();
        
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
                                pm_nbor_list[pnbor_offset[i] + n] = pperm[p]; 
                                ++n;
                            }
                        }
                    }
                }
            }
        });
    }

          unsigned int * offsetPtr ()       { return m_nbor_offsets.dataPtr(); }
    const unsigned int * offsetPtr () const { return m_nbor_offsets.dataPtr(); }

          unsigned int * nborPtr ()       { return m_nbor_list.dataPtr(); }
    const unsigned int * nborPtr () const { return m_nbor_list.dataPtr(); }

    int numParticles () { return m_nbor_offsets.size() - 1; }

    void print ();

protected:

    // This is the neighbor list data structure
    amrex::Gpu::DeviceVector<unsigned int> m_nbor_offsets;
    amrex::Gpu::DeviceVector<unsigned int> m_nbor_list;
};

}

#endif
