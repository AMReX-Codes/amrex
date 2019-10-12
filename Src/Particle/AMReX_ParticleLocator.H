#ifndef AMREX_PARTICLE_LOCATOR_H_
#define AMREX_PARTICLE_LOCATOR_H_

#include <AMReX_Particles.H>
#include <AMReX_CudaContainers.H>

namespace amrex
{

struct assignGrid
{
    const Box* m_boxes_ptr;
    const unsigned int* m_poffset;
    const unsigned int* m_pperm;

    Dim3 m_lo;
    Dim3 m_hi;
    Dim3 m_bin_size;
    
    assignGrid (const Box* a_boxes_ptr, const unsigned int* a_poffset, const unsigned int* a_pperm,
                const IntVect& a_bins_lo, const IntVect& a_bins_hi, const IntVect& a_bin_size)
        : m_boxes_ptr(a_boxes_ptr),
          m_poffset(a_poffset), m_pperm(a_pperm),
          m_lo(a_bins_lo.dim3()), m_hi(a_bins_hi.dim3()), m_bin_size(a_bin_size.dim3())
        {
            // clamp bin size to 1 for AMREX_SPACEDIM < 3
            m_bin_size.x = amrex::max(m_bin_size.x, 1);
            m_bin_size.y = amrex::max(m_bin_size.y, 1);
            m_bin_size.z = amrex::max(m_bin_size.z, 1);
        }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const IntVect& iv) const noexcept
    {
        const auto lo = iv.dim3();
        int ix = (lo.x - m_lo.x) / m_bin_size.x;
        int iy = (lo.y - m_lo.y) / m_bin_size.y;
        int iz = (lo.z - m_lo.z) / m_bin_size.z;

        int nx = m_hi.x-m_lo.x+1;
        int ny = m_hi.y-m_lo.y+1;
        int nz = m_hi.z-m_lo.z+1;

        for (int ii = amrex::max(ix-1, 0); ii <= amrex::min(ix, nx-1); ++ii) {
            for (int jj = amrex::max(iy-1, 0); jj <= amrex::min(iy, ny-1); ++jj) {
                for (int kk = amrex::max(iz-1, 0); kk <= amrex::min(iz, nz-1); ++kk) {
                    int index = (ii * ny + jj) * nz + kk;
                    for (int p = m_poffset[index]; p < m_poffset[index+1]; ++p) {
                        const Box& box = m_boxes_ptr[m_pperm[p]];
                        if (box.contains(iv)) return m_pperm[p];
                    }
                }
            }
        }

        return -1;        
    }
};

class ParticleLocator
{
public:

    ParticleLocator () : m_defined(false) {}

    void build (const BoxArray& ba)
    {
        m_defined = true;
        m_ba = ba;
        int num_boxes = ba.size();
        m_host_boxes.resize(0);
        for (int i = 0; i < num_boxes; ++i) m_host_boxes.push_back(ba[i]);

        m_device_boxes.resize(num_boxes);
        Gpu::thrust_copy(m_host_boxes.begin(), m_host_boxes.end(), m_device_boxes.begin());

        // compute the lo, hi and the max box size in each direction
        ReduceOps<AMREX_D_DECL(ReduceOpMin, ReduceOpMin, ReduceOpMin),
                  AMREX_D_DECL(ReduceOpMax, ReduceOpMax, ReduceOpMax),
                  AMREX_D_DECL(ReduceOpMax, ReduceOpMax, ReduceOpMax)> reduce_op;
        ReduceData<AMREX_D_DECL(int, int, int), 
                   AMREX_D_DECL(int, int, int), 
                   AMREX_D_DECL(int, int, int)> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const auto boxes_ptr = m_device_boxes.dataPtr();
        reduce_op.eval(num_boxes, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            const Box& box = boxes_ptr[i];
            IntVect lo = box.smallEnd();
            IntVect hi = box.bigEnd();
            IntVect si = box.length();
            return {AMREX_D_DECL(lo[0], lo[1], lo[2]), 
                    AMREX_D_DECL(hi[0], hi[1], hi[2]), 
                    AMREX_D_DECL(si[0], si[1], si[2])};
        }); 
        
        ReduceTuple hv = reduce_data.value();

        m_bins_lo  = IntVect(AMREX_D_DECL(amrex::get<0>(hv), 
                                          amrex::get<1>(hv), 
                                          amrex::get<2>(hv)));
        m_bins_hi  = IntVect(AMREX_D_DECL(amrex::get< AMREX_SPACEDIM  >(hv), 
                                          amrex::get< AMREX_SPACEDIM+1>(hv), 
                                          amrex::get< AMREX_SPACEDIM+2>(hv)));
        m_bin_size = IntVect(AMREX_D_DECL(amrex::get<2*AMREX_SPACEDIM>(hv), 
                                          amrex::get<2*AMREX_SPACEDIM+1>(hv), 
                                          amrex::get<2*AMREX_SPACEDIM+2>(hv)));

        m_bins_hi = (m_bins_hi + m_bin_size - IntVect::TheUnitVector()) / m_bin_size; 

        int num_bins = AMREX_D_TERM((m_bins_hi[0] - m_bins_lo[0] + 1),
                                   *(m_bins_hi[1] - m_bins_lo[1] + 1),
                                   *(m_bins_hi[2] - m_bins_lo[2] + 1));

        m_cells.resize(num_boxes);

        m_counts.resize(0);
        m_counts.resize(num_bins, 0);

        m_offsets.resize(0);
        m_offsets.resize(num_bins+1, num_boxes);

        m_permutation.resize(num_boxes);

        const auto lo       = m_bins_lo.dim3();
        const auto hi       = m_bins_hi.dim3();
        const auto bin_size = m_bin_size.dim3();

        unsigned int* pcell = m_cells.dataPtr();
        unsigned int* pcount = m_counts.dataPtr();
        unsigned int* pperm = m_permutation.dataPtr();

        amrex::ParallelFor(num_boxes, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            const Box& box = boxes_ptr[i];
            const auto blo = amrex::lbound(box);
            int ix = (blo.x - lo.x) / bin_size.x;
            int iy = (blo.y - lo.y) / bin_size.y;
            int iz = (blo.z - lo.z) / bin_size.z;
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;            
            unsigned int uix = amrex::min(nx-1,amrex::max(0,ix));
            unsigned int uiy = amrex::min(ny-1,amrex::max(0,iy));
            unsigned int uiz = amrex::min(nz-1,amrex::max(0,iz));
            pcell[i] = (uix * ny + uiy) * nz + uiz;
            Gpu::Atomic::Add(&pcount[pcell[i]], 1u);
        });

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());
        
        Gpu::thrust_copy(m_offsets.begin(), m_offsets.end()-1, m_counts.begin());

        constexpr unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

        AMREX_FOR_1D(num_boxes, i,
        {
            unsigned int index = Gpu::Atomic::Inc(&pcount[pcell[i]], max_unsigned_int);
            pperm[index] = i;
        });
    }

    assignGrid getGridAssignor () const noexcept 
    {
        AMREX_ASSERT(m_defined);
        return assignGrid(m_device_boxes.dataPtr(), m_offsets.dataPtr(), m_permutation.dataPtr(),
                          m_bins_lo, m_bins_hi, m_bin_size);
    }

    bool isValid (const BoxArray& ba) const noexcept
    {
        if (m_defined) return BoxArray::SameRefs(m_ba, ba);
        return false;
    }
        
    
protected:

    bool m_defined;
    
    BoxArray m_ba;
    
    IntVect m_bins_lo;
    IntVect m_bins_hi;
    IntVect m_bin_size;

    Gpu::HostVector<Box> m_host_boxes;
    Gpu::DeviceVector<Box> m_device_boxes;

    Gpu::ManagedDeviceVector<unsigned int> m_counts;
    Gpu::ManagedDeviceVector<unsigned int> m_cells;
    Gpu::ManagedDeviceVector<unsigned int> m_offsets;
    Gpu::ManagedDeviceVector<unsigned int> m_permutation;
};

}

#endif
