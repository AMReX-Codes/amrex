#ifndef AMREX_PARTICLE_LOCATOR_H_
#define AMREX_PARTICLE_LOCATOR_H_

#include <AMReX_Particles.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_Tuple.H>

namespace amrex
{

struct AssignGrid
{
    const Box* m_boxes_ptr;
    const unsigned int* m_poffset;
    const unsigned int* m_pperm;

    Dim3 m_lo;
    Dim3 m_hi;
    Dim3 m_bin_size;
    Dim3 m_num_bins;

	Box m_domain;
    GpuArray<Real, AMREX_SPACEDIM> m_plo;
    GpuArray<Real, AMREX_SPACEDIM> m_dxi;
	
    AMREX_GPU_HOST_DEVICE
    AssignGrid () {}

    AssignGrid (const Box* a_boxes_ptr, const unsigned int* a_poffset, const unsigned int* a_pperm,
                const IntVect& a_bins_lo, const IntVect& a_bins_hi, const IntVect& a_bin_size,
                const IntVect& a_num_bins, const Geometry& a_geom)
        : m_boxes_ptr(a_boxes_ptr),
          m_poffset(a_poffset), m_pperm(a_pperm),
          m_lo(a_bins_lo.dim3()), m_hi(a_bins_hi.dim3()), m_bin_size(a_bin_size.dim3()),
          m_num_bins(a_num_bins.dim3()), m_domain(a_geom.Domain()),
          m_plo(a_geom.ProbLoArray()), m_dxi(a_geom.InvCellSizeArray())
        {
            // clamp bin size and num_bins to 1 for AMREX_SPACEDIM < 3
            m_bin_size.x = amrex::max(m_bin_size.x, 1);
            m_bin_size.y = amrex::max(m_bin_size.y, 1);
            m_bin_size.z = amrex::max(m_bin_size.z, 1);

            m_num_bins.x = amrex::max(m_num_bins.x, 1);
            m_num_bins.y = amrex::max(m_num_bins.y, 1);
            m_num_bins.z = amrex::max(m_num_bins.z, 1);
        }

    template <typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const P& p) const noexcept
    {
        const auto iv = getParticleCell(p, m_plo, m_dxi, m_domain);
        return this->operator()(iv);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const IntVect& iv) const noexcept
    {
        const auto lo = iv.dim3();
        int ix = (lo.x - m_lo.x) / m_bin_size.x;
        int iy = (lo.y - m_lo.y) / m_bin_size.y;
        int iz = (lo.z - m_lo.z) / m_bin_size.z;

        int nx = m_num_bins.x;
        int ny = m_num_bins.y;
        int nz = m_num_bins.z;

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

    void build (const BoxArray& ba, const Geometry& geom)
    {
        m_defined = true;
        m_ba = ba;
        m_geom = geom;
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

        m_num_bins = (m_bins_hi - m_bins_lo + m_bin_size) / m_bin_size; 
               
        int num_bins_total = AMREX_D_TERM(m_num_bins[0],*m_num_bins[1],*m_num_bins[2]);

        m_cells.resize(num_boxes);

        m_counts.resize(0);
        m_counts.resize(num_bins_total, 0);

        m_offsets.resize(0);
        m_offsets.resize(num_bins_total+1, num_boxes);

        m_permutation.resize(num_boxes);

        const auto lo       = m_bins_lo.dim3();
        const auto hi       = m_bins_hi.dim3();
        const auto bin_size = m_bin_size.dim3();

        auto nbins    = m_num_bins.dim3();
        nbins.x = amrex::max(nbins.x, 1);
        nbins.y = amrex::max(nbins.y, 1);
        nbins.z = amrex::max(nbins.z, 1);
        
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
            int nx = nbins.x;
            int ny = nbins.y;
            int nz = nbins.z;
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

    void setGeometry (const Geometry& a_geom) noexcept
    {
        AMREX_ASSERT(m_defined);
        m_geom = a_geom;
    }
    
    AssignGrid getGridAssignor () const noexcept 
    {
        AMREX_ASSERT(m_defined);
        return AssignGrid(m_device_boxes.dataPtr(), m_offsets.dataPtr(), m_permutation.dataPtr(),
                          m_bins_lo, m_bins_hi, m_bin_size, m_num_bins, m_geom);
    }

    bool isValid (const BoxArray& ba) const noexcept
    {
        if (m_defined) return BoxArray::SameRefs(m_ba, ba);
        return false;
    }
            
protected:

    bool m_defined;
    
    BoxArray m_ba;
    Geometry m_geom;

    IntVect m_bins_lo;
    IntVect m_bins_hi;
    IntVect m_bin_size;
    IntVect m_num_bins;

    Gpu::HostVector<Box> m_host_boxes;
    Gpu::DeviceVector<Box> m_device_boxes;

    Gpu::ManagedDeviceVector<unsigned int> m_counts;
    Gpu::ManagedDeviceVector<unsigned int> m_cells;
    Gpu::ManagedDeviceVector<unsigned int> m_offsets;
    Gpu::ManagedDeviceVector<unsigned int> m_permutation;
};

struct AmrAssignGrid
{
    const AssignGrid* m_funcs;
    std::size_t m_size;
    
    AmrAssignGrid(const AssignGrid* a_funcs, std::size_t a_size)
        : m_funcs(a_funcs), m_size(a_size)
        {}

    template <typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE   
    GpuTuple<int, int> operator() (const P& p) const noexcept
    {
        for (int lev = m_size-1; lev >= 0; --lev)
        {
            int grid = m_funcs[lev](p);
            if (grid >= 0) return makeTuple(grid, lev);
        }
        return makeTuple(-1, -1);
    }
};

class AmrParticleLocator
{
    Vector<ParticleLocator> m_locators;
    Gpu::ManagedDeviceVector<AssignGrid> m_grid_assignors;
    bool m_defined = false;

public:

    AmrParticleLocator() {}
    
    AmrParticleLocator(const Vector<BoxArray>& a_ba,
                       const Vector<Geometry>& a_geom)
    {
        build(a_ba, a_geom);
    }

    AmrParticleLocator(const ParGDBBase* a_gdb)
    {
        build(a_gdb);
    }
        
    void build (const Vector<BoxArray>& a_ba,
                const Vector<Geometry>& a_geom)
    {
        m_defined = true;
        int num_levels = a_ba.size();
        m_locators.resize(num_levels);
        m_grid_assignors.resize(num_levels);
        for (int lev = 0; lev < num_levels; ++lev)
        {
            m_locators[lev].build(a_ba[lev], a_geom[lev]);
            m_grid_assignors[lev] = m_locators[lev].getGridAssignor();
        }        
    }
    
    void build (const ParGDBBase* a_gdb)
    {
        Vector<BoxArray> ba;
        Vector<Geometry> geom;
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
        {
            ba.push_back(a_gdb->ParticleBoxArray(lev));
            geom.push_back(a_gdb->Geom(lev));
        }
        build(ba, geom);
    }

    bool isValid (const Vector<BoxArray>& a_ba) const
    {
        if ( (m_locators.size() == 0) or !m_defined) return false;
        bool all_valid = true;
        int num_levels = m_locators.size();
        for (int lev = 0; lev < num_levels; ++lev)
            all_valid = all_valid && m_locators[lev].isValid(a_ba[lev]);
        return all_valid;
    }
    
    bool isValid (const ParGDBBase* a_gdb) const
    {
        Vector<BoxArray> ba;
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
            ba.push_back(a_gdb->ParticleBoxArray(lev));
        return this->isValid(ba);
    }

    void setGeometry (const ParGDBBase* a_gdb)
    {
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
        {
            m_locators[lev].setGeometry(a_gdb->Geom(lev));
            m_grid_assignors[lev] = m_locators[lev].getGridAssignor();
        }        
    }
    
    AmrAssignGrid getGridAssignor () const noexcept 
    {
        AMREX_ASSERT(m_defined);
        return AmrAssignGrid(m_grid_assignors.dataPtr(), m_locators.size());
    }
};

}

#endif
