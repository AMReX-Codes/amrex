#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_

#include <AMReX_Particles.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_DenseBins.H>

namespace amrex
{

template <class ParticleType>
struct Neighbors
{
    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (int start, int stop, const unsigned int * nbor_list_ptr, ParticleType* pstruct)
            : m_index(start), m_stop(stop), m_nbor_list_ptr(nbor_list_ptr), m_pstruct(pstruct)
        {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () { ++m_index;; }

        AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& rhs) const { return m_index < m_stop; }

        AMREX_GPU_HOST_DEVICE
        ParticleType& operator* () const { return m_pstruct[m_nbor_list_ptr[m_index]];  }

    private:
        int m_index;
        int m_stop;
        const unsigned int* m_nbor_list_ptr;
        ParticleType* m_pstruct;
    };

    struct const_iterator
    {
        AMREX_GPU_HOST_DEVICE
        const_iterator (int start, int stop, const unsigned int * nbor_list_ptr, const ParticleType* pstruct)
            : m_index(start), m_stop(stop), m_nbor_list_ptr(nbor_list_ptr), m_pstruct(pstruct)
        {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () { ++m_index;; }

        AMREX_GPU_HOST_DEVICE
        bool operator!= (const_iterator const& rhs) const { return m_index < m_stop; }

        AMREX_GPU_HOST_DEVICE
        const ParticleType& operator* () const { return m_pstruct[m_nbor_list_ptr[m_index]];  }

    private:
        int m_index;
        int m_stop;
        const unsigned int* m_nbor_list_ptr;
        const ParticleType* m_pstruct;
    };

    AMREX_GPU_HOST_DEVICE
    iterator begin () noexcept {
        return iterator(m_nbor_offsets_ptr[m_i], m_nbor_offsets_ptr[m_i+1],
                        m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    iterator end () noexcept {
        return iterator(m_nbor_offsets_ptr[m_i+1], m_nbor_offsets_ptr[m_i+1],
                        m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    const_iterator begin () const noexcept {
        return const_iterator(m_nbor_offsets_ptr[m_i], m_nbor_offsets_ptr[m_i+1],
                              m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    const_iterator end () const noexcept {
        return const_iterator(m_nbor_offsets_ptr[m_i+1], m_nbor_offsets_ptr[m_i+1],
                              m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    const_iterator cbegin () const noexcept {
        return const_iterator(m_nbor_offsets_ptr[m_i], m_nbor_offsets_ptr[m_i+1],
                              m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    const_iterator cend () const noexcept {
        return const_iterator(m_nbor_offsets_ptr[m_i+1], m_nbor_offsets_ptr[m_i+1],
                              m_nbor_list_ptr, m_pstruct);
    }

    AMREX_GPU_HOST_DEVICE
    Neighbors (int i, const unsigned int *nbor_offsets_ptr, const unsigned int *nbor_list_ptr,
               ParticleType* pstruct)
        : m_i(i),
          m_nbor_offsets_ptr(nbor_offsets_ptr),
          m_nbor_list_ptr(nbor_list_ptr),
          m_pstruct(pstruct)
    {}

private:

    const int m_i;
    const unsigned int * m_nbor_offsets_ptr;
    const unsigned int * m_nbor_list_ptr;
    ParticleType * m_pstruct;
};

template <class ParticleType>
struct NeighborData
{
    NeighborData (const Gpu::DeviceVector<unsigned int>& offsets,
                  const Gpu::DeviceVector<unsigned int>& list,
                  ParticleType* pstruct)
        : m_nbor_offsets_ptr(offsets.dataPtr()),
          m_nbor_list_ptr(list.dataPtr()),
          m_pstruct(pstruct)
    {}

    AMREX_GPU_HOST_DEVICE
    amrex::Neighbors<ParticleType> getNeighbors (int i) const
    { return Neighbors<ParticleType>(i, m_nbor_offsets_ptr, m_nbor_list_ptr, m_pstruct); }

    const unsigned int * m_nbor_offsets_ptr;
    const unsigned int * m_nbor_list_ptr;
    ParticleType * m_pstruct;
};

template <class ParticleType>
class NeighborList
{
public:

    template <class PTile, class CheckPair>
    void build (PTile& ptile,
                const amrex::Box& bx, const amrex::Geometry& geom,
                CheckPair&& check_pair, int num_cells=1)
    {
        auto& vec = ptile.GetArrayOfStructs()();
        m_pstruct = vec.dataPtr();

        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        const size_t np_real  = ptile.numRealParticles();
        const size_t np_total = vec.size();
        const ParticleType* pstruct_ptr = vec.dataPtr();

        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        m_bins.build(np_total, pstruct_ptr, bx,
                     [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> IntVect
                     {
                         AMREX_D_TERM(AMREX_ASSERT((p.pos(0)-plo[0])*dxi[0] - lo.x >= 0.0);,
                                      AMREX_ASSERT((p.pos(1)-plo[1])*dxi[1] - lo.y >= 0.0);,
                                      AMREX_ASSERT((p.pos(2)-plo[2])*dxi[2] - lo.z >= 0.0));

                         return IntVect(AMREX_D_DECL(static_cast<int>((p.pos(0)-plo[0])*dxi[0] - lo.x),
                                                     static_cast<int>((p.pos(1)-plo[1])*dxi[1] - lo.y),
                                                     static_cast<int>((p.pos(2)-plo[2])*dxi[2] - lo.z)));
                     });

        // first pass - count the number of neighbors for each particle
        m_nbor_counts.resize(np_real+1);
        m_nbor_offsets.resize(np_real+1);

        auto pnbor_counts = m_nbor_counts.dataPtr();
        auto pperm = m_bins.permutationPtr();
        auto poffset = m_bins.offsetsPtr();
        auto pnbor_offset = m_nbor_offsets.dataPtr();

        AMREX_FOR_1D ( np_real, i,
        {
            int ix = (pstruct_ptr[i].pos(0)-plo[0])*dxi[0] - lo.x;
            int iy = (pstruct_ptr[i].pos(1)-plo[1])*dxi[1] - lo.y;
            int iz = (pstruct_ptr[i].pos(2)-plo[2])*dxi[2] - lo.z;

            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;

            int count = 0;

            for (int ii = amrex::max(ix-num_cells, 0); ii <= amrex::min(ix+num_cells, nx-1); ++ii) {
                for (int jj = amrex::max(iy-num_cells, 0); jj <= amrex::min(iy+num_cells, ny-1); ++jj) {
                    for (int kk = amrex::max(iz-num_cells, 0); kk <= amrex::min(iz+num_cells, nz-1); ++kk) {
                        int index = (ii * ny + jj) * nz + kk;
                        for (auto p = poffset[index]; p < poffset[index+1]; ++p) {
                            if (pperm[p] == i) continue;
                            if (check_pair(pstruct_ptr[i], pstruct_ptr[pperm[p]]))
                                count += 1;
                        }
                    }
                }
            }

            pnbor_counts[i] = count;
        });

        Gpu::exclusive_scan(m_nbor_counts.begin(), m_nbor_counts.end(), m_nbor_offsets.begin());

        // Now we can allocate and build our neighbor list
        unsigned int total_nbors;
#ifdef AMREX_USE_GPU
        Gpu::dtoh_memcpy(&total_nbors,m_nbor_offsets.dataPtr()+np_real,sizeof(unsigned int));
#else
        std::memcpy(&total_nbors,m_nbor_offsets.dataPtr()+np_real,sizeof(unsigned int));
#endif

        m_nbor_list.resize(total_nbors);

        auto pm_nbor_list = m_nbor_list.dataPtr();

        AMREX_FOR_1D ( np_real, i,
        {
            int ix = (pstruct_ptr[i].pos(0)-plo[0])*dxi[0] - lo.x;
            int iy = (pstruct_ptr[i].pos(1)-plo[1])*dxi[1] - lo.y;
            int iz = (pstruct_ptr[i].pos(2)-plo[2])*dxi[2] - lo.z;

            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;

            int n = 0;
            for (int ii = amrex::max(ix-num_cells, 0); ii <= amrex::min(ix+num_cells, nx-1); ++ii) {
                for (int jj = amrex::max(iy-num_cells, 0); jj <= amrex::min(iy+num_cells, ny-1); ++jj) {
                    for (int kk = amrex::max(iz-num_cells, 0); kk <= amrex::min(iz+num_cells, nz-1); ++kk) {
                        int index = (ii * ny + jj) * nz + kk;
                        for (auto p = poffset[index]; p < poffset[index+1]; ++p) {
                            if (pperm[p] == i) continue;
                            if (check_pair(pstruct_ptr[i], pstruct_ptr[pperm[p]])) {
                                pm_nbor_list[pnbor_offset[i] + n] = pperm[p];
                                ++n;
                            }
                        }
                    }
                }
            }
        });
    }

    NeighborData<ParticleType> data ()
    {
        return NeighborData<ParticleType>(m_nbor_offsets, m_nbor_list, m_pstruct);
    }

    int numParticles () { return m_nbor_offsets.size() - 1; }

    Gpu::DeviceVector<unsigned int>&       GetOffsets ()       { return m_nbor_offsets; }
    const Gpu::DeviceVector<unsigned int>& GetOffsets () const { return m_nbor_offsets; }

    Gpu::DeviceVector<unsigned int>&       GetCounts ()       { return m_nbor_counts; }
    const Gpu::DeviceVector<unsigned int>& GetCounts () const { return m_nbor_counts; }

    Gpu::DeviceVector<unsigned int>&       GetList ()       { return m_nbor_list; }
    const Gpu::DeviceVector<unsigned int>& GetList () const { return m_nbor_list; }

    void print ()
    {
        BL_PROFILE("NeighborList::print");

        Gpu::HostVector<unsigned int> host_nbor_offsets(m_nbor_offsets.size());
        Gpu::HostVector<unsigned int> host_nbor_list(m_nbor_list.size());

        Gpu::copy(Gpu::deviceToHost, m_nbor_offsets.begin(), m_nbor_offsets.end(), host_nbor_offsets.begin());
        Gpu::copy(Gpu::deviceToHost, m_nbor_list.begin(), m_nbor_list.end(), host_nbor_list.begin());

        for (int i = 0; i < numParticles(); ++i) {
            amrex::Print() << "Particle " << i << " could collide with: ";
            for (int j = host_nbor_offsets[i]; j < host_nbor_offsets[i+1]; ++j) {
                amrex::Print() << host_nbor_list[j] << " ";
            }
            amrex::Print() << "\n";
        }
    }

protected:

    ParticleType* m_pstruct;

    // This is the neighbor list data structure
    Gpu::DeviceVector<unsigned int> m_nbor_offsets;
    Gpu::DeviceVector<unsigned int> m_nbor_list;
    Gpu::DeviceVector<unsigned int> m_nbor_counts;

    DenseBins<ParticleType> m_bins;
};

}

#endif
