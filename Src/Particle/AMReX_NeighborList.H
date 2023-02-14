#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_
#include <AMReX_Config.H>

#include <AMReX_Particles.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_DenseBins.H>

namespace amrex
{

namespace
{
    // SelectActualNeighbors
    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& dst_tile,
                          N1 i, N2 j)
        noexcept -> decltype(check_pair(src_tile.m_aos[i], dst_tile.m_aos[j]))
    {
        return check_pair(src_tile.m_aos[i], dst_tile.m_aos[j]);
    }

    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& /*dst_tile*/,
                          N1 i, N2 j)
        noexcept -> decltype(check_pair(src_tile.m_aos, i, j))
    {
        return check_pair(src_tile.m_aos, i, j);
    }

    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& /*dst_tile*/,
                          N1 i, N2 j)
        noexcept -> decltype(check_pair(src_tile, i, j))
    {
        return check_pair(src_tile, i, j);
    }

    // NeighborList Build
    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2, typename N3, typename N4, typename N5>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& dst_tile,
                          N1 i, N2 j, N3 /*type*/, N4 /*ghost_i*/, N5 /*ghost_pid*/)
        noexcept -> decltype(check_pair(src_tile.m_aos[i], dst_tile.m_aos[j]))
    {
        return check_pair(src_tile.m_aos[i], dst_tile.m_aos[j]);
    }

    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2, typename N3, typename N4, typename N5>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& /*dst_tile*/,
                          N1 i, N2 j, N3 /*type*/, N4 /*ghost_i*/, N5 /*ghost_pid*/)
        noexcept -> decltype(check_pair(src_tile.m_aos, i, j))
    {
        return check_pair(src_tile.m_aos, i, j);
    }

    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2, typename N3, typename N4, typename N5>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& /*dst_tile*/,
                          N1 i, N2 j, N3 /*type*/, N4 /*ghost_i*/, N5 /*ghost_pid*/)
        noexcept -> decltype(check_pair(src_tile, i, j))
    {
        return check_pair(src_tile, i, j);
    }

    template <typename F,
              typename SrcData, typename DstData,
              typename N1, typename N2, typename N3, typename N4, typename N5>
    AMREX_GPU_HOST_DEVICE
    auto call_check_pair (F const& check_pair,
                          const SrcData& src_tile, const DstData& /*dst_tile*/,
                          N1 i, N2 j, N3 type, N4 ghost_i, N5 ghost_pid)
        noexcept -> decltype(check_pair(src_tile, i, j, type, ghost_i, ghost_pid))
    {
        return check_pair(src_tile, i, j, type, ghost_i, ghost_pid);
    }
}

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
        bool operator!= (iterator const& /*rhs*/) const { return m_index < m_stop; }

        AMREX_GPU_HOST_DEVICE
        ParticleType& operator* () const { return m_pstruct[m_nbor_list_ptr[m_index]];  }

        AMREX_GPU_HOST_DEVICE
        unsigned int index () const { return m_nbor_list_ptr[m_index]; }

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
        bool operator!= (const_iterator const& /*rhs*/) const { return m_index < m_stop; }

        AMREX_GPU_HOST_DEVICE
        const ParticleType& operator* () const { return m_pstruct[m_nbor_list_ptr[m_index]];  }

        AMREX_GPU_HOST_DEVICE
        unsigned int index () const { return m_nbor_list_ptr[m_index]; }

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

template<typename A, typename B,
         typename std::enable_if<std::is_same<std::remove_cv_t<A>,
                                              std::remove_cv_t<B> >::value, int>::type = 0>
bool isSame (A const* pa, B const* pb)
{
    return pa == pb;
}

template<typename A, typename B,
         typename std::enable_if<!std::is_same<std::remove_cv_t<A>,
                                               std::remove_cv_t<B> >::value, int>::type = 0>
bool isSame (A const* /*pa*/, B const* /*pb*/)
{
    return false;
}

template <class ParticleType>
class NeighborList
{
public:

    template <class PTile, class CheckPair>
    void build (PTile& ptile,
                const amrex::Box& bx, const amrex::Geometry& geom,
                CheckPair&& check_pair, int num_cells=1)
    {
        Gpu::DeviceVector<int> off_bins_v;
        Gpu::DeviceVector<Dim3>      lo_v;
        Gpu::DeviceVector<Dim3>      hi_v;
        Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> dxi_v;
        Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> plo_v;
        off_bins_v.push_back(0);
        off_bins_v.push_back(bx.numPts());
        lo_v.push_back(lbound(bx));
        hi_v.push_back(ubound(bx));
        dxi_v.push_back(geom.InvCellSizeArray());
        plo_v.push_back(geom.ProbLoArray());

        build(ptile, ptile, check_pair, off_bins_v, dxi_v, plo_v,
              lo_v, hi_v, num_cells, 1, nullptr );
    }

    template <class PTile, class CheckPair>
    void build (PTile& ptile,
                CheckPair&& check_pair,
                const Gpu::DeviceVector<int>& off_bins_v,
                const Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>>& dxi_v,
                const Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>>& plo_v,
                const Gpu::DeviceVector<Dim3>& lo_v,
                const Gpu::DeviceVector<Dim3>& hi_v,
                int  num_cells=1,
                int  num_bin_types=1,
                int* bin_type_array=nullptr)
    {
        build(ptile, ptile, check_pair, off_bins_v, dxi_v, plo_v,
              lo_v, hi_v, num_cells, num_bin_types, bin_type_array );
    }

    template <class SrcTile, class TargetTile, class CheckPair>
    void build (SrcTile& src_tile,
                TargetTile& target_tile,
                CheckPair&& check_pair,
                const Gpu::DeviceVector<int>& off_bins_v,
                const Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>>& dxi_v,
                const Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>>& plo_v,
                const Gpu::DeviceVector<Dim3>& lo_v,
                const Gpu::DeviceVector<Dim3>& hi_v,
                int  num_cells=1,
                int  num_bin_types=1,
                int* bin_type_array=nullptr)
    {
        BL_PROFILE("NeighborList::build()");

        bool is_same = isSame(&src_tile, &target_tile);


        // Bin particles to their respective grid(s)
        //---------------------------------------------------------------------------------------------------------
        auto& aos = target_tile.GetArrayOfStructs();
        const auto dst_ptile_data = target_tile.getConstParticleTileData();

        m_pstruct = aos().dataPtr();
        const auto pstruct_ptr = aos().dataPtr();

        const size_t np_total = aos.size();
        const size_t np_real  = src_tile.numRealParticles();

        auto off_bins_p = off_bins_v.data();
        auto dxi_p      = dxi_v.data();
        auto plo_p      = plo_v.data();
        auto  lo_p      =  lo_v.data();
        auto  hi_p      =  hi_v.data();
        BinMapper bm(off_bins_p, dxi_p, plo_p, lo_p, hi_p, bin_type_array);

        // Get tot bin count on host
        int tot_bins;
#ifdef AMREX_USE_GPU
        Gpu::dtoh_memcpy( &tot_bins, off_bins_p + num_bin_types, sizeof(int) );
#else
        std::memcpy( &tot_bins, off_bins_p + num_bin_types, sizeof(int) );
#endif

        m_bins.build(np_total, pstruct_ptr, tot_bins, bm);


        // First pass: count the number of neighbors for each particle
        //---------------------------------------------------------------------------------------------------------
        const size_t np_size  = (num_bin_types > 1) ? np_total : np_real;
        m_nbor_counts.resize( np_size+1, 0);
        m_nbor_offsets.resize(np_size+1);

        auto pnbor_counts = m_nbor_counts.dataPtr();
        auto pnbor_offset = m_nbor_offsets.dataPtr();

        auto pperm   = m_bins.permutationPtr();
        auto poffset = m_bins.offsetsPtr();

        const auto src_ptile_data  = src_tile.getConstParticleTileData();
        const auto src_pstruct_ptr = src_tile.GetArrayOfStructs()().dataPtr();

        AMREX_FOR_1D ( np_size, i,
        {
            int count    = 0;
            int type_i   = bin_type_array ? bin_type_array[i] : 0;
            bool ghost_i = (i >= np_real);
            for (int type(type_i); type<num_bin_types; ++type) {
                int off_bins = off_bins_p[type];

              IntVect iv(AMREX_D_DECL(
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(0)-plo_p[type][0])*dxi_p[type][0])) - lo_p[type].x,
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(1)-plo_p[type][1])*dxi_p[type][1])) - lo_p[type].y,
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(2)-plo_p[type][2])*dxi_p[type][2])) - lo_p[type].z));
              auto iv3 = iv.dim3();

              int ix = iv3.x;
              int iy = iv3.y;
              int iz = iv3.z;

              int nx = hi_p[type].x-lo_p[type].x+1;
              int ny = hi_p[type].y-lo_p[type].y+1;
              int nz = hi_p[type].z-lo_p[type].z+1;

              for (int ii = amrex::max(ix-num_cells, 0); ii <= amrex::min(ix+num_cells, nx-1); ++ii) {
                for (int jj = amrex::max(iy-num_cells, 0); jj <= amrex::min(iy+num_cells, ny-1); ++jj) {
                  for (int kk = amrex::max(iz-num_cells, 0); kk <= amrex::min(iz+num_cells, nz-1); ++kk) {
                    int index = (ii * ny + jj) * nz + kk + off_bins;
                    for (auto p = poffset[index]; p < poffset[index+1]; ++p) {
                      const auto& pid = pperm[p];
                      bool  ghost_pid = (pid >= np_real);
                      if (is_same && (pid == i)) continue;
                      if (call_check_pair(check_pair,
                                          src_ptile_data, dst_ptile_data,
                                          i, pid, type, ghost_i, ghost_pid)) {
                          count += 1;
                      }
                    } // p
                  } // kk
                } // jj
              } // ii
            } // type

            pnbor_counts[i] = count;
        });


        // Second-pass: build the offsets (partial sums) and neighbor list
        //--------------------------------------------------------------------------------------------------------
        Gpu::exclusive_scan(m_nbor_counts.begin(), m_nbor_counts.end(), m_nbor_offsets.begin());

        // Now we can allocate and build our neighbor list
        unsigned int total_nbors;
#ifdef AMREX_USE_GPU
        Gpu::dtoh_memcpy(&total_nbors,m_nbor_offsets.dataPtr()+np_size,sizeof(unsigned int));
#else
        std::memcpy(&total_nbors,m_nbor_offsets.dataPtr()+np_size,sizeof(unsigned int));
#endif

        m_nbor_list.resize(total_nbors);
        auto pm_nbor_list = m_nbor_list.dataPtr();

        AMREX_FOR_1D ( np_size, i,
        {
          int n        = 0;
          int type_i   = bin_type_array ? bin_type_array[i] : 0;
          bool ghost_i = (i >= np_real);
          for (int type(type_i); type<num_bin_types; ++type) {
              int off_bins = off_bins_p[type];

            IntVect iv(AMREX_D_DECL(
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(0)-plo_p[type][0])*dxi_p[type][0])) - lo_p[type].x,
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(1)-plo_p[type][1])*dxi_p[type][1])) - lo_p[type].y,
                static_cast<int>(amrex::Math::floor((src_pstruct_ptr[i].pos(2)-plo_p[type][2])*dxi_p[type][2])) - lo_p[type].z));
            auto iv3 = iv.dim3();

            int ix = iv3.x;
            int iy = iv3.y;
            int iz = iv3.z;

            int nx = hi_p[type].x-lo_p[type].x+1;
            int ny = hi_p[type].y-lo_p[type].y+1;
            int nz = hi_p[type].z-lo_p[type].z+1;

            for (int ii = amrex::max(ix-num_cells, 0); ii <= amrex::min(ix+num_cells, nx-1); ++ii) {
              for (int jj = amrex::max(iy-num_cells, 0); jj <= amrex::min(iy+num_cells, ny-1); ++jj) {
                for (int kk = amrex::max(iz-num_cells, 0); kk <= amrex::min(iz+num_cells, nz-1); ++kk) {
                  int index = (ii * ny + jj) * nz + kk + off_bins;
                  for (auto p = poffset[index]; p < poffset[index+1]; ++p) {
                    const auto& pid = pperm[p];
                    bool  ghost_pid = (pid >= np_real);
                    if (is_same && (pid == i)) continue;
                    if (call_check_pair(check_pair,
                                        src_ptile_data, dst_ptile_data,
                                        i, pid, type, ghost_i, ghost_pid)) {
                        pm_nbor_list[pnbor_offset[i] + n] = pid;
                        ++n;
                    }
                  }// p
                }// kk
              }// jj
            }// ii
          } // type
        });
        Gpu::Device::streamSynchronize();
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

        Gpu::copyAsync(Gpu::deviceToHost, m_nbor_offsets.begin(), m_nbor_offsets.end(), host_nbor_offsets.begin());
        Gpu::copyAsync(Gpu::deviceToHost, m_nbor_list.begin(), m_nbor_list.end(), host_nbor_list.begin());
        Gpu::streamSynchronize();

        for (int i = 0; i < numParticles(); ++i) {
            amrex::Print() << "Particle " << i << " could collide with: ";
            for (unsigned int j = host_nbor_offsets[i]; j < host_nbor_offsets[i+1]; ++j) {
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
