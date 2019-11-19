#ifndef AMREX_PARTICLETILE_H_
#define AMREX_PARTICLETILE_H_

#include <AMReX_Particle.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>

#include <array>

namespace amrex {

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTileData
{
    static constexpr int NAR = NArrayReal;
    static constexpr int NAI = NArrayInt;    
    using ParticleType = Particle<NStructReal, NStructInt>;
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;

    long m_size;
    ParticleType* AMREX_RESTRICT m_aos;
    GpuArray<ParticleReal* AMREX_RESTRICT, NArrayReal> m_rdata;
    GpuArray<int* AMREX_RESTRICT, NArrayInt> m_idata;    

    int m_num_runtime_real;
    int m_num_runtime_int;
    ParticleReal* AMREX_RESTRICT * AMREX_RESTRICT m_runtime_rdata;
    int* AMREX_RESTRICT * AMREX_RESTRICT m_runtime_idata;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void packParticleData (char* buffer, int src_index, int dst_index,
                           const int* comm_real, const int * comm_int, long psize) const noexcept
    {
        AMREX_ASSERT(src_index < m_size);
        auto dst = buffer + dst_index*psize;
        memcpy(dst, m_aos + src_index, sizeof(ParticleType));
        dst += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            if (comm_real[i])
            {
                memcpy(dst, m_rdata[i] + src_index, sizeof(ParticleReal));
                dst += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            if (comm_real[NArrayReal+i])
            {
                memcpy(dst, m_runtime_rdata[i] + src_index, sizeof(ParticleReal));
                dst += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            if (comm_int[i])
            {
                memcpy(dst, m_idata[i] + src_index, sizeof(int));
                dst += sizeof(int);
            }
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            if (comm_int[NArrayInt+i])
            {
                memcpy(dst, m_runtime_idata[i] + src_index, sizeof(int));
                dst += sizeof(int);
            }
        }        
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void unpackParticleData (const char* buffer, int src_index, int dst_index,
                             const int* comm_real, const int * comm_int, long psize) const noexcept
    {
        AMREX_ASSERT(dst_index < m_size);
        auto src = buffer + src_index*psize;
        memcpy(m_aos + dst_index, src, sizeof(ParticleType));
        src += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            if (comm_real[i])
            {
                memcpy(m_rdata[i] + dst_index, src, sizeof(ParticleReal));
                src += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            if (comm_real[NArrayReal+i])
            {
                memcpy(m_runtime_rdata[i] + dst_index, src, sizeof(ParticleReal));
                src += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            if (comm_int[i])
            {
                memcpy(m_idata[i] + dst_index, src, sizeof(int));
                src += sizeof(int);
            }
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            if (comm_int[NArrayInt+i])
            {
                memcpy(m_runtime_idata[i] + dst_index, src, sizeof(int));
                src += sizeof(int);
            }
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SuperParticleType getSuperParticle (int index) const noexcept
    {
        AMREX_ASSERT(index < m_size);
        SuperParticleType sp;
        for (int i = 0; i < NStructReal+AMREX_SPACEDIM; ++i)
            sp.m_rdata.arr[i] = m_aos[index].m_rdata.arr[i];
        for (int i = 0; i < NArrayReal; ++i)
            sp.m_rdata.arr[NStructReal+AMREX_SPACEDIM+i] = m_rdata[i][index];
        for (int i = 0; i < NStructInt+2; ++i)
            sp.m_idata.arr[i] = m_aos[index].m_idata.arr[i];
        for (int i = 0; i < NArrayInt; ++i)
            sp.m_idata.arr[NStructInt+2+i] = m_idata[i][index];
        return sp;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setSuperParticle (const SuperParticleType& sp, int index) const noexcept
    {
        for (int i = 0; i < NStructReal+AMREX_SPACEDIM; ++i)
            m_aos[index].m_rdata.arr[i] = sp.m_rdata.arr[i];
        for (int i = 0; i < NArrayReal; ++i)
            m_rdata[i][index] = sp.m_rdata.arr[NStructReal+AMREX_SPACEDIM+i];
        for (int i = 0; i < NStructInt+2; ++i)
            m_aos[index].m_idata.arr[i] = sp.m_idata.arr[i];
        for (int i = 0; i < NArrayInt; ++i)
            m_idata[i][index] = sp.m_idata.arr[NStructInt+2+i];
    }
};

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ConstParticleTileData
{
    static constexpr int NAR = NArrayReal;
    static constexpr int NAI = NArrayInt;
    using ParticleType = Particle<NStructReal, NStructInt>;
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;

    long m_size;
    const ParticleType* AMREX_RESTRICT m_aos;
    GpuArray<const ParticleReal* AMREX_RESTRICT, NArrayReal> m_rdata;
    GpuArray<const int* AMREX_RESTRICT, NArrayInt > m_idata;    

    int m_num_runtime_real;
    int m_num_runtime_int;
    const ParticleReal* AMREX_RESTRICT * AMREX_RESTRICT m_runtime_rdata;
    const int* AMREX_RESTRICT * AMREX_RESTRICT m_runtime_idata;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void packParticleData(char* buffer, int src_index, int dst_index,
                          const int* comm_real, const int * comm_int, long psize) const noexcept
    {
        AMREX_ASSERT(src_index < m_size);
        auto dst = buffer + dst_index*psize;
        memcpy(dst, m_aos + src_index, sizeof(ParticleType));
        dst += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            if (comm_real[i])
            {
                memcpy(dst, m_rdata[i] + src_index, sizeof(ParticleReal));
                dst += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            if (comm_real[NArrayReal+i])
            {
                memcpy(dst, m_runtime_rdata[i] + src_index, sizeof(ParticleReal));
                dst += sizeof(ParticleReal);
            }
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            if (comm_int[i])
            {
                memcpy(dst, m_idata[i] + src_index, sizeof(int));
                dst += sizeof(int);
            }
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            if (comm_int[NArrayInt+i])
            {
                memcpy(dst, m_runtime_idata[i] + src_index, sizeof(int));
                dst += sizeof(int);
            }
        }        
    }

    AMREX_GPU_HOST_DEVICE
    SuperParticleType getSuperParticle (int index) const
    {
        AMREX_ASSERT(index < m_size);
        SuperParticleType sp;
        for (int i = 0; i < NStructReal+AMREX_SPACEDIM; ++i)
            sp.m_rdata.arr[i] = m_aos[index].m_rdata.arr[i];
        for (int i = 0; i < NArrayReal; ++i)
            sp.m_rdata.arr[NStructReal+AMREX_SPACEDIM+i] = m_rdata[i][index];
        for (int i = 0; i < NStructInt+2; ++i)
            sp.m_idata.arr[i] = m_aos[index].m_idata.arr[i];
        for (int i = 0; i < NArrayInt; ++i)
            sp.m_idata.arr[NStructInt+2+i] = m_idata[i][index];
        return sp;
    }
};

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTile
{
    using ParticleType = Particle<NStructReal, NStructInt>;
    static constexpr int NAR = NArrayReal;
    static constexpr int NAI = NArrayInt;

    using AoS = ArrayOfStructs<NStructReal, NStructInt>;
    using ParticleVector = typename AoS::ParticleVector;

    using SoA = StructOfArrays<NArrayReal, NArrayInt>;
    using RealVector = typename SoA::RealVector;
    using IntVector = typename SoA::IntVector;

    using ParticleTileDataType = ParticleTileData<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ConstParticleTileDataType = ConstParticleTileData<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    ParticleTile()
        : m_defined(false)
        {}

    void define (int a_num_runtime_real, int a_num_runtime_int)
    {
        m_defined = true;
        GetStructOfArrays().define(a_num_runtime_real, a_num_runtime_int);
        m_runtime_r_ptrs.resize(a_num_runtime_real);
        m_runtime_i_ptrs.resize(a_num_runtime_int);
        m_runtime_r_cptrs.resize(a_num_runtime_real);
        m_runtime_i_cptrs.resize(a_num_runtime_int);
    }
    
    AoS&       GetArrayOfStructs ()       { return m_aos_tile; }
    const AoS& GetArrayOfStructs () const { return m_aos_tile; }

    SoA&       GetStructOfArrays ()       { return m_soa_tile; }
    const SoA& GetStructOfArrays () const { return m_soa_tile; }

    bool empty () const { return m_aos_tile.empty(); }
    
    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */

    std::size_t size () const { return m_aos_tile.size(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    int numParticles () const { return m_aos_tile.numParticles(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    int numRealParticles () const { return m_aos_tile.numRealParticles(); }

    /**
    * \brief Returns the number of neighbor particles (excluding reals)
    *
    */
    int numNeighborParticles () const { return m_aos_tile.numNeighborParticles(); }    

    /**
    * \brief Returns the total number of particles, real and neighbor
    *
    */
    int numTotalParticles () const { return m_aos_tile.numTotalParticles() ; }

    void setNumNeighbors (int num_neighbors) 
    {
        m_soa_tile.setNumNeighbors(num_neighbors);
        m_aos_tile.setNumNeighbors(num_neighbors);
    }

    int getNumNeighbors () 
    {
        AMREX_ASSERT( m_soa_tile.getNumNeighbors() == m_aos_tile.getNumNeighbors() );
        return m_aos_tile.getNumNeighbors();
    }

    void resize (std::size_t count)
    {
        m_aos_tile.resize(count);
        m_soa_tile.resize(count);
    }

    ///
    /// Add one particle to this tile.
    ///
    void push_back (const ParticleType& p) { m_aos_tile().push_back(p); }

    ///
    /// Add a Real value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_real (int comp, ParticleReal v) { 
        m_soa_tile.GetRealData(comp).push_back(v);
    }

    ///
    /// Add Real values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_real (const std::array<ParticleReal, NArrayReal>& v) { 
        for (int i = 0; i < NArrayReal; ++i) {
            m_soa_tile.GetRealData(i).push_back(v[i]);
        }
    }

    ///
    /// Add a range of Real values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, const ParticleReal* beg, const ParticleReal* end) {
        auto it = m_soa_tile.GetRealData(comp).end();
        m_soa_tile.GetRealData(comp).insert(it, beg, end);
    }

    ///
    /// Add npar copies of the Real value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, std::size_t npar, ParticleReal v) {
        auto new_size = m_soa_tile.GetRealData(comp).size() + npar;
        m_soa_tile.GetRealData(comp).resize(new_size, v);
    }

    ///
    /// Add an int value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_int (int comp, int v) { 
        m_soa_tile.GetIntData(comp).push_back(v);
    }
    
    ///
    /// Add int values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_int (const std::array<int, NArrayInt>& v) { 
        for (int i = 0; i < NArrayInt; ++i) {
            m_soa_tile.GetIntData(i).push_back(v[i]);
        }
    }
    
    ///
    /// Add a range of int values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, const int* beg, const int* end) {
        auto it = m_soa_tile.GetIntData(comp).end();
        m_soa_tile.GetIntData(comp).insert(it, beg, end);
    }
    
    ///
    /// Add npar copies of the int value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, std::size_t npar, int v) {
        auto new_size = m_soa_tile.GetIntData(comp).size() + npar;
        m_soa_tile.GetIntData(comp).resize(new_size, v);
    }

    int NumRealComps () const noexcept { return m_soa_tile.NumRealComps(); }

    int NumIntComps () const noexcept { return m_soa_tile.NumIntComps(); }

    void swap (ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt>& other)
    {
        m_aos_tile().swap(other.GetArrayOfStructs()());
        for (int j = 0; j < NumRealComps(); ++j)
        {
            auto& rdata = GetStructOfArrays().GetRealData(j);
            rdata.swap(other.GetStructOfArrays().GetRealData(j));
        }
        
        for (int j = 0; j < NumIntComps(); ++j)
        {
            auto& idata = GetStructOfArrays().GetIntData(j);
            idata.swap(other.GetStructOfArrays().GetIntData(j));
        }
    }

    ParticleTileDataType getParticleTileData ()
    {
        for (int i = 0; i < m_runtime_r_ptrs.size(); ++i) {
            m_runtime_r_ptrs[i] = m_soa_tile.GetRealData(NArrayReal + i).dataPtr();
        }
        for (int i = 0; i < m_runtime_i_ptrs.size(); ++i) 
            m_runtime_i_ptrs[i] = m_soa_tile.GetIntData(NArrayInt + i).dataPtr();

        ParticleTileDataType ptd;
        ptd.m_aos = m_aos_tile().dataPtr();
        for (int i = 0; i < NArrayReal; ++i)
            ptd.m_rdata[i] = m_soa_tile.GetRealData(i).dataPtr();
        for (int i = 0; i < NArrayInt; ++i)
            ptd.m_idata[i] = m_soa_tile.GetIntData(i).dataPtr();
        ptd.m_size = size();
        ptd.m_num_runtime_real = m_runtime_r_ptrs.size();
        ptd.m_num_runtime_int = m_runtime_i_ptrs.size();
        ptd.m_runtime_rdata = m_runtime_r_ptrs.dataPtr();
        ptd.m_runtime_idata = m_runtime_i_ptrs.dataPtr();
        return ptd;
    }

    ConstParticleTileDataType getConstParticleTileData () const
    {
        for (int i = 0; i < m_runtime_r_ptrs.size(); ++i) {
            m_runtime_r_cptrs[i] = m_soa_tile.GetRealData(NArrayReal + i).dataPtr();
        }
        for (int i = 0; i < m_runtime_i_ptrs.size(); ++i) 
            m_runtime_i_cptrs[i] = m_soa_tile.GetIntData(NArrayInt + i).dataPtr();

        ConstParticleTileDataType ptd;
        ptd.m_aos = m_aos_tile().dataPtr();
        for (int i = 0; i < NArrayReal; ++i)
            ptd.m_rdata[i] = m_soa_tile.GetRealData(i).dataPtr();
        for (int i = 0; i < NArrayInt; ++i)
            ptd.m_idata[i] = m_soa_tile.GetIntData(i).dataPtr();
        ptd.m_size = size();
        ptd.m_num_runtime_real = m_runtime_r_cptrs.size();
        ptd.m_num_runtime_int = m_runtime_i_cptrs.size();
        ptd.m_runtime_rdata = m_runtime_r_cptrs.dataPtr();
        ptd.m_runtime_idata = m_runtime_i_cptrs.dataPtr();
        return ptd;
    }

private:

    AoS m_aos_tile;
    SoA m_soa_tile;

    bool m_defined;

    Gpu::DeviceVector<ParticleReal*> m_runtime_r_ptrs;
    Gpu::DeviceVector<int*> m_runtime_i_ptrs;

    mutable Gpu::DeviceVector<const ParticleReal*> m_runtime_r_cptrs;
    mutable Gpu::DeviceVector<const int*> m_runtime_i_cptrs;
};

} // namespace amrex;

#endif // AMREX_PARTICLETILE_H_
