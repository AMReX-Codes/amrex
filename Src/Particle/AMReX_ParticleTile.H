#ifndef AMREX_PARTICLETILE_H_
#define AMREX_PARTICLETILE_H_

#include <AMReX_Particle.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>
#include <AMReX_IndexSequence.H>

#include <tuple>
#include <array>
#include <cassert>

#ifdef AMREX_USE_CUDA
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_ptr.h>
#endif

namespace amrex {

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTileData
{
    using ParticleType = Particle<NStructReal, NStructInt>;
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;

    long m_size;
    ParticleType* AMREX_RESTRICT m_aos;
    GpuArray<ParticleReal* AMREX_RESTRICT, NArrayReal> m_rdata;
    GpuArray<int* AMREX_RESTRICT, NArrayInt> m_idata;    

    long m_psize;
    int m_num_runtime_real;
    int m_num_runtime_int;
    ParticleReal** m_runtime_r_data;
    int** m_runtime_i_data;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void packParticleData(char* buffer, int src_index, int dst_index) const noexcept
    {
        AMREX_ASSERT(src_index < m_size);
        auto dst = buffer + dst_index*m_psize;
        memcpy(dst, m_aos + src_index, sizeof(ParticleType));
        dst += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            memcpy(dst, m_rdata[i] + src_index, sizeof(ParticleReal));
            dst += sizeof(ParticleReal);
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            memcpy(dst, m_runtime_r_data[i] + src_index, sizeof(ParticleReal));
            dst += sizeof(ParticleReal);
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            memcpy(dst, m_idata[i] + src_index, sizeof(int));
            dst += sizeof(int);
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            memcpy(dst, m_runtime_i_data[i] + src_index, sizeof(int));
            dst += sizeof(int);
        }        
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void unpackParticleData(const char* buffer, int src_index, int dst_index) const noexcept
    {
        AMREX_ASSERT(dst_index < m_size);
        auto src = buffer + src_index*m_psize;
        memcpy(m_aos + dst_index, src, sizeof(ParticleType));
        src += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            memcpy(m_rdata[i] + dst_index, src, sizeof(ParticleReal));
            src += sizeof(ParticleReal);
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            memcpy(m_runtime_r_data[i] + dst_index, src, sizeof(ParticleReal));
            src += sizeof(ParticleReal);
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            memcpy(m_idata[i] + dst_index, src, sizeof(int));
            src += sizeof(int);
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            memcpy(m_runtime_i_data[i] + dst_index, src, sizeof(int));
            src += sizeof(int);
        }
    }
};

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ConstParticleTileData
{
    using ParticleType = Particle<NStructReal, NStructInt>;
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;

    long m_size;
    const ParticleType* AMREX_RESTRICT m_aos;
    GpuArray<const ParticleReal* AMREX_RESTRICT, NArrayReal> m_rdata;
    GpuArray<const int* AMREX_RESTRICT, NArrayInt > m_idata;    

    long m_psize;
    int m_num_runtime_real;
    int m_num_runtime_int;
    ParticleReal* const* m_runtime_r_data;
    int* const* m_runtime_i_data;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void packParticleData(char* buffer, int src_index, int dst_index) const noexcept
    {
        AMREX_ASSERT(src_index < m_size);
        auto dst = buffer + dst_index*m_psize;
        memcpy(dst, m_aos + src_index, sizeof(ParticleType));
        dst += sizeof(ParticleType);
        for (int i = 0; i < NArrayReal; ++i)
        {
            memcpy(dst, m_rdata[i] + src_index, sizeof(ParticleReal));
            dst += sizeof(ParticleReal);
        }
        for (int i = 0; i < m_num_runtime_real; ++i)
        {
            memcpy(dst, m_runtime_r_data[i] + src_index, sizeof(ParticleReal));
            dst += sizeof(ParticleReal);
        }
        for (int i = 0; i < NArrayInt; ++i)
        {
            memcpy(dst, m_idata[i] + src_index, sizeof(int));
            dst += sizeof(int);
        }
        for (int i = 0; i < m_num_runtime_int; ++i)
        {
            memcpy(dst, m_runtime_i_data[i] + src_index, sizeof(int));
            dst += sizeof(int);
        }        
    }
};

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTile
{
    using ParticleType = Particle<NStructReal, NStructInt>;

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
        for (int i = 0; i < a_num_runtime_real; ++i) 
            m_runtime_r_ptrs[i] = m_soa_tile.GetRealData(NArrayReal + i).dataPtr();
        for (int i = 0; i < a_num_runtime_int; ++i) 
            m_runtime_i_ptrs[i] = m_soa_tile.GetIntData(NArrayInt + i).dataPtr();
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

    void resize(std::size_t count)
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

    ParticleTileDataType getParticleTileData ()
    {
        ParticleTileDataType ptd;
        ptd.m_aos = m_aos_tile().dataPtr();
        for (int i = 0; i < NArrayReal; ++i)
            ptd.m_rdata[i] = m_soa_tile.GetRealData(i).dataPtr();
        for (int i = 0; i < NArrayInt; ++i)
            ptd.m_idata[i] = m_soa_tile.GetIntData(i).dataPtr();
        ptd.m_size = numParticles();
        ptd.m_num_runtime_real = m_runtime_r_ptrs.size();
        ptd.m_num_runtime_int = m_runtime_i_ptrs.size();
        ptd.m_runtime_r_data = m_runtime_r_ptrs.dataPtr();
        ptd.m_runtime_i_data = m_runtime_i_ptrs.dataPtr();
        ptd.m_psize = sizeof(ParticleType)
            + (NArrayReal + ptd.m_num_runtime_real)*sizeof(ParticleReal) 
            + (NArrayInt + ptd.m_num_runtime_int)*sizeof(int);
        return ptd;
    }

    ConstParticleTileDataType getConstParticleTileData () const
    {
        ConstParticleTileDataType ptd;
        ptd.m_aos = m_aos_tile().dataPtr();
        for (int i = 0; i < NArrayReal; ++i)
            ptd.m_rdata[i] = m_soa_tile.GetRealData(i).dataPtr();
        for (int i = 0; i < NArrayInt; ++i)
            ptd.m_idata[i] = m_soa_tile.GetIntData(i).dataPtr();
        ptd.m_size = numParticles();
        ptd.m_num_runtime_real = m_runtime_r_ptrs.size();
        ptd.m_num_runtime_int = m_runtime_i_ptrs.size();
        ptd.m_runtime_r_data = m_runtime_r_ptrs.dataPtr();
        ptd.m_runtime_i_data = m_runtime_i_ptrs.dataPtr();
        ptd.m_psize = sizeof(ParticleType)
            + (NArrayReal + ptd.m_num_runtime_real)*sizeof(ParticleReal) 
            + (NArrayInt  + ptd.m_num_runtime_int)*sizeof(int);
        return ptd;
    }

private:

    AoS m_aos_tile;
    SoA m_soa_tile;

    bool m_defined;

#ifdef AMREX_USE_CUDA
    Gpu::DeviceVector<ParticleReal*> m_runtime_r_ptrs;
    Gpu::DeviceVector<int*> m_runtime_i_ptrs;    
#endif
};

} // namespace amrex;

#endif // AMREX_PARTICLETILE_H_
