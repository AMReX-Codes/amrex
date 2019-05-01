#ifndef AMREX_STRUCTOFARRAYS_H_
#define AMREX_STRUCTOFARRAYS_H_

#include <array>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_CudaContainers.H>

namespace amrex {

template <int NReal, int NInt>
struct StructOfArrays {

    using RealVector = Gpu::ManagedDeviceVector<Real>;
    using IntVector = Gpu::ManagedDeviceVector<int>;    

    StructOfArrays()
        : m_num_neighbor_particles(0),
          m_defined(false)
        {}

    void define(int a_num_runtime_real, int a_num_runtime_int)
    {
        m_defined = true;
        m_runtime_rdata.resize(a_num_runtime_real);
        m_runtime_idata.resize(a_num_runtime_int );
    }
    
    std::array<RealVector, NReal>& GetRealData () { return m_rdata; }
    std::array< IntVector,  NInt>& GetIntData  () { return m_idata; }

    const std::array<RealVector, NReal>& GetRealData () const { return m_rdata; }
    const std::array< IntVector,  NInt>& GetIntData  () const { return m_idata; }

    RealVector& GetRealData (const int index) {
        AMREX_ASSERT(index < NReal + m_runtime_rdata.size());
        if (index < NReal) return m_rdata[index];
        else {
            AMREX_ASSERT(m_defined);
            return m_runtime_rdata[index - NReal];
        }
    }
    
    const RealVector& GetRealData (const int index) const {
        AMREX_ASSERT(index < NReal + m_runtime_rdata.size());
        if (index < NReal) return m_rdata[index];
        else {
            AMREX_ASSERT(m_defined);
            return m_runtime_rdata[index - NReal];
        }
    }

    IntVector& GetIntData (const int index) {
        AMREX_ASSERT(index < NInt + m_runtime_idata.size());
        if (index < NInt) return m_idata[index];
        else {
            AMREX_ASSERT(m_defined);
            return m_runtime_idata[index - NInt];
        }
   }
    
    const IntVector& GetIntData (const int index) const {
        AMREX_ASSERT(index < NInt + m_runtime_idata.size());
        if (index < NInt) return m_idata[index];
        else {
            AMREX_ASSERT(m_defined);
            return m_runtime_idata[index - NInt];
        }
    }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    std::size_t size () const
    {
        if (NReal > 0)
            return m_rdata[0].size();
        else if (NInt > 0) 
            return m_idata[0].size();
        else if (m_runtime_rdata.size() > 0)
            return m_runtime_rdata[0].size();
        else if (m_runtime_idata.size() > 0)
            return m_runtime_idata[0].size();
        else
            return 0;
    }
    
    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    int numParticles () const { return numRealParticles(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    int numRealParticles () const { return numTotalParticles()-m_num_neighbor_particles; }

    /**
    * \brief Returns the number of neighbor particles (excluding reals)
    *
    */
    int numNeighborParticles () const { return m_num_neighbor_particles; }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    int numTotalParticles () const { return size(); }

    void setNumNeighbors (int num_neighbors)
    { 
        auto nrp = numRealParticles();
        m_num_neighbor_particles = num_neighbors;
        resize(nrp + num_neighbors);
    }

    int getNumNeighbors () { return m_num_neighbor_particles; }

    void resize(size_t count)
    {
        for (int i = 0; i < NReal; ++i) m_rdata[i].resize(count);
        for (int i = 0; i < NInt;  ++i) m_idata[i].resize(count); 
        for (int i = 0; i < (int) m_runtime_rdata.size(); ++i) m_runtime_rdata[i].resize(count);
        for (int i = 0; i < (int) m_runtime_idata.size(); ++i) m_runtime_idata[i].resize(count);
    }

    GpuArray<Real*, NReal> realarray ()
    {
        GpuArray<Real*, NReal> arr;
        for (int i = 0; i < NReal; ++i)
        {
            arr[i] = m_rdata[i].dataPtr();
        }
        return arr;
    }

    GpuArray<int*, NInt> intarray ()
    {
        GpuArray<int*, NInt> arr;
        for (int i = 0; i < NInt; ++i)
        {
            arr[i] = m_rdata[i].dataPtr();
        }
        return arr;
    }

    int m_num_neighbor_particles;

private:
    std::array<RealVector, NReal> m_rdata;
    std::array< IntVector,  NInt> m_idata;

    std::vector<RealVector> m_runtime_rdata;
    std::vector<IntVector > m_runtime_idata;
    
    bool m_defined;
};

} // namespace amrex

#endif // AMREX_STRUCTOFARRAYS_H_
