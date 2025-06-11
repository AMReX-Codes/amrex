#ifndef AMREX_STRUCTOFARRAYS_H_
#define AMREX_STRUCTOFARRAYS_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <algorithm>
#include <array>
#include <string>
#include <vector>


namespace amrex {

template <int NReal, int NInt,
          template<class> class Allocator=DefaultAllocator,
          bool use64BitIdCpu=false>
struct StructOfArrays {

    using IdCPU = amrex::PODVector<uint64_t, Allocator<uint64_t> >;
    using RealVector = amrex::PODVector<ParticleReal, Allocator<ParticleReal> >;
    using IntVector = amrex::PODVector<int, Allocator<int> >;

    void define (
        int a_num_runtime_real,
        int a_num_runtime_int,
        std::vector<std::string>* soa_rdata_names=nullptr,
        std::vector<std::string>* soa_idata_names=nullptr
    )
    {
        m_defined = true;
        m_runtime_rdata.resize(a_num_runtime_real);
        m_runtime_idata.resize(a_num_runtime_int );
        m_rdata_names = soa_rdata_names;
        m_idata_names = soa_idata_names;
    }

    [[nodiscard]] int NumRealComps () const noexcept { return NReal + m_runtime_rdata.size(); }

    [[nodiscard]] int NumIntComps () const noexcept { return NInt + m_runtime_idata.size(); }

    [[nodiscard]] IdCPU& GetIdCPUData () { return m_idcpu; }
    [[nodiscard]] std::array<RealVector, NReal>& GetRealData () { return m_rdata; }
    [[nodiscard]] std::array< IntVector,  NInt>& GetIntData  () { return m_idata; }

    /** Get access to the particle id/cpu Arrays */
    [[nodiscard]] const IdCPU& GetIdCPUData () const { return m_idcpu; }
    /** Get access to the particle Real Arrays (only compile-time components) */
    [[nodiscard]] const std::array<RealVector, NReal>& GetRealData () const { return m_rdata; }
    /** Get access to the particle Int Arrays (only compile-time components) */
    [[nodiscard]] const std::array< IntVector,  NInt>& GetIntData  () const { return m_idata; }

    /** Get the names for the real SoA components **/
    [[nodiscard]] std::vector<std::string> GetRealNames () const
    {
        if (m_rdata_names) {
            return *m_rdata_names;
        }
        else {
            return std::vector<std::string>();
        }
    }

    /** Get the names for the int SoA components **/
    [[nodiscard]] std::vector<std::string> GetIntNames () const
    {
        if (m_idata_names) {
            return *m_idata_names;
        }
        else {
            return std::vector<std::string>();
        }
    }

    /** Get access to a particle Real component Array (compile-time and runtime component)
     *
     * @param index component with 0...NReal-1 compile-time and NReal... runtime arguments
     */
    [[nodiscard]] RealVector& GetRealData (const int index) {
        AMREX_ASSERT(index >= 0 && index < NReal + static_cast<int>(m_runtime_rdata.size()));
        if constexpr (NReal == 0) {
            AMREX_ASSERT(m_defined);
            return m_runtime_rdata[index];
        } else {
            if (index < NReal) {
                return m_rdata[index];
            } else {
                AMREX_ASSERT(m_defined);
                return m_runtime_rdata[index - NReal];
            }
        }
    }

    /** Get access to a particle Real component Array (compile-time and runtime component)
     *
     * @param index component with 0...NReal-1 compile-time and NReal... runtime arguments
     */
    [[nodiscard]] const RealVector& GetRealData (const int index) const {
        AMREX_ASSERT(index >= 0 && index < NReal + static_cast<int>(m_runtime_rdata.size()));
        if constexpr (NReal == 0) {
            AMREX_ASSERT(m_defined);
            return m_runtime_rdata[index];
        } else {
            if (index < NReal) {
                return m_rdata[index];
            } else {
                AMREX_ASSERT(m_defined);
                return m_runtime_rdata[index - NReal];
            }
        }
    }

    /** Get access to a particle Real component Array (compile-time and runtime component)
     *
     * @param name named component component with 0...NReal-1 compile-time and NReal... runtime arguments
     */
    [[nodiscard]] RealVector& GetRealData (std::string const & name) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_rdata_names != nullptr, "SoA Real names were not defined.");
        auto const pos = std::find(m_rdata_names->begin(), m_rdata_names->end(), name);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos != m_rdata_names->end(), "Soa Real name='" + name + "' was not found components");

        int const index = std::distance(m_rdata_names->begin(), pos);
        return GetRealData(index);
    }

    /** Get access to a particle Real component Array (compile-time and runtime component)
     *
     * @param name named component component with 0...NReal-1 compile-time and NReal... runtime arguments
     */
    [[nodiscard]] const RealVector& GetRealData (std::string const & name) const {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_rdata_names != nullptr, "SoA Real names were not defined.");
        auto const pos = std::find(m_rdata_names->begin(), m_rdata_names->end(), name);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos != m_rdata_names->end(), "Soa Real name='" + name + "' was not found components");

        int const index = std::distance(m_rdata_names->begin(), pos);
        return GetRealData(index);
    }

    /** Get access to a particle Int component Array (compile-time and runtime component)
     *
     * @param index component with 0...NInt-1 compile-time and NInt... runtime arguments
     */
    [[nodiscard]] IntVector& GetIntData (const int index) {
        AMREX_ASSERT(index >= 0 && index < NInt + static_cast<int>(m_runtime_idata.size()));
        if constexpr (NInt == 0) {
            AMREX_ASSERT(m_defined);
            return m_runtime_idata[index];
        } else {
            if (index < NInt) {
                return m_idata[index];
            } else {
                AMREX_ASSERT(m_defined);
                return m_runtime_idata[index - NInt];
            }
        }
   }

   /** Get access to a particle Int component Array (compile-time and runtime component)
    *
    * @param index component with 0...NInt-1 compile-time and NInt... runtime arguments
    * @return
    */
   [[nodiscard]] const IntVector& GetIntData (const int index) const {
        AMREX_ASSERT(index >= 0 && index < NInt + static_cast<int>(m_runtime_idata.size()));
        if constexpr (NInt == 0) {
            AMREX_ASSERT(m_defined);
            return m_runtime_idata[index];
        } else {
            if (index < NInt) {
                return m_idata[index];
            } else {
                AMREX_ASSERT(m_defined);
                return m_runtime_idata[index - NInt];
            }
        }
    }

    /** Get access to a particle Int component Array (compile-time and runtime component)
     *
     * @param index component with 0...NInt-1 compile-time and NInt... runtime arguments
     * @return
     */
    [[nodiscard]] IntVector& GetIntData (std::string const & name) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_idata_names != nullptr, "SoA Int names were not defined.");
        auto const pos = std::find(m_idata_names->begin(), m_idata_names->end(), name);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos != m_idata_names->end(), "Soa Int name='" + name + "' was not found components");

        int const index = std::distance(m_idata_names->begin(), pos);
        return GetIntData(index);
    }

    /** Get access to a particle Int component Array (compile-time and runtime component)
     *
     * @param index component with 0...NInt-1 compile-time and NInt... runtime arguments
     * @return
     */
    [[nodiscard]] const IntVector& GetIntData (std::string const & name) const {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_idata_names != nullptr, "SoA Int names were not defined.");
        auto const pos = std::find(m_idata_names->begin(), m_idata_names->end(), name);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos != m_idata_names->end(), "Soa Int name='" + name + "' was not found components");

        int const index = std::distance(m_idata_names->begin(), pos);
        return GetIntData(index);
    }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    [[nodiscard]] std::size_t size () const
    {
        if constexpr (use64BitIdCpu == true) {
            return m_idcpu.size();
        } else if constexpr (NReal > 0) {
            return m_rdata[0].size();
        } else if constexpr (NInt > 0) {
            return m_idata[0].size();
        } else {
            if (!m_runtime_rdata.empty()) {
                return m_runtime_rdata[0].size();
            } else if (!m_runtime_idata.empty()) {
                return m_runtime_idata[0].size();
            } else {
                return 0;
            }
        }
    }

    /**
    * \brief Returns whether the SoA is empty (i.e. has size() == 0)
    *
    */
    [[nodiscard]] bool empty () const { return this->size() == 0; }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    [[nodiscard]] int numParticles () const { return numRealParticles(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    [[nodiscard]] int numRealParticles () const { return numTotalParticles()-m_num_neighbor_particles; }

    /**
    * \brief Returns the number of neighbor particles (excluding reals)
    *
    */
    [[nodiscard]] int numNeighborParticles () const { return m_num_neighbor_particles; }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    [[nodiscard]] int numTotalParticles () const { return size(); }

    void setNumNeighbors (int num_neighbors)
    {
        auto nrp = numRealParticles();
        m_num_neighbor_particles = num_neighbors;
        resize(nrp + num_neighbors);
    }

    [[nodiscard]] int getNumNeighbors () const { return m_num_neighbor_particles; }

    void resize (size_t count)
    {
        if constexpr (use64BitIdCpu == true) {
            m_idcpu.resize(count);
        }
        if constexpr (NReal > 0) {
            for (int i = 0; i < NReal; ++i) { m_rdata[i].resize(count); }
        }
        if constexpr (NInt > 0) {
            for (int i = 0; i < NInt;  ++i) { m_idata[i].resize(count); }
        }
        for (int i = 0; i < int(m_runtime_rdata.size()); ++i) { m_runtime_rdata[i].resize(count); }
        for (int i = 0; i < int(m_runtime_idata.size()); ++i) { m_runtime_idata[i].resize(count); }
    }

    [[nodiscard]] uint64_t* idcpuarray () {
        if constexpr (use64BitIdCpu == true) {
            return m_idcpu.dataPtr();
        } else {
            return nullptr;
        }
    }

    [[nodiscard]] GpuArray<ParticleReal*, NReal> realarray ()
    {
        GpuArray<Real*, NReal> arr;
        for (int i = 0; i < NReal; ++i)
        {
            arr[i] = m_rdata[i].dataPtr();
        }
        return arr;
    }

    [[nodiscard]] GpuArray<int*, NInt> intarray ()
    {
        GpuArray<int*, NInt> arr;
        for (int i = 0; i < NInt; ++i)
        {
            arr[i] = m_idata[i].dataPtr();
        }
        return arr;
    }

    int m_num_neighbor_particles{0};

private:
    // compile-time data
    IdCPU m_idcpu;
    std::array<RealVector, NReal> m_rdata;
    std::array< IntVector,  NInt> m_idata;

    // runtime data
    std::vector<RealVector> m_runtime_rdata;
    std::vector<IntVector > m_runtime_idata;

    // names of both compile-time and runtime Real and Int data
    std::vector<std::string>* m_rdata_names = nullptr;
    std::vector<std::string>* m_idata_names = nullptr;

    //! whether the runtime components are sized correctly
    bool m_defined{false};
};

} // namespace amrex

#endif // AMREX_STRUCTOFARRAYS_H_
