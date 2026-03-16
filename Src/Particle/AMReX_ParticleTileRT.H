#ifndef AMREX_PARTICLETILERT_H_
#define AMREX_PARTICLETILERT_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_Particle.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

#include <array>
#include <string>
#include <type_traits>
#include <vector>


namespace amrex {

// Forward Declaration
template <class RType, class IType>
struct RTSoAParticle;

/// \ingroup amrex_particles

template <class T>
struct ArrayView {
    using size_type = std::size_t;

    T* m_data = nullptr;
    size_type m_capacity = 0;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator[] (const size_type i) const {
        AMREX_ASSERT(i < m_capacity);
        return m_data[i];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* data () const {
        return m_data;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* dataPtr () const {
        return m_data;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* begin () const {
        return m_data;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* end () const {
        return m_data + m_capacity;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    explicit operator T* () const {
        return m_data;
    }
};

template <class T>
ArrayView(T* data, std::size_t capacity) -> ArrayView<T>;

template <class RType, class IType>
struct ParticleTileDataRT
{
    using Self = ParticleTileDataRT<RType, IType>;
    using size_type = std::size_t;

    using ParticleType = RTSoAParticle<RType, IType>;
    using RealType = RType;
    using IntType = IType;

    static constexpr bool is_particle_tile_data = true;
    static constexpr bool is_const = std::is_const_v<RType>;

    size_type m_capacity = 0;
    std::conditional_t<is_const,
        const uint64_t* AMREX_RESTRICT,
        uint64_t* AMREX_RESTRICT> m_idcpu = nullptr;
    RType* AMREX_RESTRICT m_rdata = nullptr;
    IType* AMREX_RESTRICT m_idata = nullptr;
    int m_n_real = 0;
    int m_n_int = 0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr ParticleTileDataRT () = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr ParticleTileDataRT (size_type a_capacity,
        decltype(m_idcpu) a_idcpu, decltype(m_rdata) a_rdata,
        decltype(m_idata) a_idata, int a_n_real, int a_n_int) noexcept :
        m_capacity{a_capacity},
        m_idcpu{a_idcpu},
        m_rdata{a_rdata},
        m_idata{a_idata},
        m_n_real{a_n_real},
        m_n_int{a_n_int} {}

    template <class aRType, class aIType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr ParticleTileDataRT (ParticleTileDataRT<aRType, aIType> const& rhs) noexcept :
        m_capacity{rhs.m_capacity},
        m_idcpu{rhs.m_idcpu},
        m_rdata{rhs.m_rdata},
        m_idata{rhs.m_idata},
        m_n_real{rhs.m_n_real},
        m_n_int{rhs.m_n_int} {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType& pos (const int dir, const size_type index) const &
    {
        AMREX_ASSERT(index < m_capacity && 0 <= dir && dir < m_n_real);
        return this->m_rdata[dir * m_capacity + index];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) id (const size_type index) const &
    {
        AMREX_ASSERT(index < m_capacity);
        if constexpr (is_const) {
            return ConstParticleIDWrapper(this->m_idcpu[index]);
        } else {
            return ParticleIDWrapper(this->m_idcpu[index]);
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) cpu (const size_type index) const &
    {
        AMREX_ASSERT(index < m_capacity);
        if constexpr (is_const) {
            return ConstParticleCPUWrapper(this->m_idcpu[index]);
        } else {
            return ParticleCPUWrapper(this->m_idcpu[index]);
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) idcpu (const size_type index) const &
    {
        AMREX_ASSERT(index < m_capacity);
        return this->m_idcpu[index];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType * rdata (const int comp_index) const &
    {
        AMREX_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return this->m_rdata + comp_index * m_capacity;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IType * idata (const int comp_index) const &
    {
        AMREX_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return this->m_idata + comp_index * m_capacity;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType& rdata (const int comp_index, const size_type partice_index) const &
    {
        AMREX_ASSERT(partice_index < m_capacity && 0 <= comp_index && comp_index < m_n_real);
        return this->m_rdata[comp_index * m_capacity + partice_index];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IType& idata (const int comp_index, const size_type partice_index) const &
    {
        AMREX_ASSERT(partice_index < m_capacity && 0 <= comp_index && comp_index < m_n_int);
        return this->m_idata[comp_index * m_capacity + partice_index];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) operator[] (const size_type index) const
    {
        AMREX_ASSERT(index < m_capacity);
        return RTSoAParticle<RType, IType>(*this, index);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void packParticleData (char* buffer, size_type src_index, std::size_t dst_offset,
                           const int* comm_real, const int * comm_int) const noexcept
    {
        AMREX_ASSERT(src_index < m_capacity);
        auto* dst = buffer + dst_offset;

        memcpy(dst, m_idcpu + src_index, sizeof(uint64_t));
        dst += sizeof(uint64_t);

        for (int i = 0; i < m_n_real; ++i) {
            if (comm_real[i]) {
                memcpy(dst, m_rdata + i * m_capacity + src_index, sizeof(RType));
                dst += sizeof(RType);
            }
        }

        for (int i = 0; i < m_n_int; ++i) {
            if (comm_int[i]) {
                memcpy(dst, m_idata + i * m_capacity + src_index, sizeof(IType));
                dst += sizeof(IType);
            }
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void unpackParticleData (const char* buffer, Long src_offset, size_type dst_index,
                             const int* comm_real, const int* comm_int) const noexcept
    {
        AMREX_ASSERT(dst_index < m_capacity);
        const auto* src = buffer + src_offset;

        memcpy(m_idcpu + dst_index, src, sizeof(uint64_t));
        src += sizeof(uint64_t);

        for (int i = 0; i < m_n_real; ++i) {
            if (comm_real[i]) {
                memcpy(m_rdata + i * m_capacity + dst_index, src, sizeof(RType));
                src += sizeof(RType);
            }
        }

        for (int i = 0; i < m_n_int; ++i) {
            if (comm_int[i]) {
                memcpy(m_idata + i * m_capacity + dst_index, src, sizeof(IType));
                src += sizeof(IType);
            }
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Particle<0, 0> getSuperParticle (size_type index) const noexcept
    {
        return {{
            { AMREX_D_DECL( pos(0, index), pos(1, index), pos(2, index) ) },
            idcpu(index)
        }};
    }
};

// SOA Particle Structure
template <class RType, class IType>
struct RTSoAParticle
{
    using PTD = ParticleTileDataRT<RType, IType>;
    using ConstType = RTSoAParticle<std::add_const_t<RType>, std::add_const_t<IType>>;
    static constexpr bool is_constsoa_particle = PTD::is_const;
    static constexpr int NReal = 0;
    static constexpr int NInt = 0;
    static constexpr bool is_soa_particle = true;
    static constexpr bool is_rtsoa_particle = true;
    using size_type = typename PTD::size_type;
    using RealType = RType;
    using IntType = IType;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RTSoAParticle (PTD const& ptd, size_type i) noexcept :
        m_particle_tile_data(ptd),
        m_index(i) {}

    template <class aRType, class aIType>
    friend struct RTSoAParticle;

    template <class aRType, class aIType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RTSoAParticle (RTSoAParticle<aRType, aIType> const& rhs) noexcept :
        m_particle_tile_data(rhs.m_particle_tile_data),
        m_index(rhs.m_index) {}

    // functions to get id and cpu in the SOA data

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) cpu () const & { return this->m_particle_tile_data.cpu(m_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) id () const & { return this->m_particle_tile_data.id(m_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) idcpu () const & { return this->m_particle_tile_data.idcpu(m_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType& rdata (const int comp_index) const & {
        return this->m_particle_tile_data.rdata(comp_index, m_index);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IType& idata (const int comp_index) const & {
        return this->m_particle_tile_data.idata(comp_index, m_index);
    }

    // functions to get positions of the particle in the SOA data

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect pos () const & {
        return RealVect(AMREX_D_DECL(
            this->m_particle_tile_data.pos(0, m_index),
            this->m_particle_tile_data.pos(1, m_index),
            this->m_particle_tile_data.pos(2, m_index)
        ));
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType& pos (int position_index) const & {
        return this->m_particle_tile_data.pos(position_index, m_index);
    }

    static Long NextID ();

    /**
    * \brief This version can only be used inside omp critical.
    */
    static Long UnprotectedNextID ();

    /**
    * \brief Reset on restart.
    *
    * \param nextid
    */
    static void NextID (Long nextid);

private:

    static_assert(std::is_trivially_copyable<PTD>(), "ParticleTileData is not trivially copyable");

    PTD m_particle_tile_data;
    size_type m_index;
};

struct NextIDRTSoA {
    inline static Long the_next_id = 1;
};

template <class RType, class IType>
Long
RTSoAParticle<RType, IType>::NextID ()
{
    Long next;
// we should be able to test on _OPENMP < 201107 for capture (version 3.1)
// but we must work around a bug in gcc < 4.9
#if defined(AMREX_USE_OMP) && defined(_OPENMP) && _OPENMP < 201307
#pragma omp critical (amrex_particle_nextid)
#elif defined(AMREX_USE_OMP)
#pragma omp atomic capture
#endif
    next = NextIDRTSoA::the_next_id++;

    if (next > LongParticleIds::LastParticleID) {
        amrex::Abort("RTSoAParticle<RType, IType>::NextID() -- too many particles");
    }

    return next;
}

template <class RType, class IType>
Long
RTSoAParticle<RType, IType>::UnprotectedNextID ()
{
    Long next = NextIDRTSoA::the_next_id++;
    if (next > LongParticleIds::LastParticleID) {
        amrex::Abort("RTSoAParticle<RType, IType>::NextID() -- too many particles");
    }
    return next;
}

template <class RType, class IType>
void
RTSoAParticle<RType, IType>::NextID (Long nextid)
{
    NextIDRTSoA::the_next_id = nextid;
}

template<class RType, class IType>
using ConstParticleTileDataRT = ParticleTileDataRT<const RType, const IType>;

template <class RType=ParticleReal, class IType=int>
struct ParticleTileRT
{
    using ParticleType = RTSoAParticle<RType, IType>;
    using RealType = RType;
    using IntType = IType;

    using ParticleTileDataType = ParticleTileDataRT<RType, IType>;
    using ConstParticleTileDataType = ParticleTileDataRT<const RType, const IType>;

    using size_type = typename ParticleTileDataType::size_type;

    using AoS = ThisParticleTileHasNoAoS;
    using SoA = ParticleTileRT<RType, IType>;

    using RealVector = ArrayView<RType>;
    using IntVector = ArrayView<IType>;

    ParticleTileRT () = default;

#ifndef _WIN32  // workaround windows compiler bug
    ~ParticleTileRT () = default;

    ParticleTileRT (ParticleTileRT const&) = delete;
    ParticleTileRT (ParticleTileRT &&) noexcept = default;

    ParticleTileRT& operator= (ParticleTileRT const&) = delete;
    ParticleTileRT& operator= (ParticleTileRT &&) noexcept = default;
#endif

    void define (
        int a_num_real_comps,
        int a_num_int_comps,
        std::vector<std::string>* a_rdata_names=nullptr,
        std::vector<std::string>* a_idata_names=nullptr,
        Arena* a_arena=nullptr
    )
    {
        if (m_defined) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_arena == m_arena,
                "ParticleTileRT redefined with different memory arena");
            if (a_num_real_comps != m_n_real || a_num_int_comps != m_n_int) {
                realloc_and_move(m_size, m_capacity, a_num_real_comps, a_num_int_comps);
            }
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_arena != nullptr,
                "ParticleTileRT defined with no memory arena! "
                "Make sure to call setArena() on the ParticleContainer before initialization or "
                "to pass an Arena to ParticleTile::define()");
            m_defined = true;
            m_arena = a_arena;
            m_n_real = a_num_real_comps;
            m_n_int = a_num_int_comps;
        }

        if (a_rdata_names != nullptr) {
            m_real_names = a_rdata_names;
        }
        if (a_idata_names != nullptr) {
            m_int_names = a_idata_names;
        }
    }

    [[nodiscard]] bool empty () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return size() == 0;
    }

    /**
    * \brief Returns the number of particles
    */
    [[nodiscard]] size_type size () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    /**
    * \brief Returns the number of particles
    */
    [[nodiscard]] size_type numParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    /**
    * \brief Returns the number of particles
    */
    [[nodiscard]] size_type numRealParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    /**
    * \brief Returns the number of neighbor particles.
    * For ParticleTileRT this is always zero.
    */
    [[nodiscard]] size_type numNeighborParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return 0;
    }

    /**
    * \brief Returns the number of particles
    */
    [[nodiscard]] size_type numTotalParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    /**
    * \brief Returns the number of neighbor particles.
    * For ParticleTileRT this is always zero.
    */
    [[nodiscard]] size_type getNumNeighbors () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return 0;
    }

    [[nodiscard]] int NumRealComps () const noexcept
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_n_real;
    }

    [[nodiscard]] int NumIntComps () const noexcept
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_n_int;
    }

    [[nodiscard]] int NumRuntimeRealComps () const noexcept
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_n_real;
    }

    [[nodiscard]] int NumRuntimeIntComps () const noexcept
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_n_int;
    }

    /** Get the names for the real SoA components **/
    [[nodiscard]] std::vector<std::string> GetRealNames () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        if (m_real_names != nullptr) {
            return *m_real_names;
        } else {
            return std::vector<std::string>();
        }
    }

    /** Get the names for the int SoA components **/
    [[nodiscard]] std::vector<std::string> GetIntNames () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        if (m_int_names != nullptr) {
            return *m_int_names;
        } else {
            return std::vector<std::string>();
        }
    }

    /** Get access to the particle id/cpu Arrays */
    [[nodiscard]] ArrayView<uint64_t> GetIdCPUData () {
        AMREX_ALWAYS_ASSERT(m_defined);
        return {m_idcpu_data.data(), m_capacity};
    }

    /** Get access to the particle id/cpu Arrays */
    [[nodiscard]] ArrayView<const uint64_t> GetIdCPUData () const {
        AMREX_ALWAYS_ASSERT(m_defined);
        return {m_idcpu_data.data(), m_capacity};
    }

    /** Get access to a particle Real component Array
     *
     * \param comp_index component with 0...NumRealComps-1
     */
    [[nodiscard]] ArrayView<RType> GetRealData (int comp_index) {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return {m_real_data.data() + comp_index * m_capacity, m_capacity};
    }

    /** Get access to a particle Real component Array
     *
     * \param comp_index component with 0...NumRealComps-1
     */
    [[nodiscard]] ArrayView<const RType> GetRealData (int comp_index) const {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return {m_real_data.data() + comp_index * m_capacity, m_capacity};
    }

    /** Get access to a particle Int component Array
     *
     * \param comp_index component with 0...NumIntComps-1
     */
    [[nodiscard]] ArrayView<IType> GetIntData (int comp_index) {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return {m_int_data.data() + comp_index * m_capacity, m_capacity};
    }

    /** Get access to a particle Int component Array
     *
     * \param comp_index component with 0...NumIntComps-1
     */
    [[nodiscard]] ArrayView<const IType> GetIntData (int comp_index) const {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return {m_int_data.data() + comp_index * m_capacity, m_capacity};
    }

    /** Get access to a particle Real component Array
     *
     * \param name name of the component
     */
    [[nodiscard]] ArrayView<RType> GetRealData (std::string const & name) {
        return GetRealData(get_idx_from_str(name, m_real_names));
    }

    /** Get access to a particle Real component Array
     *
     * \param name name of the component
     */
    [[nodiscard]] ArrayView<const RType> GetRealData (std::string const & name) const {
        return GetRealData(get_idx_from_str(name, m_real_names));
    }

    /** Get access to a particle Int component Array
     *
     * \param name name of the component
     */
    [[nodiscard]] ArrayView<IType> GetIntData (std::string const & name) {
        return GetIntData(get_idx_from_str(name, m_int_names));
    }

    /** Get access to a particle Int component Array
     *
     * \param name name of the component
     */
    [[nodiscard]] ArrayView<const IType> GetIntData (std::string const & name) const {
        return GetIntData(get_idx_from_str(name, m_int_names));
    }

    void resize (size_type new_size, GrowthStrategy strategy = GrowthStrategy::Poisson)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (m_capacity < new_size) {
            size_type new_grown_capacity =
                grow_podvector_capacity(strategy, new_size, m_capacity, sizeof(RType));
            realloc_and_move(new_size, new_grown_capacity, m_n_real, m_n_int);
        } else {
            m_size = new_size;
        }
    }

    void reserve (size_type new_capacity, GrowthStrategy strategy = GrowthStrategy::Poisson)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (m_capacity < new_capacity) {
            size_type new_grown_capacity =
                grow_podvector_capacity(strategy, new_capacity, m_capacity, sizeof(RType));
            realloc_and_move(m_size, new_grown_capacity, m_n_real, m_n_int);
        }
    }

    static void reserve (std::map<ParticleTileRT<RType, IType>*, int> const& addsizes)
    {
        // we don't do anything special here for now
        for (auto [p,s] : addsizes) {
            p->reserve(p->size()+s);
        }
    }

    void realloc_and_move (size_type new_size, size_type new_capacity,
                           int new_n_real, int new_n_int)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (new_capacity <= new_size) {
            new_capacity = new_size;
        }

        decltype(m_idcpu_data) new_idcpu_data;
        decltype(m_real_data) new_real_data;
        decltype(m_int_data) new_int_data;

        new_idcpu_data.setArena(m_arena);
        new_real_data.setArena(m_arena);
        new_int_data.setArena(m_arena);

        new_idcpu_data.resize(new_capacity, GrowthStrategy::Exact);
        new_real_data.resize(new_capacity * new_n_real, GrowthStrategy::Exact);
        new_int_data.resize(new_capacity * new_n_int, GrowthStrategy::Exact);

        auto p_idcpu = new_idcpu_data.data();
        auto p_real = new_real_data.data();
        auto p_int = new_int_data.data();

        auto old_ptd = getParticleTileData();
        auto copy_size = std::min(m_size, new_size);
        auto copy_n_real = std::min(m_n_real, new_n_real);
        auto copy_n_int = std::min(m_n_int, new_n_int);

        m_size = new_size;
        m_capacity = new_capacity;
        m_n_real = new_n_real;
        m_n_int = new_n_int;
        m_idcpu_data.swap(new_idcpu_data);
        m_real_data.swap(new_real_data);
        m_int_data.swap(new_int_data);

        auto new_ptd = getParticleTileData();

        if (new_size == 0) {
            return;
        }

        if (m_arena->isDeviceAccessible()) {
            ParallelFor(copy_size,
                [=] AMREX_GPU_DEVICE (size_type i) noexcept {
                    new_ptd.idcpu(i) = old_ptd.idcpu(i);

                    for (int n = 0;  n < copy_n_real; ++n) {
                        new_ptd.rdata(n, i) = old_ptd.rdata(n, i);
                    }

                    for (int n = 0;  n < copy_n_int; ++n) {
                        new_ptd.idata(n, i) = old_ptd.idata(n, i);
                    }
                });
            Gpu::streamSynchronize();
        } else {
            for (size_type i = 0; i < copy_size; ++i) {
                new_ptd.idcpu(i) = old_ptd.idcpu(i);

                for (int n = 0;  n < copy_n_real; ++n) {
                    new_ptd.rdata(n, i) = old_ptd.rdata(n, i);
                }

                for (int n = 0;  n < copy_n_int; ++n) {
                    new_ptd.idata(n, i) = old_ptd.idata(n, i);
                }
            }
        }
    }

    void shrink_to_fit ()
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (m_size != m_capacity) {
            realloc_and_move(m_size, m_size, m_n_real, m_n_int);
        }
    }

    size_type capacity () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_capacity * (sizeof(uint64_t) + m_n_real * sizeof(RType) + m_n_int * sizeof(IType));
    }

    void swap (ParticleTileRT<RType, IType>& other) noexcept
    {
        std::swap(m_defined, other.m_defined);
        std::swap(m_arena, other.m_arena);
        std::swap(m_capacity, other.m_capacity);
        std::swap(m_size, other.m_size);
        std::swap(m_n_real, other.m_n_real);
        std::swap(m_n_int, other.m_n_int);
        m_idcpu_data.swap(other.m_idcpu_data);
        m_real_data.swap(other.m_real_data);
        m_int_data.swap(other.m_int_data);
        std::swap(m_real_names, other.m_real_names);
        std::swap(m_int_names, other.m_int_names);
    }

    ParticleTileDataType getParticleTileData ()
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return ParticleTileDataType{
            m_capacity,
            m_idcpu_data.data(),
            m_real_data.data(),
            m_int_data.data(),
            m_n_real,
            m_n_int
        };
    }

    ConstParticleTileDataType getConstParticleTileData () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return ConstParticleTileDataType{
            m_capacity,
            m_idcpu_data.data(),
            m_real_data.data(),
            m_int_data.data(),
            m_n_real,
            m_n_int
        };
    }

    auto&       GetStructOfArrays ()       { return *this; }
    const auto& GetStructOfArrays () const { return *this; }

    [[nodiscard]] Arena* arena() const {
        return m_arena;
    }

private:

    static int get_idx_from_str (std::string const & name, std::vector<std::string>* name_list) {
        AMREX_ALWAYS_ASSERT(name_list != nullptr);
        auto const pos = std::find(name_list->begin(), name_list->end(), name);
        AMREX_ALWAYS_ASSERT(pos != name_list->end());
        return static_cast<int>(std::distance(name_list->begin(), pos));
    }

    bool m_defined = false;

    Arena* m_arena = nullptr;

    size_type m_capacity = 0;
    size_type m_size = 0;

    int m_n_real = 0;
    int m_n_int = 0;

    amrex::PODVector<uint64_t, PolymorphicArenaAllocator<uint64_t>> m_idcpu_data;
    amrex::PODVector<RType, PolymorphicArenaAllocator<RType>> m_real_data;
    amrex::PODVector<IType, PolymorphicArenaAllocator<IType>> m_int_data;

    std::vector<std::string>* m_real_names = nullptr;
    std::vector<std::string>* m_int_names = nullptr;
};

} // namespace amrex

#endif // AMREX_PARTICLETILERT_H_
