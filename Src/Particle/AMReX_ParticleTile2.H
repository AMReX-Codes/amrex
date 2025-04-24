#ifndef AMREX_PARTICLETILE2_H_
#define AMREX_PARTICLETILE2_H_
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
struct SoAParticle2;

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

struct SoASuperParticle2 {
    amrex::Real m_pos[AMREX_SPACEDIM];

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real& pos(const int dir) {
        return m_pos[dir];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real pos(const int dir) const {
        return m_pos[dir];
    }
};

template <class RType, class IType>
struct ParticleTileData2
{
    using Self = ParticleTileData2<RType, IType>;
    using size_type = std::size_t;

    static constexpr bool is_particle_tile_data = true;
    static constexpr bool is_const = std::is_const_v<RType>;

    size_type m_capacity;
    std::conditional_t<is_const, const uint64_t* AMREX_RESTRICT, uint64_t* AMREX_RESTRICT> m_idcpu;
    RType* AMREX_RESTRICT m_rdata;
    IType* AMREX_RESTRICT m_idata;
    int m_n_real;
    int m_n_int;

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
    ArrayView<RType> rdata (const int comp_index) const &
    {
        AMREX_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return {this->m_rdata + comp_index * m_capacity, m_capacity};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ArrayView<IType> idata (const int comp_index) const &
    {
        AMREX_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return {this->m_idata + comp_index * m_capacity, m_capacity};
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
        return SoAParticle2<RType, IType>(*this, index);
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SoASuperParticle2 getSuperParticle (size_type index) const noexcept
    {
        return { AMREX_D_DECL( pos(0, index), pos(1, index), pos(2, index) )};
    }
};

template <class RType, class IType>
struct SoAParticle2
{
    using PTD = ParticleTileData2<RType, IType>;
    using ConstType = SoAParticle2<std::add_const_t<RType>, std::add_const_t<IType>>;
    static constexpr bool is_constsoa_particle = PTD::is_const;
    static constexpr int NReal = 0;
    static constexpr int NInt = 0;
    static constexpr bool is_soa_particle = true;
    static constexpr bool is_soa2_particle = true;
    using size_type = typename PTD::size_type;
    using RealType = RType;
    using IntType = IType;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SoAParticle2 (PTD const& ptd, size_type i) :
        m_particle_tile_data(ptd), m_index(i) {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) cpu () const & { return this->m_particle_tile_data.cpu(m_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) id () const & { return this->m_particle_tile_data.id(m_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) idcpu () const & { return this->m_particle_tile_data.idcpu(m_index); }

    //functions to get positions of the particle in the SOA data

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect pos () const & {
        return RealVect(AMREX_D_DECL(
            this->m_particle_tile_data.pos(0, m_index),
            this->m_particle_tile_data.pos(1, m_index),
            this->m_particle_tile_data.pos(2, m_index)
        ));
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RType& pos (int position_index) const &
    {
        return this->m_particle_tile_data.pos(position_index, m_index);
    }

private :

    static_assert(std::is_trivially_copyable<PTD>(), "ParticleTileData is not trivially copyable");

    PTD m_particle_tile_data;
    size_type m_index;
};

template<class RType, class IType>
using ConstParticleTileData2 = ParticleTileData2<const RType, const IType>;

template <class RType=ParticleReal, class IType=int>
struct ParticleTile2
{
    using ParticleType = SoAParticle2<RType, IType>;
    using RealType = RType;

    using ParticleTileDataType = ParticleTileData2<RType, IType>;
    using ConstParticleTileDataType = ParticleTileData2<const RType, const IType>;

    using size_type = typename ParticleTileDataType::size_type;

    using AoS = ThisParticleTileHasNoAoS;
    using SoA = ParticleTile2<RType, IType>;

    using RealVector = ThisParticleTileHasNoRealVector;
    using IntVector = ThisParticleTileHasNoIntVector;

    ParticleTile2 () = default;

#ifndef _WIN32  // workaround windows compiler bug
    ~ParticleTile2 () = default;

    ParticleTile2 (ParticleTile2 const&) = delete;
    ParticleTile2 (ParticleTile2 &&) noexcept = default;

    ParticleTile2& operator= (ParticleTile2 const&) = delete;
    ParticleTile2& operator= (ParticleTile2 &&) noexcept = default;
#endif

    void define (
        Arena* a_arena,
        int a_num_real_comps,
        int a_num_int_comps,
        std::vector<std::string>* a_rdata_names=nullptr,
        std::vector<std::string>* a_idata_names=nullptr)
    {
        if (m_defined) {
            AMREX_ALWAYS_ASSERT(a_arena == nullptr || a_arena == m_arena);
            if (a_num_real_comps != m_n_real || a_num_int_comps != m_n_int) {
                realloc_and_move(m_size, m_capacity, a_num_real_comps, a_num_int_comps);
            }
        } else {
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
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    size_type size () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    size_type numParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return numRealParticles();
    }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    size_type numRealParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size - m_num_neighbor_particles;
    }

    /**
    * \brief Returns the number of neighbor particles (excluding reals)
    *
    */
    size_type numNeighborParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_num_neighbor_particles;
    }

    /**
    * \brief Returns the total number of particles, real and neighbor
    *
    */
    size_type numTotalParticles () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_size;
    }

    void setNumNeighbors (size_type num_neighbors)
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        auto nrp = numRealParticles();
        m_num_neighbor_particles = num_neighbors;
        resize(nrp + num_neighbors);
    }

    size_type getNumNeighbors () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        return m_num_neighbor_particles;
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

    [[nodiscard]] ArrayView<uint64_t> GetIdCPUData () {
        AMREX_ALWAYS_ASSERT(m_defined);
        return {m_idcpu_data.data(), m_capacity};
    }

    [[nodiscard]] ArrayView<const uint64_t> GetIdCPUData () const {
        AMREX_ALWAYS_ASSERT(m_defined);
        return {m_idcpu_data.data(), m_capacity};
    }

    [[nodiscard]] ArrayView<RType> GetRealData (int comp_index) {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return {m_real_data.data() + comp_index * m_capacity, m_capacity};
    }

    [[nodiscard]] ArrayView<const RType> GetRealData (int comp_index) const {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_real);
        return {m_real_data.data() + comp_index * m_capacity, m_capacity};
    }

    [[nodiscard]] ArrayView<IType> GetIntData (int comp_index) {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return {m_int_data.data() + comp_index * m_capacity, m_capacity};
    }

    [[nodiscard]] ArrayView<const IType> GetIntData (int comp_index) const {
        AMREX_ALWAYS_ASSERT(m_defined);
        AMREX_ALWAYS_ASSERT(0 <= comp_index && comp_index < m_n_int);
        return {m_int_data.data() + comp_index * m_capacity, m_capacity};
    }

    [[nodiscard]] ArrayView<RType> GetRealData (std::string const & name) {
        return GetRealData(get_idx_from_str(name, m_real_names));
    }

    [[nodiscard]] ArrayView<const RType> GetRealData (std::string const & name) const {
        return GetRealData(get_idx_from_str(name, m_real_names));
    }

    [[nodiscard]] ArrayView<IType> GetIntData (std::string const & name) {
        return GetIntData(get_idx_from_str(name, m_int_names));
    }

    [[nodiscard]] ArrayView<const IType> GetIntData (std::string const & name) const {
        return GetIntData(get_idx_from_str(name, m_int_names));
    }

    void resize (size_type new_size)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (new_size <= m_capacity && new_size * 10 > m_capacity) {
            m_size = new_size;
        } else {
            realloc_and_move(new_size, new_size, m_n_real, m_n_int);
        }
    }

    void resize_geometric_grow (size_type new_size)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (new_size <= m_capacity && new_size * 10 > m_capacity) {
            m_size = new_size;
        } else {
            realloc_and_move(new_size, std::max(new_size, m_capacity + m_capacity / 4),
                             m_n_real, m_n_int);
        }
    }

    void reserve (size_type new_capacity)
    {
        AMREX_ALWAYS_ASSERT(m_defined);

        if (new_capacity > m_capacity) {
            realloc_and_move(m_size, new_capacity, m_n_real, m_n_int);
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

        new_idcpu_data.resize(new_capacity);
        new_real_data.resize(new_capacity * new_n_real);
        new_int_data.resize(new_capacity * new_n_int);

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

    void swap (ParticleTile2<RType, IType>& other) noexcept
    {
        std::swap(m_defined, other.m_defined);
        std::swap(m_arena, other.m_arena);
        std::swap(m_capacity, other.m_capacity);
        std::swap(m_size, other.m_size);
        std::swap(m_num_neighbor_particles, other.m_num_neighbor_particles);
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
        ParticleTileDataType ptd;

        ptd.m_capacity = m_capacity;
        ptd.m_idcpu = m_idcpu_data.data();
        ptd.m_rdata = m_real_data.data();
        ptd.m_idata = m_int_data.data();
        ptd.m_n_real = m_n_real;
        ptd.m_n_int = m_n_int;

        return ptd;
    }

    ConstParticleTileDataType getConstParticleTileData () const
    {
        AMREX_ALWAYS_ASSERT(m_defined);
        ConstParticleTileDataType ptd;

        ptd.m_capacity = m_capacity;
        ptd.m_idcpu = m_idcpu_data.data();
        ptd.m_rdata = m_real_data.data();
        ptd.m_idata = m_int_data.data();
        ptd.m_n_real = m_n_real;
        ptd.m_n_int = m_n_int;

        return ptd;
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
    size_type m_num_neighbor_particles = 0;

    int m_n_real = 0;
    int m_n_int = 0;

    amrex::PODVector<uint64_t, PolymorphicArenaAllocator<uint64_t>> m_idcpu_data;
    amrex::PODVector<RType, PolymorphicArenaAllocator<RType>> m_real_data;
    amrex::PODVector<IType, PolymorphicArenaAllocator<IType>> m_int_data;

    std::vector<std::string>* m_real_names = nullptr;
    std::vector<std::string>* m_int_names = nullptr;
};

} // namespace amrex

#endif // AMREX_PARTICLETILE2_H_
