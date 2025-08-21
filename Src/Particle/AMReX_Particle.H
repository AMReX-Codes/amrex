#ifndef AMREX_PARTICLE_H_
#define AMREX_PARTICLE_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

#include <type_traits>


namespace amrex {

    /** Special flags used for 64-bit Long particle Ids */
    namespace LongParticleIds {
        constexpr Long GhostParticleID = 549755813887L; // 2**39-1
        constexpr Long VirtualParticleID = GhostParticleID - 1;
        constexpr Long LastParticleID = GhostParticleID - 2;
        constexpr Long DoSplitParticleID = GhostParticleID - 3;
        constexpr Long NoSplitParticleID = GhostParticleID - 4;
    }

    using namespace LongParticleIds;

    /** Flags used to set the entire uint64_t idcpu
        to special values at once.
     */
    namespace ParticleIdCpus {
        constexpr std::uint64_t Invalid = 16777216; // corresponds to id = -1, cpu = 0
    }

    using namespace ParticleIdCpus;

    namespace particle_impl {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Long unpack_id (const uint64_t idcpu) noexcept {
        Long r = 0;

        uint64_t sign = idcpu >> 63;  // extract leftmost sign bit
        uint64_t val  = ((idcpu >> 24) & 0x7FFFFFFFFF);  // extract next 39 id bits

        Long lval = static_cast<Long>(val);  // bc we take -
        r = (sign) ? lval : -lval;
        return r;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int unpack_cpu (const uint64_t idcpu) noexcept {
        return static_cast<int>(idcpu & 0x00FFFFFF);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void pack_id (uint64_t& idcpu, const Long id) noexcept {
        // zero out the 40 leftmost bits, which store the sign and the abs of the id;
        idcpu &= 0x00FFFFFF;

        uint64_t val;
        uint64_t sign = id >= 0;
        if (sign)
        {
            // 2**39-1, the max value representable in this fashion
            AMREX_ASSERT(id <= 549755813887L);
            val = id;
        }
        else
        {
            // -2**39-1, the min value representable in this fashion
            AMREX_ASSERT(id >= -549755813887L);
            val = -id;
        }

        (idcpu) |= (sign << 63);  // put the sign in the leftmost bit
        (idcpu) |= (val << 24);   // put the val in the next 39
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void pack_cpu (uint64_t& idcpu, const int cpu) noexcept {
        // zero out the first 24 bits, which are used to store the cpu number
        idcpu &= (~ 0x00FFFFFF);

        AMREX_ASSERT(cpu >= 0);
        AMREX_ASSERT(cpu <= 16777215);  // 2**24-1, the max representable number

        idcpu |= cpu;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_invalid (uint64_t& idcpu) noexcept {
        // RHS mask: 0111...
        idcpu &= ~(uint64_t(1) << 63);
    }

    template<typename T_SIMD, typename T_Mask>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_invalid (T_SIMD & idcpu, T_Mask const & mask) noexcept {
        // RHS mask: 0111...
        if constexpr (std::is_same_v<T_Mask, bool>) {
            // scalar
            static_assert(std::is_same_v<T_SIMD, uint64_t>, "make_invalid: forgot template on wrapper?");
            if (mask) { make_invalid(idcpu); }
        } else {
#ifdef AMREX_USE_SIMD
            // SIMD
            static_assert(!std::is_same_v<T_SIMD, uint64_t>, "make_invalid: forgot template on wrapper?");
            // cvt / .__cvt() because of: https://github.com/VcDevel/std-simd/issues/41#issuecomment-3185164952
            // see also https://github.com/mattkretz/vir-simd/issues/45
            simd::stdx::where(simd::stdx::cvt(mask), idcpu) &= ~(T_SIMD{1} << 63);
#else
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "make_invalid: AMReX_SIMD was not enabled!");
#endif
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_valid (uint64_t& idcpu) noexcept {
        // RHS mask: 1000...
        idcpu |= uint64_t(1) << 63;
    }

    template<typename T_SIMD, typename T_Mask>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_valid (T_SIMD & idcpu, T_Mask const & mask) noexcept {
        // RHS mask: 1000...
        if constexpr (std::is_same_v<T_Mask, bool>) {
            // scalar
            static_assert(std::is_same_v<T_SIMD, uint64_t>, "make_valid: forgot template on wrapper?");
            if (mask) { make_valid(idcpu); }
        } else {
#ifdef AMREX_USE_SIMD
            // SIMD
            static_assert(!std::is_same_v<T_SIMD, uint64_t>, "make_valid: forgot template on wrapper?");
            // cvt / .__cvt() because of: https://github.com/VcDevel/std-simd/issues/41#issuecomment-3185164952
            // see also https://github.com/mattkretz/vir-simd/issues/45
            simd::stdx::where(simd::stdx::cvt(mask), idcpu) |= T_SIMD{1} << 63;
#else
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "make_valid: AMReX_SIMD was not enabled!");
#endif
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_valid (const uint64_t idcpu) noexcept {
        // the leftmost bit is our id's valid sign
        return idcpu >> 63;
    }
    } // namespace particle_impl

template<typename T = uint64_t>
struct ParticleIDWrapper
{
    /*
    static_assert(std::is_same<T, uint64_t>::value ||
                  std::is_same<T, simd::SIMDIdCpu<>>::value,  // note: this one will not work because users can change the SIMDIdCpu width
                  "ParticleIDWrapper only supports uint64_t or SIMD thereof"
    );
    */

    T* m_idata;

    ~ParticleIDWrapper () noexcept = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper (T& idata) noexcept
        : m_idata(&idata)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper (const ParticleIDWrapper& rhs) = default;

    ParticleIDWrapper (ParticleIDWrapper&&) = delete;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper& operator= (const ParticleIDWrapper& pidw) noexcept
    {
        return this->operator=(Long(pidw)); // NOLINT
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper& operator= (ParticleIDWrapper&& pidw) noexcept
    {
        return this->operator=(Long(pidw)); // NOLINT
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper& operator= (const Long id) noexcept
    {
        particle_impl::pack_id(*m_idata, id);
        return *this;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    operator Long () const noexcept
    {
        return particle_impl::unpack_id(*m_idata);
    }

    /** Mark the particle as invalid
     *
     * Swaps the is_valid (sign) bit to invalid.
     * This is NOT identical to id = -id, but it is equally reversible via make_valid().
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_invalid () const noexcept
    {
        particle_impl::make_invalid(*m_idata);
    }

    /** Mark the particle as invalid
     *
     * Swaps the is_valid (sign) bit to invalid.
     * This is NOT identical to id = -id, but it is equally reversible via make_valid().
     */
    template<typename T_Bool>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_invalid (T_Bool const & mask) const noexcept
    {
        particle_impl::make_invalid(*m_idata, mask);
    }

    /** Mark the particle as valid
     *
     * Swaps the is_valid (sign) bit to valid.
     * This is NOT identical to id = -id, but it is equally reversible via make_invalid().
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_valid () const noexcept
    {
        particle_impl::make_valid(*m_idata);
    }

    /** Mark the particle as valid
     *
     * Swaps the is_valid (sign) bit to valid.
     * This is NOT identical to id = -id, but it is equally reversible via make_invalid().
     */
    template<typename T_Bool>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void make_valid (T_Bool const & mask) const noexcept
    {
        particle_impl::make_valid(*m_idata, mask);
    }

    /** Check the particle is valid, via the sign of the id.
     *
     * Returns true if the particle is valid (the id is positive), otherwise false (invalid particle).
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_valid () const noexcept
    {
        return particle_impl::is_valid(*m_idata);
    }
};

struct ParticleCPUWrapper
{
    uint64_t* m_idata;

    ~ParticleCPUWrapper () noexcept = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper (uint64_t& idata) noexcept
        : m_idata(&idata)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper (const ParticleCPUWrapper& rhs) = default;

    ParticleCPUWrapper (ParticleCPUWrapper&&) = delete;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper& operator= (const ParticleCPUWrapper& pcpuw) noexcept
    {
        return this->operator=(int(pcpuw));  // NOLINT
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper& operator= (ParticleCPUWrapper&& pcpuw) noexcept
    {
        return this->operator=(int(pcpuw)); // NOLINT
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper& operator= (const int cpu) noexcept
    {
        particle_impl::pack_cpu(*m_idata, cpu);
        return *this;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    operator int () const noexcept
    {
        return particle_impl::unpack_cpu(*m_idata);
    }
};

struct ConstParticleIDWrapper
{
    const uint64_t* m_idata;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ConstParticleIDWrapper (const uint64_t& idata) noexcept
        : m_idata(&idata)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    operator Long () const noexcept
    {
        return particle_impl::unpack_id(*m_idata);
    }

    /** Check the sign of the id.
     *
     * Returns true if the id is positive, otherwise false (invalid particle).
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_valid () const noexcept
    {
        return particle_impl::is_valid(*m_idata);
    }
};

struct ConstParticleCPUWrapper
{
    const uint64_t* m_idata;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ConstParticleCPUWrapper (const uint64_t& idata) noexcept
        : m_idata(&idata)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    operator int () const noexcept {
        return particle_impl::unpack_cpu(*m_idata);
    }
};

/** Set the idcpu value at once, based on a particle id and cpuid
 *
 * This can be used in initialization and assignments,
 * to avoid writing twice into the same memory bank.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t SetParticleIDandCPU (Long id, int cpu) noexcept{
    std::uint64_t idcpu = 0;
    ParticleIDWrapper{idcpu} = id;
    ParticleCPUWrapper{idcpu} = cpu;
    return idcpu;
}

template <typename T, int NReal, int NInt>
struct ParticleBase
{
    T m_pos[AMREX_SPACEDIM];
    T m_rdata[NReal];
    uint64_t m_idcpu = 0;
    int m_idata[NInt];
};

template <typename T, int NInt>
struct ParticleBase<T,0,NInt>
{
    T m_pos[AMREX_SPACEDIM];
    uint64_t m_idcpu = 0;
    int m_idata[NInt];
};

template <typename T, int NReal>
struct ParticleBase<T,NReal,0>
{
    T m_pos[AMREX_SPACEDIM];
    T m_rdata[NReal];
    uint64_t m_idcpu = 0;
};

template <typename T>
struct ParticleBase<T,0,0>
{
    T m_pos[AMREX_SPACEDIM];
    uint64_t m_idcpu = 0;
};


struct SoAParticleBase
{
    static constexpr int NReal=0;
    static constexpr int NInt=0;
    static constexpr bool is_soa_particle = true;
};

/** \brief The struct used to store particles.
 *
 * \tparam T_NReal The number of extra Real components
 * \tparam T_NInt The number of extra integer components
 */
template <int T_NReal, int T_NInt=0>
struct alignas(sizeof(double)) Particle
    : ParticleBase<ParticleReal,T_NReal,T_NInt>
{
    static constexpr bool is_soa_particle = false;
    using StorageParticleType = Particle;
    using ConstType = Particle const;

    //! \brief number of extra Real components in the particle struct
    static constexpr int NReal = T_NReal;

    //! \brief number of extra integer components in the particle struct
    static constexpr int NInt = T_NInt;

    //! The floating point type used for the particles.
    using RealType = ParticleReal;
    using IntType = int;

    static Long the_next_id;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleCPUWrapper cpu () & { return ParticleCPUWrapper(this->m_idcpu); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleIDWrapper<> id () & { return ParticleIDWrapper(this->m_idcpu); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ConstParticleCPUWrapper cpu () const & { return ConstParticleCPUWrapper(this->m_idcpu); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ConstParticleIDWrapper id () const & { return ConstParticleIDWrapper(this->m_idcpu); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void atomicSetID (const Long id) {
        uint64_t tmp = 0;
        ParticleIDWrapper wrapper(tmp);
        wrapper = id;
#if defined(AMREX_USE_OMP)
#pragma omp atomic write
        this->m_idcpu = *(wrapper.m_idata);
#else
        auto *old_ptr = reinterpret_cast<unsigned long long*>(&(this->m_idcpu));
        amrex::Gpu::Atomic::Exch(old_ptr, (unsigned long long) (*wrapper.m_idata));
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t& idcpu () & { return this->m_idcpu; }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const uint64_t& idcpu () const & { return this->m_idcpu; }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect pos () const & {return RealVect(AMREX_D_DECL(this->m_pos[0], this->m_pos[1], this->m_pos[2]));}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType& pos (int index) &
    {
        AMREX_ASSERT(index < AMREX_SPACEDIM);
        return this->m_pos[index];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType  pos (int index) const &
    {
        AMREX_ASSERT(index < AMREX_SPACEDIM);
        return this->m_pos[index];
    }

    template <int U = T_NReal, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType& rdata (int index) &
    {
        AMREX_ASSERT(index < NReal);
        return this->m_rdata[index];
    }

    template <int U = T_NReal, std::enable_if_t<U == 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType& rdata (int /*index*/) &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->pos(0);  // bc we must return something
    }

    template <int U = T_NReal, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const RealType& rdata (int index) const &
    {
        AMREX_ASSERT(index < NReal);
        return this->m_rdata[index];
    }

    template <int U = T_NReal, std::enable_if_t<U == 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType rdata (int /*index*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->pos(0);  // because we must return something
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType rdata (int /*index*/) && = delete;

    template <int U = T_NReal, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect rvec (AMREX_D_DECL(int indx, int indy, int indz)) const &
    {
        AMREX_ASSERT(AMREX_D_TERM(indx < NReal, && indy < NReal, && indz < NReal));
        return RealVect(AMREX_D_DECL(this->m_rdata[indx],
                                     this->m_rdata[indy],
                                     this->m_rdata[indz]));
    }

    template <int U = T_NReal, std::enable_if_t<U == 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect rvec (AMREX_D_DECL(int /*indx*/, int /*indy*/, int /*indz*/)) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return RealVect(AMREX_D_DECL(0.0, 0.0, 0.0)); // bc we must return something
    }

    template <int U = T_NReal, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect rvec (const IntVect& indices) const &
    {
        AMREX_ASSERT(indices.allLT(NReal));
        return RealVect(AMREX_D_DECL(this->m_rdata[indices[0]],
                                     this->m_rdata[indices[1]],
                                     this->m_rdata[indices[2]]));
    }

    template <int U = T_NReal, std::enable_if_t<U == 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVect rvec (const IntVect& /*indices*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return RealVect(AMREX_D_DECL(0.0, 0.0, 0.0)); // bc we must return something
    }

    template <int U = T_NInt, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int& idata (int index) &
    {
        AMREX_ASSERT(index < NInt);
        return this->m_idata[index];
    }

    template <int U = T_NInt, std::enable_if_t<U == 0, int> = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t& idata (int /*index*/) &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->m_idcpu;  //bc we must return something
    }

    template <int U = T_NInt, std::enable_if_t<U != 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int& idata (int index) const &
    {
        AMREX_ASSERT(index < NInt);
        return this->m_idata[index];
    }

    template <int U = T_NInt, std::enable_if_t<U == 0, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int idata (int /*index*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->m_idcpu;  //bc we must return something
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealType idata (int /*index*/) && = delete;

    /**
    * \brief Returns the next particle ID for this processor.
    * Particle IDs start at 1 and are never reused.
    * The pair, consisting of the ID and the CPU on which the particle is "born",
    * is a globally unique identifier for a particle.  The maximum of this value
    * across all processors must be checkpointed and then restored on restart
    * so that we don't reuse particle IDs.
    */
    [[nodiscard]] static Long NextID ();

    /**
    * \brief This version can only be used inside omp critical.
    */
    [[nodiscard]] static Long UnprotectedNextID ();

    /**
    * \brief Reset on restart.
    *
    * \param nextid
    */
    static void NextID (Long nextid);
};

template <int NReal, int NInt> Long Particle<NReal, NInt>::the_next_id = 1;

template <int NReal, int NInt>
Long
Particle<NReal, NInt>::NextID ()
{
    Long next;
// we should be able to test on _OPENMP < 201107 for capture (version 3.1)
// but we must work around a bug in gcc < 4.9
#if defined(AMREX_USE_OMP) && defined(_OPENMP) && _OPENMP < 201307
#pragma omp critical (amrex_particle_nextid)
#elif defined(AMREX_USE_OMP)
#pragma omp atomic capture
#endif
    next = the_next_id++;

    if (next > LongParticleIds::LastParticleID) {
        amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");
    }

    return next;
}

template <int NReal, int NInt>
Long
Particle<NReal, NInt>::UnprotectedNextID ()
{
    Long next = the_next_id++;
    if (next > LongParticleIds::LastParticleID) {
        amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");
    }
    return next;
}

template <int NReal, int NInt>
void
Particle<NReal, NInt>::NextID (Long nextid)
{
    the_next_id = nextid;
}

template <int NReal, int NInt>
std::ostream&
operator<< (std::ostream& os, const Particle<NReal, NInt>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        os << p.pos(i) << ' ';
    }

    for (int i = 0; i < NReal; i++) {
        os << p.rdata(i) << ' ';
    }

    for (int i = 0; i < NInt; i++) {
        os << p.idata(i) << ' ';
    }

    if (!os.good()) {
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
    }

    return os;
}

template <int NReal>
std::ostream&
operator<< (std::ostream& os, const Particle<NReal, 0>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        os << p.pos(i) << ' ';
    }

    for (int i = 0; i < NReal; i++) {
        os << p.rdata(i) << ' ';
    }

    if (!os.good()) {
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
    }

    return os;
}

template <int NInt>
std::ostream&
operator<< (std::ostream& os, const Particle<0, NInt>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        os << p.pos(i) << ' ';
    }

    for (int i = 0; i < NInt; i++) {
        os << p.idata(i) << ' ';
    }

    if (!os.good()) {
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
    }

    return os;
}

template <int NReal=0, int NInt=0>
std::ostream&
operator<< (std::ostream& os, const Particle<0, 0>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        os << p.pos(i) << ' ';
    }

    if (!os.good()) {
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
    }

    return os;
}

} // namespace amrex

#endif // AMREX_PARTICLE_H_
