#ifndef AMREX_PARTICLE_H_
#define AMREX_PARTICLE_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>

#include <string>

namespace amrex {

namespace
{
    constexpr Long GhostParticleID    = 549755813887L; // 2**39-1
    constexpr Long VirtualParticleID  = GhostParticleID-1;
    constexpr Long LastParticleID     = GhostParticleID-2;
    constexpr Long DoSplitParticleID  = GhostParticleID-3;
    constexpr Long NoSplitParticleID  = GhostParticleID-4;
}

struct ParticleIDWrapper
{
    uint64_t& m_idata;

    AMREX_GPU_HOST_DEVICE
    ParticleIDWrapper (uint64_t& idata) noexcept
        : m_idata(idata)
    {}

    AMREX_GPU_HOST_DEVICE
    ParticleIDWrapper (const ParticleIDWrapper& rhs) = default;

    AMREX_GPU_HOST_DEVICE
    ParticleIDWrapper& operator= (const ParticleIDWrapper& pidw) noexcept
    {
        return this->operator=(Long(pidw));
    }

    AMREX_GPU_HOST_DEVICE
    ParticleIDWrapper& operator= (const Long id) noexcept
    {
        // zero out the 40 leftmost bits, which store the sign and the abs of the id;
        m_idata &= 0x00FFFFFF;

        uint64_t val;
        uint64_t sign = id >= 0;
        if (sign)
        {
            // 2**39-1, the max value representible in this fashion
            AMREX_ASSERT(id <= 549755813887L);
            val = id;
        }
        else
        {
            // -2**39-1, the min value representible in this fashion
            AMREX_ASSERT(id >= -549755813887L);
            val = -id;
        }

        m_idata |= (sign << 63);  // put the sign in the leftmost bit
        m_idata |= (val << 24);   // put the val in the next 39
        return *this;
    }

    AMREX_GPU_HOST_DEVICE
    operator Long () const noexcept
    {
        Long r = 0;

        uint64_t sign = m_idata >> 63;  // extract leftmost sign bit
        uint64_t val  = ((m_idata >> 24) & 0x7FFFFFFFFF);  // extract next 39 id bits

        Long lval = static_cast<Long>(val);  // bc we take -
        r = (sign) ? lval : -lval;
        return r;
    }
};

struct ParticleCPUWrapper
{
    uint64_t& m_idata;

    AMREX_GPU_HOST_DEVICE
    ParticleCPUWrapper (uint64_t& idata) noexcept
        : m_idata(idata)
    {}

    AMREX_GPU_HOST_DEVICE
    ParticleCPUWrapper (const ParticleCPUWrapper& rhs) = default;

    AMREX_GPU_HOST_DEVICE
    ParticleCPUWrapper& operator= (const ParticleCPUWrapper& pcpuw) noexcept
    {
        return this->operator=(int(pcpuw));
    }

    AMREX_GPU_HOST_DEVICE
    ParticleCPUWrapper& operator= (const int cpu) noexcept
    {
        // zero out the first 24 bits, which are used to store the cpu number
        m_idata &= (~ 0x00FFFFFF);

        AMREX_ASSERT(cpu >= 0);
        AMREX_ASSERT(cpu <= 16777215);  // 2**24-1, the max representable number

        m_idata |= cpu;
        return *this;
    }

    AMREX_GPU_HOST_DEVICE
    operator int () const noexcept
    {
        return (m_idata & 0x00FFFFFF);
    }
};

struct ConstParticleIDWrapper
{
    const uint64_t& m_idata;

    AMREX_GPU_HOST_DEVICE
    ConstParticleIDWrapper (const uint64_t& idata) noexcept
        : m_idata(idata)
    {}

    AMREX_GPU_HOST_DEVICE
    operator Long () const noexcept
    {
        Long r = 0;

        uint64_t sign = m_idata >> 63;  // extract leftmost sign bit
        uint64_t val  = ((m_idata >> 24) & 0x7FFFFFFFFF);  // extract next 39 id bits

        Long lval = static_cast<Long>(val);  // bc we take -
        r = (sign) ? lval : -lval;
        return r;
    }
};

struct ConstParticleCPUWrapper
{
    const uint64_t& m_idata;

    AMREX_GPU_HOST_DEVICE
    ConstParticleCPUWrapper (const uint64_t& idata) noexcept
        : m_idata(idata)
    {}

    AMREX_GPU_HOST_DEVICE
    operator int () const noexcept { return (m_idata & 0x00FFFFFF); }
};


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


/** \brief The struct used to store particles.
 *
 * \tparam T_NReal The number of extra Real components
 * \tparam T_NInt The number of extra integer components
 */
template <int T_NReal, int T_NInt=0>
struct Particle
    : ParticleBase<ParticleReal,T_NReal,T_NInt>
{
    //! \brief number of extra Real components in the particle struct
    static constexpr int NReal = T_NReal;

    //! \brief number of extra integer components in the particle struct
    static constexpr int NInt = T_NInt;

    //! The floating point type used for the particles.
    using RealType = ParticleReal;

    static Long the_next_id;

    AMREX_GPU_HOST_DEVICE ParticleCPUWrapper cpu () & { return ParticleCPUWrapper(this->m_idcpu); }
    AMREX_GPU_HOST_DEVICE ParticleIDWrapper id () & { return ParticleIDWrapper(this->m_idcpu); }
    AMREX_GPU_HOST_DEVICE ConstParticleCPUWrapper cpu () const & { return ConstParticleCPUWrapper(this->m_idcpu); }
    AMREX_GPU_HOST_DEVICE ConstParticleIDWrapper id () const & { return ConstParticleIDWrapper(this->m_idcpu); }

    AMREX_GPU_HOST_DEVICE RealVect pos () const & {return RealVect(AMREX_D_DECL(this->m_pos[0], this->m_pos[1], this->m_pos[2]));}

    AMREX_GPU_HOST_DEVICE RealType& pos (int index) &
    {
        AMREX_ASSERT(index < AMREX_SPACEDIM);
        return this->m_pos[index];
    }

    AMREX_GPU_HOST_DEVICE RealType  pos (int index) const &
    {
        AMREX_ASSERT(index < AMREX_SPACEDIM);
        return this->m_pos[index];
    }

    template <int U = T_NReal, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealType& rdata (int index) &
    {
        AMREX_ASSERT(index < NReal);
        return this->m_rdata[index];
    }

    template <int U = T_NReal, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealType& rdata (int /*index*/) &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->pos(0);  // bc we must return something
    }

    template <int U = T_NReal, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE const RealType& rdata (int index) const &
    {
        AMREX_ASSERT(index < NReal);
        return this->m_rdata[index];
    }

    template <int U = T_NReal, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealType rdata (int /*index*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->pos(0);  // because we must return something
    }

    AMREX_GPU_HOST_DEVICE RealType rdata (int /*index*/) && = delete;

    template <int U = T_NReal, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealVect  rvec (AMREX_D_DECL(int indx, int indy, int indz)) const &
    {
        AMREX_ASSERT(AMREX_D_TERM(indx < NReal, && indy < NReal, && indz < NReal));
        return RealVect(AMREX_D_DECL(this->m_rdata[indx],
                                     this->m_rdata[indy],
                                     this->m_rdata[indz]));
    }

    template <int U = T_NReal, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealVect  rvec (AMREX_D_DECL(int /*indx*/, int /*indy*/, int /*indz*/)) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return RealVect(AMREX_D_DECL(0.0, 0.0, 0.0)); // bc we must return something
    }

    template <int U = T_NReal, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealVect  rvec (const IntVect& indices) const &
    {
        AMREX_ASSERT(indices.max() < NReal);
        return RealVect(AMREX_D_DECL(this->m_rdata[indices[0]],
                                     this->m_rdata[indices[1]],
                                     this->m_rdata[indices[2]]));
    }

    template <int U = T_NReal, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE RealVect  rvec (const IntVect& /*indices*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return RealVect(AMREX_D_DECL(0.0, 0.0, 0.0)); // bc we must return something
    }

    template <int U = T_NInt, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE int& idata (int index) &
    {
        AMREX_ASSERT(index < NInt);
        return this->m_idata[index];
    }

    template <int U = T_NInt, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE uint64_t& idata (int /*index*/) &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->m_idcpu;  //bc we must return something
    }

    template <int U = T_NInt, typename std::enable_if<U != 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE const int& idata (int index) const &
    {
        AMREX_ASSERT(index < NInt);
        return this->m_idata[index];
    }

    template <int U = T_NInt, typename std::enable_if<U == 0, int>::type = 0>
    AMREX_GPU_HOST_DEVICE int idata (int /*index*/) const &
    {
        AMREX_ALWAYS_ASSERT(false);
        return this->m_idcpu;  //bc we must return something
    }

    AMREX_GPU_HOST_DEVICE RealType idata (int /*index*/) && = delete;

    /**
    * \brief Returns the next particle ID for this processor.
    * Particle IDs start at 1 and are never reused.
    * The pair, consisting of the ID and the CPU on which the particle is "born",
    * is a globally unique identifier for a particle.  The maximum of this value
    * across all processors must be checkpointed and then restored on restart
    * so that we don't reuse particle IDs.
    */
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

    if (next > LastParticleID)
        amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");

    return next;
}

template <int NReal, int NInt>
Long
Particle<NReal, NInt>::UnprotectedNextID ()
{
    Long next = the_next_id++;
    if (next > LastParticleID)
        amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");
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

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << p.pos(i) << ' ';

    for (int i = 0; i < NReal; i++)
        os << p.rdata(i) << ' ';

    for (int i = 2; i < NInt; i++)
        os << p.idata(i) << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");

    return os;
}

template <int NReal>
std::ostream&
operator<< (std::ostream& os, const Particle<NReal, 0>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << p.pos(i) << ' ';

    for (int i = 0; i < NReal; i++)
        os << p.rdata(i) << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");

    return os;
}

template <int NInt>
std::ostream&
operator<< (std::ostream& os, const Particle<0, NInt>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << p.pos(i) << ' ';

    for (int i = 2; i < NInt; i++)
        os << p.idata(i) << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");

    return os;
}

template <int NReal=0, int NInt=0>
std::ostream&
operator<< (std::ostream& os, const Particle<0, 0>& p)
{
    os << p.id()   << ' '
       << p.cpu()  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << p.pos(i) << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");

    return os;
}

} // namespace amrex

#endif // AMREX_PARTICLE_H_
