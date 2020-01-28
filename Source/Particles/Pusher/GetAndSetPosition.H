#ifndef WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
#define WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_

#include <limits>
#include <WarpXParticleContainer.H>
#include <AMReX_REAL.H>

struct GetPosition
{
    using PType = WarpXParticleContainer::ParticleType;
    using RType = amrex::ParticleReal;

    const PType* AMREX_RESTRICT m_structs;
#if (defined WARPX_DIM_RZ)
    const RType* m_theta;
#elif (AMREX_SPACEDIM == 2)
    static constexpr RType m_snan = std::numeric_limits<RType>::quiet_NaN();
#endif
    GetPosition (const WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        const auto& aos = a_pti.GetArrayOfStructs();
        m_structs = aos().dataPtr() + a_offset;
#if (defined WARPX_DIM_RZ)
        const auto& soa = a_pti.GetStructOfArrays();
        m_theta = soa.GetRealData(PIdx::theta).dataPtr() + a_offset;
#elif (AMREX_SPACEDIM == 2)
        static constexpr RType m_snan = std::numeric_limits<RType>::quiet_NaN();
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i, RType& x, RType& y, RType& z) const noexcept
    {
#ifdef WARPX_DIM_RZ
        RType r = m_structs[i].pos(0);
        x = r*std::cos(m_theta[i]);
        y = r*std::sin(m_theta[i]);
        z = m_structs[i].pos(1);
#elif WARPX_DIM_3D
        x = m_structs[i].pos(0);
        y = m_structs[i].pos(1);
        z = m_structs[i].pos(2);
#else
        x = m_structs[i].pos(0);
        y = m_snan;
        z = m_structs[i].pos(2);
#endif
    }
};

struct SetPosition
{
    using PType = WarpXParticleContainer::ParticleType;
    using RType = amrex::ParticleReal;

    PType* AMREX_RESTRICT m_structs;
#if (defined WARPX_DIM_RZ)
    RType* AMREX_RESTRICT m_theta;
#endif
    SetPosition (WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        auto& aos = a_pti.GetArrayOfStructs();
        m_structs = aos().dataPtr() + a_offset;
#if (defined WARPX_DIM_RZ)
        auto& soa = a_pti.GetStructOfArrays();
        m_theta = soa.GetRealData(PIdx::theta).dataPtr() + a_offset;
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i, RType x, RType y, RType z) const noexcept
    {
#ifdef WARPX_DIM_RZ
        m_theta[i] = std::atan2(y, x);
        m_structs[i].pos(0) = std::sqrt(x*x + y*y);
        m_structs[i].pos(1) = z;
#elif WARPX_DIM_3D
        m_structs[i].pos(0) = x;
        m_structs[i].pos(1) = y;
        m_structs[i].pos(2) = z;
#else
        m_structs[i].pos(0) = x;
        m_structs[i].pos(2) = z;
#endif
    }
};

#ifdef WARPX_DIM_RZ

/** \brief Extract the particle's coordinates from `theta` and the attributes
 *         of the ParticleType struct `p` (which contains the radius),
 *         and store them in the variables `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void GetCartesianPositionFromCylindrical(
    amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
    const WarpXParticleContainer::ParticleType& p, const amrex::ParticleReal theta)
{
    const amrex::ParticleReal r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
}

/** \brief Set the particle's cylindrical coordinates by setting `theta`
 *        and the attributes of the ParticleType struct `p` (which stores the radius),
 *        from the values of `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SetCylindricalPositionFromCartesian(
    WarpXParticleContainer::ParticleType& p, amrex::ParticleReal& theta,
    const amrex::ParticleReal x, const amrex::ParticleReal y, const amrex::ParticleReal z )
{
    theta = std::atan2(y, x);
    p.pos(0) = std::sqrt(x*x + y*y);
    p.pos(1) = z;
}

#endif // WARPX_DIM_RZ

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
