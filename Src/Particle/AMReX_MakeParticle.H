#ifndef AMREX_MAKEPARTICLE_H_
#define AMREX_MAKEPARTICLE_H_

#include <type_traits>

template< class T >
struct is_soa_particle
     : std::integral_constant<
         bool,
         T::is_soa_particle
     > {};


template <typename T_ParticleType, class Enable = void>
struct make_particle
{
    template <typename PTD>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto&
    operator() (PTD const& ptd, int i)
    {
        // legacy Particle (AoS)
        return ptd.m_aos[i];
    }
};

template <typename T_ParticleType>
struct make_particle<T_ParticleType, std::enable_if_t<is_soa_particle<T_ParticleType>::value>>
{
    template <typename PTD>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto
    operator() (PTD const& ptd, int index)
    {
        // SoAParticle
        return T_ParticleType(ptd, index);
    }
};

#endif
