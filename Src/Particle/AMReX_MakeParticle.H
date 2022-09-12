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
    auto
    operator()(PTD, int)
    {
        // legacy Particle (AoS)
        return T_ParticleType();
    }
};

template <typename T_ParticleType>
struct make_particle<T_ParticleType, typename std::enable_if<is_soa_particle<T_ParticleType>::value>::type>
{
    template <typename PTD>
    auto
    operator()(PTD ptd, int index)
    {
        // SoAParticle
        return T_ParticleType(ptd, index);
    }
};