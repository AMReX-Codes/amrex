#ifndef DEFAULTINITIALIZATION_H_
#define DEFAULTINITIALIZATION_H_

#include "AMReX_REAL.H"
#include "AMReX_GpuContainers.H"

#include <map>
#include <string>

/**
 * \brief This set of initialization policies describes what happens
 * when we need to create a new particle due to an elementary process.
 * For example, when an ionization event creates an electron, these
 * policies control the initial values of the electron's components.
 * These can always be over-written later.
 *
 * The specific meanings are as follows:
 *     Zero         - set the component to zero
 *     One          - set the component to one
 *     SignalingNan - set the component to a signalling nan. If you run
 *                    with amrex.fpe_trap_invalid=1 AND built the code
 *                    with gcc, this will throw an exception the first time
 *                    it is used.
 *     RandomTau    - a special flag for the "tau" component used by certain
 *                    QED processes, which gets a random initial value.
 *
 */
enum struct InitializationPolicy {Zero=0, One, SignalingNan, RandomTau};

/**
 * \brief This is a functor that implements the above InitializationPolicies
 * in a switch statement. It is designed to be used in device code when
 * running on GPUs.
 */
struct DefaultInitializer
{
    const InitializationPolicy* m_policies;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int j) const noexcept
    {
        switch (m_policies[j]) {
            case InitializationPolicy::Zero : return 0.0;
            case InitializationPolicy::One  : return 1.0;
            case InitializationPolicy::SignalingNan : {
                return 1.0;
            }
            case InitializationPolicy::RandomTau : {
                return 1.0;
            }
            default : {
                amrex::Abort("Initialization Policy not recognized");
                return 1.0;
            }
        }
    }
};

/**
 * \brief This map sets the initialization policy for each particle component
 * used in WarpX.
 */
static std::map<std::string, InitializationPolicy> initialization_policies = {
    {"w",     InitializationPolicy::Zero },
    {"ux",    InitializationPolicy::Zero },
    {"uy",    InitializationPolicy::Zero },
    {"uz",    InitializationPolicy::Zero },
    {"Ex",    InitializationPolicy::Zero },
    {"Ey",    InitializationPolicy::Zero },
    {"Ez",    InitializationPolicy::Zero },
    {"Bx",    InitializationPolicy::Zero },
    {"By",    InitializationPolicy::Zero },
    {"Bz",    InitializationPolicy::Zero },
    #ifdef WARPX_DIM_RZ
    {"theta", InitializationPolicy::Zero},
    #endif
    {"tau", InitializationPolicy::RandomTau}
};

/**
 * \brief This is a factory class that creates instances of "DefaultInitializer".
 * Each particle type in WarpX has one of these associated with it. While
 * "DefaultInitializer" is a lightweight class intended to be used on the device,
 * this class is meant to be used on the host.
 */
class DefaultInitializerFactory
{
    amrex::Gpu::DeviceVector<InitializationPolicy> m_policies;
    bool m_defined;
    
public:

    DefaultInitializerFactory () noexcept : m_defined(false) {};
    
    void define (const std::map<std::string, int>& comps) noexcept
    {
        m_policies.resize(comps.size());
        for (const auto& kv : comps)
        {
            m_policies[kv.second] = initialization_policies[kv.first];
            m_defined = true;
        }
    }

    DefaultInitializer getDefaultInitializer () const noexcept
    {
        AMREX_ASSERT(m_defined);
        return DefaultInitializer{m_policies.dataPtr()};
    }

    bool isDefined () const noexcept { return m_defined; }
};

#endif
