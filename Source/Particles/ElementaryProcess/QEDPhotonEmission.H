/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef QED_PHOTON_EMISSION_H_
#define QED_PHOTON_EMISSION_H_

#include "Utils/WarpXConst.H"
#include "Particles/WarpXParticleContainer.H"

#include <QuantumSyncEngineWrapper.H>

/** @file
 *
 * This file contains the implementation of the elementary process
 * functors needed for QED photon emission (an electron or a positron
 * emits a photon).
 */

/**
 * \brief Filter functor for the QED photon emission process
 */
class PhotonEmissionFilterFunc
{
public:

    /**
    * \brief Constructor of the PhotonEmissionFilterFunc functor.
    *
    * @param[in] opt_depth_runtime_comp Index of the optical depth component
    */
    PhotonEmissionFilterFunc(int const opt_depth_runtime_comp):
        m_opt_depth_runtime_comp{opt_depth_runtime_comp}
    {}

    /**
    * \brief Functor call. This method determines if a given (electron or positron)
    * particle should undergo QED photon emission.
    *
    * @param[in] ptd particle tile data
    * @param[in] i particle index
    * @return true if a pair has to be generated, false otherwise
    */
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int const i) const noexcept
    {
        using namespace amrex;

        const amrex::ParticleReal opt_depth =
            ptd.m_runtime_rdata[m_opt_depth_runtime_comp][i];
        return (opt_depth < 0.0_rt);
    }

private:
    int m_opt_depth_runtime_comp /*!< Index of the optical depth component of the source species*/;
};

/**
 * \brief Transform functor for the QED photon emission process
 */
class PhotonEmissionTransformFunc
{

public:

    /**
    * \brief Constructor of the PhotonEmissionTransformFunc functor.
    *
    * A QuantumSynchrotronGeneratePhotonAndUpdateMomentum functor is passed by value. However, it contains
    * only few integer and real parameters and few pointers to the raw data of the
    * lookup tables. Therefore, it should be rather lightweight to copy.
    *
    * Also a QuantumSynchrotronGetOpticalDepth has to be passed, since the
    * optical depth has to be re-initialized after each photon emission.
    *
    * @param[in] opt_depth_functor functor to re-initialize the optical depth of the source particles
    * @param[in] opt_depth_runtime_comp index of the optical depth component of the source species
    * @param[in] emission_functor functor to generate photons and update momentum of the source particles
    */
    PhotonEmissionTransformFunc(
        QuantumSynchrotronGetOpticalDepth opt_depth_functor,
        int const opt_depth_runtime_comp,
        QuantumSynchrotronGeneratePhotonAndUpdateMomentum const emission_functor
        ):
        m_opt_depth_functor{opt_depth_functor},
        m_opt_depth_runtime_comp{opt_depth_runtime_comp},
        m_emission_functor{emission_functor}
        {}

    /**
    * \brief Functor call. It determines the properties of the generated photon
    * and updates the momentum of the source particle
    *
    * @param[in,out] dst target species (photons)
    * @param[in, out] src source species (either electrons or positrons)
    * @param[in] i_src particle index of the source species
    * @param[in] i_dst particle index of target species
    */
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, SrcData& src, int i_src, int i_dst) const noexcept
    {
        using namespace amrex;

        constexpr ParticleReal me = PhysConst::m_e;
        constexpr ParticleReal one_over_me = 1._prt/me;

        const ParticleReal w = src.m_rdata[PIdx::w][i_src];
        const ParticleReal ux = src.m_rdata[PIdx::ux][i_src];
        const ParticleReal uy = src.m_rdata[PIdx::uy][i_src];
        const ParticleReal uz = src.m_rdata[PIdx::uz][i_src];
        const ParticleReal ex = src.m_rdata[PIdx::Ex][i_src];
        const ParticleReal ey = src.m_rdata[PIdx::Ey][i_src];
        const ParticleReal ez = src.m_rdata[PIdx::Ez][i_src];
        const ParticleReal bx = src.m_rdata[PIdx::Bx][i_src];
        const ParticleReal by = src.m_rdata[PIdx::By][i_src];
        const ParticleReal bz = src.m_rdata[PIdx::Bz][i_src];

        // Particle momentum is stored as gamma * velocity.
        // Convert to m * gamma * velocity before applying the emission functor.
        auto px = ux*me;
        auto py = uy*me;
        auto pz = uz*me;

        auto g_w = 0.0_rt;
        auto g_px = 0.0_rt;
        auto g_py = 0.0_rt;
        auto g_pz = 0.0_rt;

        m_emission_functor.operator()<1>(
            &px, &py, &pz,
            ex, ey, ez,
            bx, by, bz,
            w,
            &g_px, &g_py, &g_pz,
            &g_w);

        // Then convert back to WarpX convention.
        src.m_rdata[PIdx::ux][i_src] = px*one_over_me;
        src.m_rdata[PIdx::uy][i_src] = py*one_over_me;
        src.m_rdata[PIdx::uz][i_src] = pz*one_over_me;

        dst.m_rdata[PIdx::w][i_dst] = g_w;
        dst.m_rdata[PIdx::ux][i_dst] = g_px*one_over_me;
        dst.m_rdata[PIdx::uy][i_dst] = g_py*one_over_me;
        dst.m_rdata[PIdx::uz][i_dst] = g_pz*one_over_me;

        //Initialize the optical depth component of the source species.
        src.m_runtime_rdata[m_opt_depth_runtime_comp][i_src] =
            m_opt_depth_functor();
    }

private:
    const QuantumSynchrotronGetOpticalDepth
        m_opt_depth_functor;  /*!< A copy of the functor to initialize the optical depth of the source species. */
    const QuantumSynchrotronGeneratePhotonAndUpdateMomentum
        m_emission_functor;  /*!< A copy of the functor to generate photons. It contains only pointers to the lookup tables.*/
    const int m_opt_depth_runtime_comp = 0;  /*!< Index of the optical depth component of source species*/
};


/**
* \brief Free function to call to remove immediately
* low energy photons by setting their ID to -1.
*
* @tparam PTile particle tile type
* @param[in,out] ptile a particle tile
* @param[in] old_size the old number of particles
* @param[in] num_added the number of photons added to the tile
* @param[in] energy_threshold the energy threshold
*/
template <typename PTile>
void cleanLowEnergyPhotons(
    PTile& ptile,
    const int old_size, const int num_added,
    const amrex::ParticleReal energy_threshold)
{
    auto pp = ptile.GetArrayOfStructs()().data() + old_size;

    const auto& soa = ptile.GetStructOfArrays();
    const auto p_ux = soa.GetRealData(PIdx::ux).data() + old_size;
    const auto  p_uy = soa.GetRealData(PIdx::uy).data() + old_size;
    const auto  p_uz = soa.GetRealData(PIdx::uz).data() + old_size;

    //The square of the energy threshold
    const auto energy_threshold2 = energy_threshold*energy_threshold;

    amrex::ParallelFor(num_added, [=] AMREX_GPU_DEVICE (int ip) noexcept
    {
        auto& p = pp[ip];

        const auto ux = p_ux[ip];
        const auto uy = p_uy[ip];
        const auto uz = p_uz[ip];

        //The square of the photon energy (in SI units)
        // ( Particle momentum is stored as gamma * velocity.)
        constexpr amrex::ParticleReal me_c = PhysConst::m_e*PhysConst::c;
        const auto phot_energy2 = (ux*ux + uy*uy + uz*uz)*me_c*me_c;

        if (phot_energy2 < energy_threshold2){
            p.id() = - 1;
        }
    });
}


#endif //QED_PHOTON_EMISSION_H_
