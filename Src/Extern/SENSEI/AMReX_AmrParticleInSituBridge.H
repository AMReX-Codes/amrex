#ifndef AMReX_AmrParticleInSituBridge_H
#define AMReX_AmrParticleInSituBridge_H

#include <AMReX_Config.H>
#ifdef AMREX_PARTICLES
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Amr.H>
#ifdef AMREX_USE_SENSEI_INSITU
#include <AMReX_AmrParticleDataAdaptor.H>
#include <AMReX_InSituBridge.H>
#include <AnalysisAdaptor.h>
#include <Profiler.h>
#include <chrono>
#endif

namespace amrex
{

/// SENSEI bridge code for simulations processing both amrex::Amr and amrex::ParticleContainer
class AmrParticleInSituBridge : public InSituBridge
{
public:
  AmrParticleInSituBridge() {}
  ~AmrParticleInSituBridge() {}

  AmrParticleInSituBridge(const AmrParticleInSituBridge&) = delete;
  void operator=(const AmrParticleInSituBridge&) = delete;

  /** invoke the in situ analysis with data from an AmrMesh and ParticleContainer
   * instance.
   */
  template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
  int update(amrex::Amr *amr,
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
    const std::map<std::string, std::vector<int>> & particles_rStructs = {},
    const std::map<std::string, int> & particles_iStructs = {},
    const std::map<std::string, std::vector<int>> & particles_rArrays = {},
    const std::map<std::string, int> & particles_iArrays = {});
};

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleInSituBridge::update(
  amrex::Amr *amr,
  amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
  const std::map<std::string, std::vector<int>> & particles_rStructs,
  const std::map<std::string, int> & particles_iStructs,
  const std::map<std::string, std::vector<int>> & particles_rArrays,
  const std::map<std::string, int> & particles_iArrays)
{
  int ret = 0;
#if defined(AMREX_USE_SENSEI_INSITU)
  if (doUpdate())
  {
    sensei::TimeEvent<64> event("AmrParticleInSituBridge::update");

    amrex::Print() << "SENSEI Particles Begin update..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    if(!particles)
    {
      SENSEI_ERROR("no particles presented at update call");
      return -1;
    }

    amrex::AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt> *data_adaptor
      = amrex::AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New();

    if (comm != MPI_COMM_NULL)
      data_adaptor->SetCommunicator(comm);

    data_adaptor->SetDataSource(amr, particles, particles_rStructs,
      particles_iStructs, particles_rArrays, particles_iArrays);

    data_adaptor->SetDataTime(amr->cumTime());
    data_adaptor->SetDataTimeStep(amr->levelSteps(0));
    ret = analysis_adaptor->Execute(data_adaptor) ? 0 : -1;
    data_adaptor->ReleaseData();
    data_adaptor->Delete();

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "SENSEI update complete (" << dt.count() << " sec)" << std::endl;
  }
#endif
  return ret;
}
}
#endif
#endif
