#ifndef AMReX_ParticleInSituBridge_H
#define AMReX_ParticleInSituBridge_H

#include <AMReX_Config.H>
#ifdef AMREX_PARTICLES
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#ifdef AMREX_USE_SENSEI_INSITU
#include <AMReX_InSituBridge.H>
#include <AMReX_ParticleDataAdaptor.H>
#include <AnalysisAdaptor.h>
#include <Profiler.h>
#include <chrono>
#endif

namespace amrex
{
/// Bridge for simulations that use amrex::ParticleContainer
/// inherits initialize and finalize from amrex::InSituBridge
class ParticleInSituBridge : public InSituBridge
{
public:
  ParticleInSituBridge() {}
  ~ParticleInSituBridge() {}

  ParticleInSituBridge(const ParticleInSituBridge&) = delete;
  void operator=(const ParticleInSituBridge&) = delete;

  // invoke the in situ analysis with data from an AmrMesh
  // instance.
  template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
  int update(
    const double time,
    const int time_step,
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> *particles,
    const std::map<std::string, std::vector<int>> & particles_rStructs = {},
    const std::map<std::string, int> & particles_iStructs = {},
    const std::map<std::string, std::vector<int>> & particles_rArrays = {},
    const std::map<std::string, int> & particles_iArrays = {});
};

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleInSituBridge::update(
  const double time,
  const int time_step,
  amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> *particles,
  const std::map<std::string, std::vector<int>> & particles_rStructs,
  const std::map<std::string, int> & particles_iStructs,
  const std::map<std::string, std::vector<int>> & particles_rArrays,
  const std::map<std::string, int> & particles_iArrays)
{
  int ret = 0;
#if defined(AMREX_USE_SENSEI_INSITU)
  if (doUpdate())
  {
    sensei::TimeEvent<64> event("ParticleInSituBridge::update");

    amrex::Print() << "SENSEI Particles Begin update..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    if (!particles)
    {
      SENSEI_ERROR("no particles presented at update call");
      return -1;
    }

    ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt> *data_adaptor
      = ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New();

    if (comm != MPI_COMM_NULL)
        data_adaptor->SetCommunicator(comm);

    data_adaptor->SetDataSource(particles, particles_rStructs,
      particles_iStructs, particles_rArrays, particles_iArrays);

    data_adaptor->SetDataTime(time);
    data_adaptor->SetDataTimeStep(time_step);
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
