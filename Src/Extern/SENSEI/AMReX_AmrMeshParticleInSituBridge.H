#ifndef AMReX_AmrMeshParticleInSituBridge_H
#define AMReX_AmrMeshParticleInSituBridge_H
#include <AMReX_Config.H>

#ifdef AMREX_PARTICLES

#include <AMReX_InSituBridge.H>
#include <AMReX_Vector.H>

#include <AMReX_AmrParticles.H>

#include <AMReX_AmrMeshParticleDataAdaptor.H>

#ifdef AMREX_USE_SENSEI_INSITU
#include <AnalysisAdaptor.h>
#include <Profiler.h>
#endif

#include <chrono>

namespace amrex
{

/// SENSEI bridge code for simulations proccessing both amrex::Amr and amrex::ParticleContainer
class AmrMeshParticleInSituBridge : public InSituBridge
{
public:
  AmrMeshParticleInSituBridge() {}
  ~AmrMeshParticleInSituBridge() {}

  AmrMeshParticleInSituBridge(const AmrMeshParticleInSituBridge&) = delete;
  void operator=(const AmrMeshParticleInSituBridge&) = delete;

  // invoke the in situ analysis with data from an AmrMesh and ParticleContainer
  // instance.
  template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
  int update(
    long step, double time,
    amrex::AmrMesh *mesh,
    const std::vector<amrex::Vector<amrex::MultiFab>*> &mesh_states,
    const std::vector<std::vector<std::string>> &mesh_names,
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
    const std::map<std::string, std::vector<int>> & particles_rStructs = {},
    const std::map<std::string, int> & particles_iStructs = {},
    const std::map<std::string, std::vector<int>> & particles_rArrays = {},
    const std::map<std::string, int> & particles_iArrays = {});
};

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrMeshParticleInSituBridge::update(
  long step, double time,
  amrex::AmrMesh *mesh,
  const std::vector<amrex::Vector<amrex::MultiFab>*> &mesh_states,
  const std::vector<std::vector<std::string>> &mesh_names,
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
    amrex::Print() << "SENSEI AmrMesh+Particles Begin update..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    sensei::TimeEvent<64> event("AmrMeshParticleInSituBridge::update");

    if(!particles)
    {
      SENSEI_ERROR("no particles presented at update call");
      return -1;
    }

    amrex::AmrMeshParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt> *data_adaptor
      = amrex::AmrMeshParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New();

    if (comm != MPI_COMM_NULL)
      data_adaptor->SetCommunicator(comm);

    data_adaptor->SetDataSource(
      mesh, mesh_states, mesh_names,
      particles, particles_rStructs, particles_iStructs,
      particles_rArrays, particles_iArrays);

    data_adaptor->SetDataTime(time);
    data_adaptor->SetDataTimeStep(step);
    ret = analysis_adaptor->Execute(data_adaptor, nullptr) ? 0 : -1;
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

#endif // AMREX_PARTICLES
#endif // AMReX_AmrMeshParticleInSituBridge_H
