#include <AMReX_InSituBridge.H>

#include <AMReX_ParmParse.H>

#ifdef BL_USE_SENSEI_INSITU
#include <chrono>
#include <DataAdaptor.h>
#include <AnalysisAdaptor.h>
#include <ConfigurableAnalysis.h>
#include <Profiler.h>
#include <AMReX_AmrDataAdaptor.H>
#include <AMReX_AmrMeshDataAdaptor.H>
#endif

namespace amrex {

InSituBridge::InSituBridge() :
#if defined(BL_USE_MPI)
    comm(MPI_COMM_NULL),
#endif
#if defined(BL_USE_SENSEI_INSITU)
    analysis_adaptor(nullptr),
#endif
    enabled(0), frequency(1), counter(0), pinMesh(0)
{
#if defined(BL_USE_SENSEI_INSITU)
    sensei::Profiler::Initialize();
    sensei::Profiler::StartEvent("InSituBridge::LifeTime");
#endif
}

InSituBridge::~InSituBridge()
{
#if defined(BL_USE_SENSEI_INSITU)
    if (analysis_adaptor)
        analysis_adaptor->Delete();
    sensei::Profiler::EndEvent("InSituBridge::LifeTime");
    sensei::Profiler::Finalize();
#endif
}

int
InSituBridge::initialize()
{
#if defined(BL_USE_SENSEI_INSITU)
    auto t0 = std::chrono::high_resolution_clock::now();
    sensei::TimeEvent<64> event("InSituBridge::initialize");

    // read config from ParmParse
    ParmParse pp("sensei");

    pp.query("enabled", enabled);

    if (!enabled)
        return 0;

    pp.query("config", config);
    pp.query("frequency", frequency);
    pp.query("pin_mesh", pinMesh);

    amrex::Print() << "SENSEI Begin initialize..." << std::endl;

    // Check for invalid values
    if (config.empty())
    {
        amrex::ErrorStream()
            << "Error: Missing SENSEI XML configuration. To correct add "
               "\"sensei.config=/path/to/file.xml\" to your inputs file."
            << std::endl;
        return -1;
    }

    if (frequency < 1)
    {
        amrex::ErrorStream()
            << "Error: sensei.frequency=" << frequency
            << " Frequency must be greater or equal to 1. "
            << std::endl;
        return -1;
    }

    // create and initialize the analysis adaptor
    sensei::ConfigurableAnalysis *aa = sensei::ConfigurableAnalysis::New();

#if defined(BL_USE_MPI)
    if (comm != MPI_COMM_NULL)
        aa->SetCommunicator(comm);
#endif

    if (aa->Initialize(config))
    {
        aa->Delete();
        aa = nullptr;
        return -1;
    }

    analysis_adaptor = aa;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "SENSEI initialize complete (" << dt.count() << " sec)" << std::endl;
#endif
    return 0;
}

bool
InSituBridge::doUpdate()
{
    bool ret = analysis_adaptor && (frequency > 0) && ((counter % frequency) == 0);
    counter += 1;
    return ret;
}

int
InSituBridge::finalize()
{
    int ret = 0;
#if defined(BL_USE_SENSEI_INSITU)
    if (!analysis_adaptor)
        return ret;

    amrex::Print() << "SENSEI Begin finalize..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    sensei::TimeEvent<64> event("InSituBridge::finalize");
    ret = analysis_adaptor->Finalize();

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "SENSEI finalize complete (" << dt.count() << " sec)" << std::endl;
#endif
    return ret;
}

}
