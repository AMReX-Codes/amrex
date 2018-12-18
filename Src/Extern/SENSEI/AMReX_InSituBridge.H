#ifndef AMReX_InSituBridge_H
#define AMReX_InSituBridge_H

#include <string>
#include <AMReX_Vector.H>

#if defined(BL_USE_SENSEI_INSITU)
#include <mpi.h>
namespace sensei {
class AnalysisAdaptor;
}
#endif

namespace amrex
{

/// Base class for coupling to in situ. implements initialze and finalize
/// derived classes implement an update method. see AmrInSituBridge and
/// AmrMeshInSituBridge.
class InSituBridge
{
public:
    InSituBridge();
    virtual ~InSituBridge();

    InSituBridge(const InSituBridge&) = delete;
    void operator=(const InSituBridge&) = delete;

#if defined(BL_USE_MPI)
    void setCommuniator(MPI_Comm acomm) { comm = acomm; }
#endif

    // Runtime control parameters.
    //
    // These may set explicitly or via ParmParse. ParmParse values
    // take precedence over set calls.
    //
    //   sensei.enabled  -- integer 0 or 1, when 0 the bridge is disabled
    //   sensei.frequency -- integer >= 1, how often the analysis executes
    //   sensei.config -- file path, path to the SENSEI XML config
    //   sensei.pin_mesh -- integer, 0 or 1, when 1 mesh origin is forced to 0,0,0
    //
    void setEnabled(int aenabled) { enabled = aenabled; }
    void setFrequency(int freq) { frequency = freq; }
    void setConfig(const std::string &aconfig) { config = aconfig; }
    void setPinMesh(int apinMesh) { pinMesh = apinMesh; }

    // initialize the adaptors
    // will make use of parm parse for any unset
    // control parameters
    virtual int initialize();

    // clean up the adaptors
    virtual int finalize();

protected:
    bool doUpdate();

protected:
#if defined(BL_USE_MPI)
    MPI_Comm comm;
#endif
#if defined(BL_USE_SENSEI_INSITU)
    sensei::AnalysisAdaptor *analysis_adaptor;
#endif
    std::string config;
    int enabled;
    int frequency;
    int counter;
    int pinMesh;
};

}

#endif
