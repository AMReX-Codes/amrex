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
    // control parameters, can be set here or by parm parse
    // variables of the same name prefixed with sensei.
    void setConfig(const std::string &aconfig) { config = aconfig; }
    void setFrequency(int freq) { frequency = freq; }

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
    int frequency;
    int counter;
};

}

#endif
