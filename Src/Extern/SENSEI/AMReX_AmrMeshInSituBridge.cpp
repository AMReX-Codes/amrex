#include <AMReX_AmrMeshInSituBridge.H>

#include <AMReX_ParmParse.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#ifdef BL_USE_SENSEI_INSITU
#include <chrono>
#include <AnalysisAdaptor.h>
#include <Profiler.h>
#include <AMReX_AmrMeshDataAdaptor.H>
#endif

namespace amrex {

int
AmrMeshInSituBridge::update(unsigned int step, double time,
    amrex::AmrMesh *mesh, const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
    const std::vector<std::vector<std::string>> &names)
{
    int ret = 0;
#if defined(BL_USE_SENSEI_INSITU)
    if (doUpdate())
    {
        amrex::Print() << "SENSEI Begin update..." << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        sensei::TimeEvent<64> event("AmrMeshInSituBridge::update");

        amrex::AmrMeshDataAdaptor *data_adaptor = amrex::AmrMeshDataAdaptor::New();
        if (comm != MPI_COMM_NULL)
            data_adaptor->SetCommunicator(comm);
        data_adaptor->SetPinMesh(pinMesh);
        data_adaptor->SetDataSource(mesh, states, names);
        data_adaptor->SetDataTime(time);
        data_adaptor->SetDataTimeStep(step);
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
