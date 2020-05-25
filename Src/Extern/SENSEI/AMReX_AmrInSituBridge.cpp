#include <AMReX_AmrInSituBridge.H>

#include <AMReX_Amr.H>
#ifdef BL_USE_SENSEI_INSITU
#include <chrono>
#include <AnalysisAdaptor.h>
#include <Profiler.h>
#include <AMReX_AmrDataAdaptor.H>
#endif

namespace amrex {

int
AmrInSituBridge::update(Amr *dataSource)
{
    int ret = 0;
#if defined(BL_USE_SENSEI_INSITU)
    if (doUpdate())
    {
        amrex::Print() << "SENSEI Begin update..." << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        sensei::TimeEvent<64> event("AMRInSituBridge::update");

        amrex::AmrDataAdaptor *data_adaptor = amrex::AmrDataAdaptor::New();
        if (comm != MPI_COMM_NULL)
            data_adaptor->SetCommunicator(comm);
        data_adaptor->SetPinMesh(pinMesh);
        data_adaptor->SetDataSource(dataSource);
        data_adaptor->SetDataTime(dataSource->cumTime());
        data_adaptor->SetDataTimeStep(dataSource->levelSteps(0));
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
