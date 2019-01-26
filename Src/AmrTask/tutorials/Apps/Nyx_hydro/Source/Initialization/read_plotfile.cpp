//
// This is the version that reads input from PlotFiles.
//
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>

#ifndef WIN32
#include <unistd.h>
#endif

#include "Nyx.H"
#include "Nyx_F.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"

using namespace amrex;

void
Nyx::ReadPlotFile (bool               first,
                   const std::string& file, bool& rhoe_infile)
{
    amrex::Print() << "Reading data from plotfile: " << file << std::endl;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(file, fileType);

    if ( ! dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData&                 amrData    = dataServices.AmrDataRef();
    const Vector<std::string> plotnames  = amrData.PlotVarNames();

    // Sanity checks
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (amrData.ProbLo()[i] != geom.ProbLo()[i])
        {
            std::cout << "AmrData.prob_lo(i) = " << amrData.ProbLo()[i] << std::endl;
            std::cout << "   geom.prob_lo(i) = " <<    geom.ProbLo()[i] << std::endl;
            amrex::Error("prob_lo from plotfile doesn't match prob_lo from inputs file");
        }
        if (amrData.ProbHi()[i] != geom.ProbHi()[i])
        {
            std::cout << "AmrData.prob_hi(i) = " << amrData.ProbHi()[i] << std::endl;
            std::cout << "   geom.prob_hi(i) = " <<    geom.ProbHi()[i] << std::endl;
            amrex::Error("prob_hi from plotfile doesn't match prob_hi from inputs file");
        }
    }

    if (amrData.ProbDomain()[level] != parent->getLevel(level).Domain())
    {
        std::cout << "AmrData.domain = " << amrData.ProbDomain()[level] << std::endl;
        std::cout << "   geom.domain = " << parent->getLevel(level).Domain() << std::endl;
        amrex::Error("prob_domain from plotfile doesn't match prob_domain from inputs file");
    }

    if (amrData.boxArray(level) != grids)
    {
        std::cout << "AmrData.boxArray = " << amrData.boxArray(level) << std::endl;
        std::cout << "   grids         = " << grids << std::endl;
        amrex::Error("boxArray from plotfile doesn't match grids ");
    }

    if (amrData.FinestLevel() != parent->finestLevel()) {
        amrex::Error("finest_level from plotfile doesn't match finest_level from inputs file");
    }

    const int Nlev = parent->finestLevel() + 1;

#ifndef NO_HYDRO
    int iURHO = -1, iUMX = -1, iUMY = -1, iUMZ = -1, iTEMP = -1, iNE = -1, iURHOE = -1;

    for (int i = 0; i < plotnames.size(); ++i)
    {
        if (plotnames[i] == "density")                  iURHO  = i;
        if (plotnames[i] == "xmom")                     iUMX   = i;
        if (plotnames[i] == "ymom")                     iUMY   = i;
        if (plotnames[i] == "zmom")                     iUMZ   = i;
        if (plotnames[i] == "rho_e")                    iURHOE = i;
        if (plotnames[i] == "Temp")                     iTEMP  = i;
        if (plotnames[i] == "Ne")                       iNE    = i;
    }

    if ( iURHO < 0 ) amrex::Abort("\"density\" is not in the plotfile");
    if ( iUMX  < 0 ) amrex::Abort("\"xmom\" is not in the plotfile");
    if ( iUMY  < 0 ) amrex::Abort("\"ymom\" is not in the plotfile");
    if ( iUMZ  < 0 ) amrex::Abort("\"zmom\" is not in the plotfile");

    if ( iURHOE < 0 )
    {
        if ( iTEMP < 0 ) amrex::Abort("\"rho_e\" nor \"Temp\" is in the plotfile");
        if ( iNE   < 0 ) amrex::Abort("\"rho_e\" nor \"Ne\" is in the plotfile");

        if (iNE != iTEMP + 1) amrex::Abort("We assume Ne = Temp +1");
    }
    else {rhoe_infile = true;}

    Real time = amrData.Time();
    parent->setCumTime(time);

    // Note that we only store this for level 0.
    nsteps_from_plotfile = amrData.LevelSteps()[0];

    //
    // Read density and momentum
    //
    for (int lev = 0; lev < Nlev; ++lev)
    {

        MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);
        const Box bx = grids.minimalBox();

        S_new.copy(amrData.GetGrids(lev,iURHO,bx),0,Density,1);
        amrData.FlushGrids(iURHO);

        S_new.copy(amrData.GetGrids(lev,iUMX,bx),0,Xmom,1);
        amrData.FlushGrids(iUMX);

        S_new.copy(amrData.GetGrids(lev,iUMY,bx),0,Ymom,1);
        amrData.FlushGrids(iUMY);

        S_new.copy(amrData.GetGrids(lev,iUMZ,bx),0,Zmom,1);
        amrData.FlushGrids(iUMZ);

        if (rhoe_infile)
        {
            S_new.copy(amrData.GetGrids(lev,iURHOE,bx),0,Eint,1);
            amrData.FlushGrids(iURHOE);
        }
    }

    amrex::Print() << "Successfully read state data" << std::endl;

    //
    // Read temperature and Ne if there is no rho_e in the file
    //
    if ( ! rhoe_infile)
    {
        for (int lev = 0; lev < Nlev; ++lev)
        {
            MultiFab& D_new = parent->getLevel(lev).get_new_data(DiagEOS_Type);
            const Box bx = grids.minimalBox();

            D_new.copy(amrData.GetGrids(lev,iTEMP,bx),0,Temp_comp,1);
            amrData.FlushGrids(iTEMP);

            D_new.copy(amrData.GetGrids(lev,iNE,bx),0,Ne_comp,1);
            amrData.FlushGrids(iNE);

            amrex::Print() << "D_new.max " << D_new.norm0() << std::endl;;
        }

        amrex::Print() << "Successfully read temperature and Ne" << std::endl;
    }
#endif
}
