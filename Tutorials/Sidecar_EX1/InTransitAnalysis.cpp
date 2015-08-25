#include <iostream>

#include <MultiFab.H>
#include <ParallelDescriptor.H>

#include <InTransitAnalysis.H>

InTransitAnalysis::InTransitAnalysis (MultiFab &mf_in, Geometry &geom_in, int time_step_in)
    : mf(&mf_in),
      geom(&geom_in),
      time_step(time_step_in)
{
};

void InTransitAnalysis::Initialize ()
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "InTransitAnalysis class is initializing ..." << std::endl;
};

void InTransitAnalysis::DoAnalysis ()
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "InTransitAnalysis class is doing analysis ..." << std::endl;

    norm0 = mf->norm0();
    probsize = geom->ProbSize();
};

void InTransitAnalysis::Finalize()
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "InTransitAnalysis class is finalizing ..." << std::endl;
        PrintResults();
    }
};

void InTransitAnalysis::PrintResults() const
{
    std::cout << "L0 norm is " << norm0 << std::endl;
    std::cout << "problem size is " << probsize << std::endl;
    std::cout << "time step is " << time_step << std::endl;
};
