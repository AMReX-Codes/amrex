#include <iostream>

#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include <Geometry.H>

#include <InTransitAnalysis.H>

void InTransitAnalysis::Initialize (MultiFab &mf_in, Geometry &geom_in, int time_step_in)
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "InTransitAnalysis class is initializing ..." << std::endl;

    mf = &mf_in;
    geom = &geom_in;
    time_step = time_step_in;
};

void InTransitAnalysis::DoAnalysis ()
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "InTransitAnalysis class is doing analysis ..." << std::endl;

    norm0 = mf->norm0();
    probsize = geom->ProbSize();
    if (ParallelDescriptor::IOProcessor())
        PrintResults();
};

void InTransitAnalysis::Finalize()
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "InTransitAnalysis class is finalizing ..." << std::endl;
};

void InTransitAnalysis::PrintResults() const
{
    std::cout << "L0 norm is " << norm0 << std::endl;
    std::cout << "problem size is " << probsize << std::endl;
    std::cout << "time step is " << time_step << std::endl;
};
