#include <iostream>

#include <MultiFab.H>

#include <InSituAnalysis.H>

InSituAnalysis::InSituAnalysis(MultiFab& mf_in, Geometry &geom_in)
    : mf(&mf_in),
      geom(&geom_in)
{};

void InSituAnalysis::DoAnalysis ()
{
    norm0 = mf->norm0();
    probsize = geom->ProbSize();
};

void InSituAnalysis::PrintResults() const
{
    std::cout << "========== BEGIN ANALYSIS ==========" << std::endl;
    std::cout << "L0 norm is " << norm0 << std::endl;
    std::cout << "problem size is " << probsize << std::endl;
    std::cout << "========== END ANALYSIS ==========" << std::endl;
};
