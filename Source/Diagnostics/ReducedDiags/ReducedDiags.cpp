#include "ReducedDiags.H"
#include "AMReX_ParmParse.H"

using namespace amrex;

ReducedDiags::ReducedDiags (std::string rd_name)
{

    /** read reduced diags frequency */
    ParmParse pp1(rd_name);
    pp1.query("frequency", m_freq);

    /** read number species */
    ParmParse pp2("particles");
    pp2.query("nspecies", m_nspecies);

}

ReducedDiags::~ReducedDiags ()
{}
