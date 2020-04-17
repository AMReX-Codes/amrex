#include "FlushFormatOpenPMD.H"
#include "WarpX.H"
#include "Utils/Interpolate.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

FlushFormatOpenPMD::FlushFormatOpenPMD (const std::string& diag_name)
{
    ParmParse pp(diag_name);
    pp.query("openpmd_backend", openpmd_backend);
    pp.query("openpmd_tspf", openpmd_tspf);
    auto & warpx = WarpX::GetInstance();
    m_OpenPMDPlotWriter = new WarpXOpenPMDPlot(
        openpmd_tspf, openpmd_backend, warpx.getPMLdirections()
        );
}

void
FlushFormatOpenPMD::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<const amrex::MultiFab*> mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    MultiParticleContainer& mpc, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_raw_rho, bool plot_raw_F) const
{
    WARPX_PROFILE("FlushFormatOpenPMD::WriteToFile()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !plot_raw_fields && !plot_raw_fields_guards && !plot_raw_rho && !plot_raw_F,
        "Cannot plot raw data with OpenPMD output format. Use plotfile instead.");

    // Set step and output directory name.
    m_OpenPMDPlotWriter->SetStep(iteration[0], prefix);

    // fields: only dumped for coarse level
    m_OpenPMDPlotWriter->WriteOpenPMDFields(
        varnames, *mf[0], geom[0], iteration[0], time);

    // particles: all (reside only on locally finest level)
    m_OpenPMDPlotWriter->WriteOpenPMDParticles(mpc);
}

FlushFormatOpenPMD::~FlushFormatOpenPMD (){
    delete m_OpenPMDPlotWriter;
}
