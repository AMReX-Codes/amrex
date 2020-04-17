#ifndef WARPX_FLUSHFORMATOPENPMD_H_
#define WARPX_FLUSHFORMATOPENPMD_H_

#include "FlushFormat.H"
#include "Diagnostics/WarpXOpenPMD.H"

/**
 * \brief This class aims at dumping diags data to disk using the OpenPMD standard.
 * In particular, function WriteToFile takes fields and particles as input arguments,
 * and writes data to file.
 */
class FlushFormatOpenPMD : public FlushFormat
{
public:

    /** Constructor takes name of diagnostics to set the output directory */
    FlushFormatOpenPMD (const std::string& diag_name);

    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<const amrex::MultiFab*> mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        MultiParticleContainer& mpc, int nlev, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F) const override final;

    ~FlushFormatOpenPMD ();

private:
    /** This is responsible for dumping to file */
    WarpXOpenPMDPlot* m_OpenPMDPlotWriter = nullptr;
};

#endif // WARPX_FLUSHFORMATOPENPMD_H_
