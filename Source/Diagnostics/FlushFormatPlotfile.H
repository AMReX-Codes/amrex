#ifndef WARPX_FLUSHFORMATPLOTFILE_H_
#define WARPX_FLUSHFORMATPLOTFILE_H_

#include "FlushFormat.H"

class FlushFormatPlotfile : public FlushFormat
{
public:
    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (const amrex::Vector<std::string> varnames,
                              const amrex::Vector<const amrex::MultiFab*> mf,
                              amrex::Vector<amrex::Geometry>& geom,
                              const amrex::Vector<int> iteration, const double time,
                              MultiParticleContainer& mpc, int nlev) const override final;

    /** Write general info of the run into the plotfile */
    void WriteJobInfo(const std::string& dir) const;
    /** Write WarpX-specific plotfile header */
    void WriteWarpXHeader(const std::string& name) const;

};

#endif // WARPX_FLUSHFORMATPLOTFILE_H_
