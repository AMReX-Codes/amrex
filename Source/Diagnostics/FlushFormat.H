#ifndef WARPX_FLUSHFORMAT_H_
#define WARPX_FLUSHFORMAT_H_

#include "Particles/MultiParticleContainer.H"

class FlushFormat
{
public:
    /** Flush fields and particles to file */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<const amrex::MultiFab*> mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        MultiParticleContainer& mpc, int nlev,
        const std::string prefix) const = 0;
};

#endif // WARPX_FLUSHFORMAT_H_
