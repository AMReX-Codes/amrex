#ifndef AMREX_PLOTFILE_DATA_H_
#define AMREX_PLOTFILE_DATA_H_

#include <AMReX_DataServices.H>
#include <string>
#include <memory>

namespace amrex {

class PlotfileData
{
public:

    PlotfileData (std::string const& plotfile_name);

    ~PlotfileData ();

    int spaceDim () const;

    Real time () const;

    int finestLevel () const;

    IntVect refRatio (int level) const;

    const BoxArray& boxArray (int level) const;

    const DistributionMapping& DistributionMap (int level) const;

    int coordSys () const;

    Box ProbDomain (int level) const;

    Array<Real,AMREX_SPACEDIM> probSize () const;
    Array<Real,AMREX_SPACEDIM> probLo () const;
    Array<Real,AMREX_SPACEDIM> probHi () const;
    Array<Real,AMREX_SPACEDIM> cellSize (int level) const;

    const Vector<std::string>& varNames () const;

    int nComp () const;
    IntVect nGrowVect () const;

private:

    std::unique_ptr<DataServices> m_ds;
    AmrData* m_ad;
};

}

#endif
