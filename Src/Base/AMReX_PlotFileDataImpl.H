#ifndef AMREX_PLOT_FILE_DATA_IMPL_H_
#define AMREX_PLOT_FILE_DATA_IMPL_H_

#include <string>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

namespace amrex {

class PlotFileDataImpl
{
public:
    PlotFileDataImpl (std::string const& plotfile_name);
    ~PlotFileDataImpl ();

    int spaceDim () const noexcept { return m_spacedim; }

    Real time () const noexcept { return m_time; }

    int finestLevel () const noexcept { return m_finest_level; }

    int refRatio (int level) const noexcept { return m_ref_ratio[level]; }

    int levelStep (int level) const noexcept { return m_level_steps[level]; }

    const BoxArray& boxArray (int level) const noexcept { return m_ba[level]; }

    const DistributionMapping& DistributionMap (int level) const noexcept { return m_dmap[level]; }

    void syncDistributionMap (PlotFileDataImpl const& src) noexcept;

    void syncDistributionMap (int level, PlotFileDataImpl const& src) noexcept;

    int coordSys () const noexcept { return m_coordsys; }

    Box probDomain (int level) const noexcept { return m_prob_domain[level]; }

    Array<Real,AMREX_SPACEDIM> probSize () const noexcept { return m_prob_size; }
    Array<Real,AMREX_SPACEDIM> probLo () const noexcept { return m_prob_lo; }
    Array<Real,AMREX_SPACEDIM> probHi () const noexcept { return m_prob_hi; }
    Array<Real,AMREX_SPACEDIM> cellSize (int level) const noexcept { return m_cell_size[level]; }

    const Vector<std::string>& varNames () const noexcept { return m_var_names; }

    int nComp () const noexcept { return m_ncomp; }
    IntVect nGrowVect (int level) const noexcept { return m_ngrow[level]; }

    MultiFab get (int level) noexcept;
    MultiFab get (int level, std::string const& varname) noexcept;

private:
    std::string m_plotfile_name;
    std::string m_file_version;
    int m_ncomp;
    Vector<std::string> m_var_names;
    int m_spacedim;
    Real m_time;
    int m_finest_level, m_nlevels;
    Array<Real,AMREX_SPACEDIM> m_prob_lo {{AMREX_D_DECL(0.,0.,0.)}};
    Array<Real,AMREX_SPACEDIM> m_prob_hi {{AMREX_D_DECL(1.,1.,1.)}};
    Array<Real,AMREX_SPACEDIM> m_prob_size {{AMREX_D_DECL(1.,1.,1.)}};
    Vector<int> m_ref_ratio;
    Vector<Box> m_prob_domain;
    Vector<int> m_level_steps;
    Vector<Array<Real,AMREX_SPACEDIM> > m_cell_size;
    int m_coordsys;
    Vector<std::string> m_mf_name;
    Vector<std::unique_ptr<VisMF> > m_vismf;
    Vector<BoxArray> m_ba;
    Vector<DistributionMapping> m_dmap;
    Vector<IntVect> m_ngrow;
};

}

#endif
