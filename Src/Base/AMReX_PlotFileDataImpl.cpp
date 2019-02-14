#include <algorithm>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

namespace amrex {

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

PlotFileDataImpl::PlotFileDataImpl (std::string const& plotfile_name)
    : m_plotfile_name(plotfile_name)
{
    // Header
    std::string File(plotfile_name+"/Header");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::istringstream is(std::string(fileCharPtr.dataPtr()), std::istringstream::in);

    is >> m_file_version;

    is >> m_ncomp;
    m_var_names.resize(m_ncomp);
    for (int i = 0; i < m_ncomp; ++i) {
        is >> m_var_names[i];
    }

    is >> m_spacedim >> m_time >> m_finest_level;
    m_nlevels = m_finest_level+1;

    for (int i = 0; i < m_spacedim; ++i) {
        is >> m_prob_lo[i];
    }
    for (int i = 0; i < m_spacedim; ++i) {
        is >> m_prob_hi[i];
        m_prob_size[i] = m_prob_hi[i] - m_prob_lo[i];
    }

    m_ref_ratio.resize(m_nlevels, 0);
    for (int i = 0; i < m_finest_level; ++i) {
        is >> m_ref_ratio[i];
    }
    GotoNextLine(is);

    m_prob_domain.resize(m_nlevels);
    for (int i = 0; i < m_nlevels; ++i) {
        is >> m_prob_domain[i];
    }

    m_level_steps.resize(m_nlevels);
    for (int i = 0; i < m_nlevels; ++i) {
        is >> m_level_steps[i];
    }

    m_cell_size.resize(m_nlevels, Array<Real,AMREX_SPACEDIM>{AMREX_D_DECL(1.,1.,1.)});
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        for (int idim = 0; idim < m_spacedim; ++idim) {
            is >> m_cell_size[ilev][idim];
        }
    }

    is >> m_coordsys;
    int bwidth;
    is >> bwidth;

    m_mf_name.resize(m_nlevels);
    m_vismf.resize(m_nlevels);
    m_ba.resize(m_nlevels);
    m_dmap.resize(m_nlevels);
    m_ngrow.resize(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        int levtmp, ngrids, levsteptmp;
        Real gtime;
        is >> levtmp >> ngrids >> gtime;
        is >> levsteptmp;
        Real glo[3], ghi[3];
        for (int igrid = 0; igrid < ngrids; ++igrid) {
            for (int idim = 0; idim < m_spacedim; ++idim) {
                is >> glo[idim] >> ghi[idim];
            }
        }
        std::string relname;
        is >> relname;
        m_mf_name[ilev] = m_plotfile_name + "/" + relname;
        if (m_ncomp > 0) {
            m_vismf[ilev].reset(new VisMF(m_mf_name[ilev]));
            m_ba[ilev] = m_vismf[ilev]->boxArray();
            m_dmap[ilev].define(m_ba[ilev]);
            m_ngrow[ilev] = m_vismf[ilev]->nGrowVect();
        }
    }
}

PlotFileDataImpl::~PlotFileDataImpl () {}

void
PlotFileDataImpl::syncDistributionMap (PlotFileDataImpl const& src)
{
    int nlevs_min = std::min(m_nlevels, src.m_nlevels);
    for (int ilev = 0; ilev < nlevs_min; ++ilev) {
        syncDistributionMap(ilev, src);
    }
}

void
PlotFileDataImpl::syncDistributionMap (int level, PlotFileDataImpl const& src)
{
    if (level <= src.finestLevel() and m_dmap[level].size() == src.DistributionMap(level).size()) {
        m_dmap[level] = src.DistributionMap(level);
    }
}

MultiFab
PlotFileDataImpl::get (int level)
{
    MultiFab mf(m_ba[level], m_dmap[level], m_ncomp, m_ngrow[level]);
    VisMF::Read(mf, m_mf_name[level]);
    return mf;
}

MultiFab
PlotFileDataImpl::get (int level, std::string const& varname)
{
    MultiFab mf(m_ba[level], m_dmap[level], 1, m_ngrow[level]);
    auto r = std::find(std::begin(m_var_names), std::end(m_var_names), varname);
    if (r == std::end(m_var_names)) {
        amrex::Abort("PlotFileDataImpl::get: varname not found "+varname);
    } else {
        int icomp = std::distance(std::begin(m_var_names), r);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            int gid = mfi.index();
            FArrayBox& dstfab = mf[mfi];
            std::unique_ptr<FArrayBox> srcfab(m_vismf[level]->readFAB(gid, icomp));
            dstfab.copy(*srcfab);
        }
    }
    return mf;
}

}
