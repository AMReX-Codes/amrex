#ifndef AMREX_EB_CHKPT_FILE_H_
#define AMREX_EB_CHKPT_FILE_H_

#include <AMReX_MultiFab.H>

namespace amrex { namespace EB2 {

class ChkptFile
{
private:
    std::string m_restart_file = "";

public:
    ChkptFile(const std::string &fname);

    void fill_from_chkpt_file(BoxArray& ba, DistributionMapping& dmap,
            MultiFab& volfrac, MultiFab& centroid, MultiFab& bndryarea,
            MultiFab& bndrycent,
            MultiFab& bndrynorm, Array<MultiFab,AMREX_SPACEDIM>& areafrac,
            Array<MultiFab,AMREX_SPACEDIM>& facecent,
            Array<MultiFab,AMREX_SPACEDIM>& edgecent,
            MultiFab& levelset, int ng) const;
};

}}

#endif
