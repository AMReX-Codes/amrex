#ifndef AMREX_EB_CHKPT_FILE_H_
#define AMREX_EB_CHKPT_FILE_H_

#include <AMReX_MultiFab.H>

namespace amrex::EB2 {

class ChkptFile
{
private:
    std::string m_restart_file = "";

    const std::string m_volfrac_name   = "volfrac";
    const std::string m_centroid_name  = "centroid";
    const std::string m_bndryarea_name = "bndryarea";
    const std::string m_bndrycent_name = "bndrycent";
    const std::string m_bndrynorm_name = "bndrynorm";
    const std::string m_levelset_name  = "levelset";

    const amrex::Vector<std::string> m_areafrac_name
        = {AMREX_D_DECL("areafrac_x", "areafrac_y", "areafrac_z")};
    const amrex::Vector<std::string> m_facecent_name
        = {AMREX_D_DECL("facecent_x", "facecent_y", "facecent_z")};
    const amrex::Vector<std::string> m_edgecent_name
        = {AMREX_D_DECL("edgecent_x", "edgecent_y", "edgecent_z")};

    void writeHeader (const BoxArray& cut_ba, const BoxArray& covered_ba, const Geometry& geom,
                      const IntVect& ngrow, bool extend_domain_face, int max_grid_size) const;

    void writeToFile (const MultiFab& mf, const std::string& mf_name) const;


public:
    ChkptFile (std::string fname);

    void read_from_chkpt_file (BoxArray& cut_grids, BoxArray& covered_grids,
                               DistributionMapping& dmap,
                               MultiFab& volfrac, MultiFab& centroid, MultiFab& bndryarea,
                               MultiFab& bndrycent, MultiFab& bndrynorm,
                               Array<MultiFab,AMREX_SPACEDIM>& areafrac,
                               Array<MultiFab,AMREX_SPACEDIM>& facecent,
                               Array<MultiFab,AMREX_SPACEDIM>& edgecent,
                               MultiFab& levelset, int ng_gfab, const Geometry& geom,
                               const IntVect& ngrow_finest, bool extend_domain_face, int max_grid_size) const;

    void write_to_chkpt_file (const BoxArray& cut_grids,
                              const BoxArray& covered_grids,
                              const MultiFab& volfrac,
                              const MultiFab& centroid, const MultiFab& bndryarea,
                              const MultiFab& bndrycent, const MultiFab& bndrynorm,
                              const Array<MultiFab,AMREX_SPACEDIM>& areafrac,
                              const Array<MultiFab,AMREX_SPACEDIM>& facecent,
                              const Array<MultiFab,AMREX_SPACEDIM>& edgecent,
                              const MultiFab& levelset, const Geometry& geom,
                              const IntVect& ngrow, bool extend_domain_face, int max_grid_size) const;
};

}

#endif
