#ifndef AMREX_EB2_INDEXSPACE_CHKPTFILE_H_
#define AMREX_EB2_INDEXSPACE_CHKPTFILE_H_
#include <AMReX_Config.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level_chkpt_file.H>

#include <string>

namespace amrex::EB2 {

class IndexSpaceChkptFile
    : public IndexSpace
{
public:

    IndexSpaceChkptFile (const ChkptFile& chkpt_file,
                         const Geometry& geom, int required_coarsening_level,
                         int max_coarsening_level, int ngrow,
                         bool build_coarse_level_by_coarsening,
                         bool extend_domain_face);

    IndexSpaceChkptFile (IndexSpaceChkptFile const&) = delete;
    IndexSpaceChkptFile (IndexSpaceChkptFile &&) = delete;
    void operator= (IndexSpaceChkptFile const&) = delete;
    void operator= (IndexSpaceChkptFile &&) = delete;

    virtual ~IndexSpaceChkptFile () {}

    virtual const Level& getLevel (const Geometry& geom) const final;
    virtual const Geometry& getGeometry (const Box& dom) const final;
    virtual const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }
    virtual void addFineLevels (int num_new_fine_levels) final;

private:

    Vector<ChkptFileLevel> m_chkpt_file_level;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
};

}

#endif
