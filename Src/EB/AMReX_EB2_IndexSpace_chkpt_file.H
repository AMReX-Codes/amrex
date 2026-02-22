#ifndef AMREX_EB2_INDEXSPACE_CHKPTFILE_H_
#define AMREX_EB2_INDEXSPACE_CHKPTFILE_H_
#include <AMReX_Config.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level_chkpt_file.H>

#include <string>

/**
 * \file AMReX_EB2_IndexSpace_chkpt_file.H
 *
 * IndexSpace helper for reconstructing EB data from checkpoint files.
 */

namespace amrex::EB2 {

//! IndexSpace that reads EB levels from a checkpoint file on disk.
class IndexSpaceChkptFile
    : public IndexSpace
{
public:

    /**
     * \brief Build an EB hierarchy described by \p chkpt_file.
     *
     * \param chkpt_file     Reader encapsulating the checkpoint contents.
     * \param geom           Finest-level geometry descriptor.
     * \param required_coarsening_level Number of coarse levels the app requires.
     * \param max_coarsening_level      Max coarse level available in the file.
     * \param ngrow          Number of EB grow cells to load on level 0.
     * \param build_coarse_level_by_coarsening Whether to coarsen from fine data.
     * \param extend_domain_face        Extend EB data to domain boundaries.
     */
    IndexSpaceChkptFile (const ChkptFile& chkpt_file,
                         const Geometry& geom, int required_coarsening_level,
                         int max_coarsening_level, int ngrow,
                         bool build_coarse_level_by_coarsening,
                         bool extend_domain_face);

    IndexSpaceChkptFile (IndexSpaceChkptFile const&) = delete;
    IndexSpaceChkptFile (IndexSpaceChkptFile &&) = delete;
    void operator= (IndexSpaceChkptFile const&) = delete;
    void operator= (IndexSpaceChkptFile &&) = delete;

    ~IndexSpaceChkptFile () override = default;

    //! Return the EB level matching \p geom.Domain().
    [[nodiscard]] const Level& getLevel (const Geometry& geom) const final;
    //! Return the Geometry whose domain equals \p dom.
    [[nodiscard]] const Geometry& getGeometry (const Box& dom) const final;
    //! Domain of the coarsest level loaded from the checkpoint.
    [[nodiscard]] const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }
    //! Insert additional fine levels generated from checkpoint data.
    void addFineLevels (int num_new_fine_levels) final;
    //! Add regular coarse levels beneath the existing stack.
    void addRegularCoarseLevels (int num_new_coarse_levels) final;

private:

    Vector<ChkptFileLevel> m_chkpt_file_level;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
};

}

#endif
