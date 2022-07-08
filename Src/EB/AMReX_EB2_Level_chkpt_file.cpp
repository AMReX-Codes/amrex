#include <AMReX_EB2_Level_chkpt_file.H>

namespace amrex { namespace EB2 {

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, ChkptFile const& chkpt_file,
        const Geometry& geom, int max_grid_size, int ngrow, bool extend_domain_face)
    : GShopLevel<ChkptFile>(is, geom)
{
    BL_PROFILE("EB2::ChkptFileLevel()-fine");

    m_ngrow = IntVect{static_cast<int>(std::ceil(ngrow/16.)) * 16};

    Box const& domain = geom.Domain();
    Box domain_grown = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            m_ngrow[idim] = 0;
        } else {
            m_ngrow[idim] = std::min(m_ngrow[idim], domain_grown.length(idim));
        }
    }
    domain_grown.grow(m_ngrow);
    Box bounding_box = (extend_domain_face) ? domain : domain_grown;
    bounding_box.surroundingNodes();

    chkpt_file.fill_from_chkpt_file(m_grids, m_dmap, m_volfrac, m_centroid, m_bndryarea,
            m_bndrycent, m_bndrynorm, m_areafrac, m_facecent, m_edgecent, m_levelset, GFab::ng);

    m_volfrac.FillBoundary(m_geom.periodicity());
    m_centroid.FillBoundary(m_geom.periodicity());
    m_bndryarea.FillBoundary(m_geom.periodicity());
    m_bndrycent.FillBoundary(m_geom.periodicity());
    m_bndrynorm.FillBoundary(m_geom.periodicity());

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].FillBoundary(m_geom.periodicity());
        m_facecent[idim].FillBoundary(m_geom.periodicity());
        m_edgecent[idim].FillBoundary(m_geom.periodicity());
    }

    m_levelset.FillBoundary(m_geom.periodicity());

    m_mgf.define(m_grids, m_dmap);
    const int ng = GFab::ng;
    MFInfo mf_info;
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:nsmallcells,nmulticuts)
#endif
    {
        for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
        {
            auto& gfab = m_mgf[mfi];
            const Box& vbx = gfab.validbox();

            const auto& levelset = m_levelset.const_array(mfi);
            const auto& ls = gfab.getLevelSet().array();

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(vbx, i, j, k,
            {
                ls(i,j,k) = levelset(i,j,k);
            });

            auto& cellflag = m_cellflag[mfi];
            gfab.buildTypes(cellflag);
        }
    }


    Vector<Box> cut_boxes;
    Vector<Box> covered_boxes;
    for (MFIter mfi(m_grids, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1));

        auto& flagfab = m_cellflag[mfi];
        const auto& flag = m_cellflag[mfi].const_array();
        const auto& vfrac = m_volfrac.const_array(mfi);
        const auto& apx = m_areafrac[0].const_array(mfi);

        if (flagfab.getType(gbx & bounding_box) == FabType::covered) {
            covered_boxes.push_back(vbx);
        } else if (flagfab.getType(gbx & bounding_box) == FabType::singlevalued) {
            cut_boxes.push_back(vbx);
        }
    }

    AllGatherBoxes(cut_boxes);
    AllGatherBoxes(covered_boxes);

    //Print() << "cut_boxes = " << cut_boxes.size() << std::endl;
    //Print() << "covered_boxes = " << covered_boxes.size() << std::endl;

    if ( cut_boxes.empty() &&
            !covered_boxes.empty())
    {
        Abort("AMReX_EB2_Level.H: Domain is completely covered");
    }

    if (!covered_boxes.empty()) {
        m_covered_grids = BoxArray(BoxList(std::move(covered_boxes)));
    }

    if (cut_boxes.empty()) {
        m_grids = BoxArray();
        m_dmap = DistributionMapping();
        m_allregular = true;
        m_ok = true;
        return;
    }

    m_grids = BoxArray(BoxList(std::move(cut_boxes)));
    m_dmap = DistributionMapping(m_grids);

    m_ok = true;
}

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
        const Geometry& geom, ChkptFileLevel& fineLevel)
: GShopLevel<ChkptFile>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}}
