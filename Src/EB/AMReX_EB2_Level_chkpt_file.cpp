#include <AMReX_EB2_Level_chkpt_file.H>

namespace amrex { namespace EB2 {

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, ChkptFile const& chkpt_file,
        Geometry const& geom, int max_grid_size, int ngrow, bool extend_domain_face)
    : GShopLevel<ChkptFile>(is, geom)
{
    BL_PROFILE("EB2::ChkptFileLevel()-fine");

    define_fine_chkpt_file(chkpt_file, geom, max_grid_size, ngrow, extend_domain_face);
}

void ChkptFileLevel::define_fine_chkpt_file(ChkptFile const& chkpt_file,
        Geometry const& geom, int max_grid_size, int ngrow, bool extend_domain_face)
{
    BL_PROFILE("EB2::ChkptFileLevel()-define-fine-chkptfile");

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

   BoxList bl(domain);
   bl.maxSize(max_grid_size);
    if (m_ngrow != 0) {
        const IntVect& domlo = domain.smallEnd();
        const IntVect& domhi = domain.bigEnd();
        for (auto& b : bl) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_ngrow[idim] != 0) {
                    if (b.smallEnd(idim) == domlo[idim]) {
                        b.growLo(idim,m_ngrow[idim]);
                    }
                    if (b.bigEnd(idim) == domhi[idim]) {
                        b.growHi(idim,m_ngrow[idim]);
                    }
                }
            }
        }
    }

    chkpt_file.fill_from_chkpt_file(m_grids, m_dmap, m_volfrac, m_centroid, m_bndryarea,
            m_bndrycent, m_bndrynorm, m_areafrac, m_facecent, m_edgecent, m_levelset, GFab::ng);

    BoxList uncovered_bl = m_grids.boxList();
    BoxList covered_bl;
    for (const auto& bx: bl) {
        if (!uncovered_bl.contains(BoxList(bx))) {
            covered_bl.push_back(bx);
        }
    }

    m_mgf.define(m_grids, m_dmap);
    const int ng = GFab::ng;
    MFInfo mf_info;
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];

        const auto& levelset = m_levelset.const_array(mfi);
        const Box& bxg2 = amrex::grow(gfab.validbox(),2);
        const Box& nodal_box = amrex::surroundingNodes(bxg2);
        const auto& ls = gfab.getLevelSet().array();

        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(nodal_box, i, j, k,
        {
            ls(i,j,k) = levelset(i,j,k);
        });

        auto& cellflag = m_cellflag[mfi];
        gfab.buildTypes(cellflag);
    }

    Vector<Box> cut_boxes;
    for (MFIter mfi(m_grids, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1));

        auto& flagfab = m_cellflag[mfi];
        const auto& flag = m_cellflag[mfi].const_array();
        const auto& vfrac = m_volfrac.const_array(mfi);
        const auto& apx = m_areafrac[0].const_array(mfi);

        if (flagfab.getType(gbx & bounding_box) == FabType::singlevalued) {
            cut_boxes.push_back(vbx);
        }
    }

    AllGatherBoxes(cut_boxes);

    //Print() << "cut_boxes = " << cut_boxes.size() << std::endl;

    if ( cut_boxes.empty() &&
            !covered_bl.isEmpty())
    {
        Abort("AMReX_EB2_Level.H: Domain is completely covered");
    }

    if (!covered_bl.isEmpty()) {
        m_covered_grids = BoxArray(covered_bl);
    }

    if (cut_boxes.empty()) {
        m_allregular = true;
        m_ok = true;
        return;
    }

    EBCellFlagFab cellflagtmp;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];
        const Box& vbx = mfi.validbox();
        const Box& bxg1 = amrex::grow(vbx,1);
        Array4<EBCellFlag> const& cell = m_cellflag.array(mfi);

        cellflagtmp.resize(m_cellflag[mfi].box());
        Elixir cellflagtmp_eli = cellflagtmp.elixir();
        Array4<EBCellFlag> const& ctmp = cellflagtmp.array();

        auto& facetype = gfab.getFaceType();
        AMREX_D_TERM(Array4<Type_t> const& fx = facetype[0].array();,
                     Array4<Type_t> const& fy = facetype[1].array();,
                     Array4<Type_t> const& fz = facetype[2].array(););


        AMREX_D_TERM(Array4<Real const> const& apx = m_areafrac[0].const_array(mfi);,
                     Array4<Real const> const& apy = m_areafrac[1].const_array(mfi);,
                     Array4<Real const> const& apz = m_areafrac[2].const_array(mfi););

        const Box& xbx = amrex::grow(amrex::surroundingNodes(vbx,0),0);
        AMREX_HOST_DEVICE_FOR_3D ( xbx, i, j, k,
        {
            if (apx(i,j,k) == 0.0_rt) {
                fx(i,j,k) = Type::covered;
            } else if (apx(i,j,k) == 1.0_rt) {
                fx(i,j,k) = Type::regular;
            }
        });

        const Box& ybx = amrex::grow(amrex::surroundingNodes(vbx,1),0);
        AMREX_HOST_DEVICE_FOR_3D ( ybx, i, j, k,
        {
            if (apy(i,j,k) == 0.0_rt) {
                fy(i,j,k) = Type::covered;
            } else if (apy(i,j,k) == 1.0_rt) {
                fy(i,j,k) = Type::regular;
            }
        });

        const Box& zbx = amrex::grow(amrex::surroundingNodes(vbx,2),0);
        AMREX_HOST_DEVICE_FOR_3D ( zbx, i, j, k,
        {
            if (apz(i,j,k) == 0.0_rt) {
                fz(i,j,k) = Type::covered;
            } else if (apz(i,j,k) == 1.0_rt) {
                fz(i,j,k) = Type::regular;
            }
        });

        AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
        {
            if (cell(i,j,k).isSingleValued()) {
                if (fx(i,j,k) == Type::covered && fx(i+1,j,k) == Type::covered &&
                    fy(i,j,k) == Type::covered && fy(i,j+1,k) == Type::covered &&
                    fz(i,j,k) == Type::covered && fz(i,j,k+1) == Type::covered)
                {
                    cell(i,j,k).setCovered();
                }
                else if (fx(i,j,k) == Type::regular && fx(i+1,j,k) == Type::regular &&
                         fy(i,j,k) == Type::regular && fy(i,j+1,k) == Type::regular &&
                         fz(i,j,k) == Type::regular && fz(i,j,k+1) == Type::regular)
                {
                    cell(i,j,k).setRegular();
                }
            }
        });


        // Build neighbors.  By default all 26 neighbors are already set.
        AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
        {
            if (cell(i,j,k).isCovered()) {
                cell(i,j,k).setDisconnected();
            } else {
                auto flg = cell(i,j,k);

                if (fx(i,j,k) == Type::covered) {
                    flg.setDisconnected(-1,0,0);
                }
                if (fx(i+1,j,k) == Type::covered) {
                    flg.setDisconnected(1,0,0);
                }
                if (fy(i,j,k) == Type::covered) {
                    flg.setDisconnected(0,-1,0);
                }
                if (fy(i,j+1,k) == Type::covered) {
                    flg.setDisconnected(0,1,0);
                }
                if (fz(i,j,k) == Type::covered) {
                    flg.setDisconnected(0,0,-1);
                }
                if (fz(i,j,k+1) == Type::covered) {
                    flg.setDisconnected(0,0,1);
                }

                // x-y
                if ((fx(i,j,k) == Type::covered || fy(i-1,j,k) == Type::covered) &&
                    (fx(i,j-1,k) == Type::covered || fy(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(-1,-1,0);
                }

                if ((fx(i+1,j,k) == Type::covered || fy(i+1,j,k) == Type::covered) &&
                    (fx(i+1,j-1,k) == Type::covered || fy(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(1,-1,0);
                }

                if ((fx(i,j,k) == Type::covered || fy(i-1,j+1,k) == Type::covered) &&
                    (fx(i,j+1,k) == Type::covered || fy(i,j+1,k) == Type::covered))
                {
                    flg.setDisconnected(-1,1,0);
                }

                if ((fx(i+1,j,k) == Type::covered || fy(i+1,j+1,k) == Type::covered) &&
                    (fx(i+1,j+1,k) == Type::covered || fy(i,j+1,k) == Type::covered))
                {
                    flg.setDisconnected(1,1,0);
                }

                // x-z
                if ((fx(i,j,k) == Type::covered || fz(i-1,j,k) == Type::covered) &&
                    (fx(i,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(-1,0,-1);
                }

                if ((fx(i+1,j,k) == Type::covered || fz(i+1,j,k) == Type::covered) &&
                    (fx(i+1,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(1,0,-1);
                }

                if ((fx(i,j,k) == Type::covered || fz(i-1,j,k+1) == Type::covered) &&
                    (fx(i,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
                {
                    flg.setDisconnected(-1,0,1);
                }

                if ((fx(i+1,j,k) == Type::covered || fz(i+1,j,k+1) == Type::covered) &&
                    (fx(i+1,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
                {
                    flg.setDisconnected(1,0,1);
                }

                // y-z
                if ((fy(i,j,k) == Type::covered || fz(i,j-1,k) == Type::covered) &&
                    (fy(i,j,k-1) == Type::covered || fz(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(0,-1,-1);
                }

                if ((fy(i,j+1,k) == Type::covered || fz(i,j+1,k) == Type::covered) &&
                    (fy(i,j+1,k-1) == Type::covered || fz(i,j,k) == Type::covered))
                {
                    flg.setDisconnected(0,1,-1);
                }

                if ((fy(i,j,k) == Type::covered || fz(i,j-1,k+1) == Type::covered) &&
                    (fy(i,j,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
                {
                    flg.setDisconnected(0,-1,1);
                }

                if ((fy(i,j+1,k) == Type::covered || fz(i,j+1,k+1) == Type::covered) &&
                    (fy(i,j+1,k+1) == Type::covered || fz(i,j,k+1) == Type::covered))
                {
                    flg.setDisconnected(0,1,1);
                }

                cell(i,j,k) = flg;
            }

            ctmp(i,j,k) = cell(i,j,k);
        });

        AMREX_HOST_DEVICE_FOR_3D ( vbx, i, j, k,
        {
            if (!cell(i,j,k).isCovered()) {
                auto tmpflg = ctmp(i,j,k);
                auto newflg = tmpflg;

                // -1, -1, -1 corner
                if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0,-1,-1)) &&
                    (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected(-1, 0,-1)) &&
                    (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected(-1,-1, 0)))
                {
                    newflg.setDisconnected(-1,-1,-1);
                }

                // 1, -1, -1 corner
                if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0,-1,-1)) &&
                    (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected( 1, 0,-1)) &&
                    (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected( 1,-1, 0)))
                {
                    newflg.setDisconnected(1,-1,-1);
                }

                // -1, 1, -1 corner
                if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0, 1,-1)) &&
                    (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected(-1, 0,-1)) &&
                    (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected(-1, 1, 0)))
                {
                    newflg.setDisconnected(-1, 1,-1);
                }

                // 1, 1, -1 corner
                if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0, 1,-1)) &&
                    (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected( 1, 0,-1)) &&
                    (tmpflg.isDisconnected( 0, 0,-1) || ctmp(i  ,j  ,k-1).isDisconnected( 1, 1, 0)))
                {
                    newflg.setDisconnected(1, 1,-1);
                }

                // -1, -1, 1 corner
                if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0,-1, 1)) &&
                    (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected(-1, 0, 1)) &&
                    (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected(-1,-1, 0)))
                {
                    newflg.setDisconnected(-1,-1, 1);
                }

                // 1, -1, 1 corner
                if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0,-1, 1)) &&
                    (tmpflg.isDisconnected( 0,-1, 0) || ctmp(i  ,j-1,k  ).isDisconnected( 1, 0, 1)) &&
                    (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected( 1,-1, 0)))
                {
                    newflg.setDisconnected(1,-1, 1);
                }

                // -1, 1, 1 corner
                if ((tmpflg.isDisconnected(-1, 0, 0) || ctmp(i-1,j  ,k  ).isDisconnected( 0, 1, 1)) &&
                    (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected(-1, 0, 1)) &&
                    (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected(-1, 1, 0)))
                {
                    newflg.setDisconnected(-1,1,1);
                }

                // 1, 1, 1 corner
                if ((tmpflg.isDisconnected( 1, 0, 0) || ctmp(i+1,j  ,k  ).isDisconnected( 0, 1, 1)) &&
                    (tmpflg.isDisconnected( 0, 1, 0) || ctmp(i  ,j+1,k  ).isDisconnected( 1, 0, 1)) &&
                    (tmpflg.isDisconnected( 0, 0, 1) || ctmp(i  ,j  ,k+1).isDisconnected( 1, 1, 0)))
                {
                    newflg.setDisconnected(1,1,1);
                }

                cell(i,j,k) = newflg;
            }
        });
    }

    m_ok = true;
}

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
        const Geometry& geom, ChkptFileLevel& fineLevel)
: GShopLevel<ChkptFile>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}}
