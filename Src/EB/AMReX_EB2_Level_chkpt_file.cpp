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
    EBCellFlagFab cellflagtmp;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:nsmallcells,nmulticuts)
#endif
    {
        for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
        {
            auto& gfab = m_mgf[mfi];

            const auto& levelset = m_levelset.const_array(mfi);
            const Box& bx = gfab.getLevelSet().box();
            const auto& ls = gfab.getLevelSet().array();

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:nsmallcells,nmulticuts)
#endif
    {
        for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
        {
            const auto& gfab = m_mgf[mfi];
            const Box& vbx = mfi.validbox();
            const Box& bxg1 = amrex::grow(vbx,1);
            Array4<EBCellFlag> const& cell = m_cellflag.array(mfi);

            // set cells in the extended region to covered if the
            // corresponding cell on the domain face is covered
            if(extend_domain_face) {

               Box gdomain = geom.Domain();
               for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                   if (geom.isPeriodic(idim)) {
                       gdomain.setSmall(idim, std::min(gdomain.smallEnd(idim), bxg1.smallEnd(idim)));
                       gdomain.setBig(idim, std::max(gdomain.bigEnd(idim), bxg1.bigEnd(idim)));
                   }
               }

               if (! gdomain.contains(bxg1)) {
                  AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
                  {
                      const auto & dlo = gdomain.loVect();
                      const auto & dhi = gdomain.hiVect();

                      // find the cell(ii,jj,kk) on the corr. domain face
                      // this would have already been set to correct value
                      bool in_extended_domain = false;
                      int ii = i;
                      int jj = j;
                      int kk = k;
                      if(i < dlo[0]) {
                          in_extended_domain = true;
                          ii = dlo[0];
                      }
                      else if(i > dhi[0]) {
                          in_extended_domain = true;
                          ii = dhi[0];
                      }

                      if(j < dlo[1]) {
                          in_extended_domain = true;
                          jj = dlo[1];
                      }
                      else if(j > dhi[1]) {
                          in_extended_domain = true;
                          jj = dhi[1];
                      }

                      if(k < dlo[2]) {
                          in_extended_domain = true;
                          kk = dlo[2];
                      }
                      else if(k > dhi[2]) {
                          in_extended_domain = true;
                          kk = dhi[2];
                      }

                      // set cell in extendable region to covered if necessary
                      if( in_extended_domain && (! cell(i,j,k).isCovered())
                          && cell(ii,jj,kk).isCovered() )
                      {
                          cell(i,j,k).setCovered();
                      }
                  });
               }
            }

            cellflagtmp.resize(m_cellflag[mfi].box());
            Array4<EBCellFlag> const& ctmp = cellflagtmp.array();

            const auto& facetype = gfab.getFaceType();
            AMREX_D_TERM(Array4<Type_t const> const& fx = facetype[0].const_array();,
                         Array4<Type_t const> const& fy = facetype[1].const_array();,
                         Array4<Type_t const> const& fz = facetype[2].const_array(););


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
    }

    m_ok = true;
}

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
        const Geometry& geom, ChkptFileLevel& fineLevel)
: GShopLevel<ChkptFile>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}}
