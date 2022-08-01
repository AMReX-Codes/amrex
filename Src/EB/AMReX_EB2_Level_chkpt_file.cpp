#include <AMReX_EB2_Level_chkpt_file.H>
#include <AMReX_EB2_C.H>

#include <AMReX_MultiFabUtil.H>

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

    ignore_unused(max_grid_size, extend_domain_face);

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

    const int ng = GFab::ng;
    chkpt_file.read_from_chkpt_file(m_grids, m_covered_grids,
            m_dmap, m_volfrac, m_centroid, m_bndryarea,
            m_bndrycent, m_bndrynorm, m_areafrac, m_facecent,
            m_edgecent, m_levelset, ng, geom);


    if ( m_grids.empty() &&
            !m_covered_grids.empty())
    {
        Abort("AMReX_EB2_Level.H: Domain is completely covered");
    }

    if (m_grids.empty()) {
        m_allregular = true;
        m_ok = true;
        return;
    }

    set_invalid_ghost_data_covered();
    set_invalid_ghost_data_extended();

    m_mgf.define(m_grids, m_dmap);
    MFInfo mf_info;
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];

        const auto& levelset = m_levelset.const_array(mfi);
        const Box& bxg2 = amrex::grow(gfab.validbox(),ng);
        const Box& nodal_box = amrex::surroundingNodes(bxg2);
        const auto& ls = gfab.getLevelSet().array();

        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(nodal_box, i, j, k,
        {
            ls(i,j,k) = levelset(i,j,k);
        });

        auto& cellflag = m_cellflag[mfi];
        gfab.buildTypes(cellflag);
    }

    finalize_cell_flags();
}

void ChkptFileLevel::finalize_cell_flags () {

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

        const Box& xbx = amrex::grow(amrex::surroundingNodes(vbx,0),1);
        AMREX_HOST_DEVICE_FOR_3D ( xbx, i, j, k,
        {
            if (apx(i,j,k) == 0.0_rt) {
                fx(i,j,k) = Type::covered;
            } else if (apx(i,j,k) == 1.0_rt) {
                fx(i,j,k) = Type::regular;
            }
        });

        const Box& ybx = amrex::grow(amrex::surroundingNodes(vbx,1),1);
        AMREX_HOST_DEVICE_FOR_3D ( ybx, i, j, k,
        {
            if (apy(i,j,k) == 0.0_rt) {
                fy(i,j,k) = Type::covered;
            } else if (apy(i,j,k) == 1.0_rt) {
                fy(i,j,k) = Type::regular;
            }
        });

#if (AMREX_SPACEDIM == 3)
        const Box& zbx = amrex::grow(amrex::surroundingNodes(vbx,2),1);
        AMREX_HOST_DEVICE_FOR_3D ( zbx, i, j, k,
        {
            if (apz(i,j,k) == 0.0_rt) {
                fz(i,j,k) = Type::covered;
            } else if (apz(i,j,k) == 1.0_rt) {
                fz(i,j,k) = Type::regular;
            }
        });
#endif


#if (AMREX_SPACEDIM == 2)
        ignore_unused(ctmp);
        AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
        {
            ignore_unused(k);
            if (cell(i,j,0).isSingleValued()) {
                if (fx(i,j,0) == Type::regular && fx(i+1,j,0) == Type::regular &&
                    fy(i,j,0) == Type::regular && fy(i,j+1,0) == Type::regular)
                {
                    cell(i,j,0).setRegular();
                }
                else if (fx(i,j,0) == Type::covered && fx(i+1,j,0) == Type::covered &&
                         fy(i,j,0) == Type::covered && fy(i,j+1,0) == Type::covered)
                {
                    cell(i,j,0).setCovered();
                }
            }
        });

        set_connection_flags(bxg1, cell, fx, fy);

#else
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

        set_connection_flags(vbx, bxg1, cell, ctmp, fx, fy, fz);

#endif

    }

    m_ok = true;
}

void ChkptFileLevel::set_invalid_ghost_data_covered ()
{
    const std::vector<IntVect>& pshifts = m_geom.periodicity().shiftIntVect();

    if (!m_covered_grids.empty()) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            std::vector<std::pair<int,Box> > isects;
            for (MFIter mfi(m_levelset); mfi.isValid(); ++mfi)
            {
                const auto& lsfab = m_levelset.array(mfi);
                const Box& ccbx = amrex::enclosedCells(mfi.fabbox());
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(ccbx+iv, isects);
                    for (const auto& is : isects) {
                        const Box& fbx = amrex::surroundingNodes(is.second-iv);
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(fbx, i, j, k,
                        {
                            lsfab(i,j,k) = 1.0;
                        });
                    }
                }
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            std::vector<std::pair<int,Box> > isects;
            for (MFIter mfi(m_volfrac); mfi.isValid(); ++mfi)
            {
                auto const& fab = m_volfrac.array(mfi);
                const Box& bx = mfi.fabbox();
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(bx+iv, isects);
                    for (const auto& is : isects) {
                        Box const& ibox = is.second-iv;
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(ibox, i, j, k,
                        {
                            fab(i,j,k) = 0.0;  // covered cells
                        });
                    }
                }
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            std::vector<std::pair<int,Box> > isects;
            for (MFIter mfi(m_areafrac[0]); mfi.isValid(); ++mfi)
            {
                const Box& ccbx = amrex::enclosedCells((m_areafrac[0])[mfi].box());
                AMREX_D_TERM(auto const& apx = m_areafrac[0].array(mfi);,
                             auto const& apy = m_areafrac[1].array(mfi);,
                             auto const& apz = m_areafrac[2].array(mfi););
                for (const auto& iv : pshifts)
                {
                    m_covered_grids.intersections(ccbx+iv, isects);
                    for (const auto& is : isects) {
                        if (Gpu::inLaunchRegion()) {
                            const Box& bx = is.second-iv;
                            AMREX_D_TERM(const Box& xbx = amrex::surroundingNodes(bx,0);,
                                         const Box& ybx = amrex::surroundingNodes(bx,1);,
                                         const Box& zbx = amrex::surroundingNodes(bx,2););
                            amrex::ParallelFor(AMREX_D_DECL(xbx,ybx,zbx),
                              [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apx(i,j,k) = 0.0;
                              }
#if (AMREX_SPACEDIM >= 2)
                            , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apy(i,j,k) = 0.0;
                              }
#if (AMREX_SPACEDIM == 3)
                            , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                              {
                                  apz(i,j,k) = 0.0;
                              }
#endif
#endif
                            );
                        } else {
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                const Box& fbx = amrex::surroundingNodes(is.second-iv,idim);
                                (m_areafrac[idim])[mfi].setVal<RunOn::Host>(0.0, fbx, 0, 1);
                            }
                        }
                    }
                }
            }
        }

        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                auto& edgecent = m_edgecent[idim];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                {
                    BoxArray const& covered_edge_grids = amrex::convert(m_covered_grids,
                                                                        edgecent.ixType());
                    std::vector<std::pair<int,Box> > isects;
                    for (MFIter mfi(edgecent); mfi.isValid(); ++mfi)
                    {
                        auto const& fab = edgecent.array(mfi);
                        const Box& bx = mfi.fabbox();
                        for (const auto& iv : pshifts)
                        {
                            covered_edge_grids.intersections(bx+iv, isects);
                            for (const auto& is : isects) {
                                Box const& ibox = is.second-iv;
                                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(ibox, i, j, k,
                                {
                                    fab(i,j,k) = Real(-1.0);  // covered edges
                                });
                            }
                        }
                    }
                }
            }
        }
    }
}

void ChkptFileLevel::set_invalid_ghost_data_extended ()
{
    Box domain_grown = m_geom.Domain();
    domain_grown.grow(m_ngrow);
    const Box& bounding_box = domain_grown;

    const std::vector<IntVect>& pshifts = m_geom.periodicity().shiftIntVect();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(m_levelset); mfi.isValid(); ++mfi)
        {
            const auto& lsfab = m_levelset.array(mfi);
            const Box& ccbx = amrex::enclosedCells(mfi.fabbox());
            Box bbx = amrex::enclosedCells(bounding_box);
            for (const auto& iv : pshifts)
            {
               if (!bbx.contains(ccbx+iv)) {
                 bbx.surroundingNodes();
                 const Box& fbx = amrex::surroundingNodes(ccbx);

                 const auto blo = amrex::lbound(bbx);
                 const auto bhi = amrex::ubound(bbx);

                 AMREX_HOST_DEVICE_PARALLEL_FOR_3D(fbx, i, j, k,
                 {
                     if (i < blo.x || i > bhi.x || j < blo.y ||
                         i > bhi.y || k < blo.z || k > bhi.z) {
                         lsfab(i,j,k) = lsfab(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                           amrex::Clamp(j,blo.y,bhi.y),
                                                           amrex::Clamp(k,blo.z,bhi.z)));
                     }
                 });
               }
            }
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(m_volfrac); mfi.isValid(); ++mfi)
        {
            auto const& fab = m_volfrac.array(mfi);
            const Box& bx = mfi.fabbox();
            for (const auto& iv : pshifts)
            {
               if (!bounding_box.contains(bx+iv)) {
                 const auto blo = amrex::lbound(bounding_box);
                 const auto bhi = amrex::ubound(bounding_box);

                 AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                 {
                     if (i < blo.x || i > bhi.x || j < blo.y ||
                         i > bhi.y || k < blo.z || k > bhi.z) {
                         fab(i,j,k) = fab(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                       amrex::Clamp(j,blo.y,bhi.y),
                                                       amrex::Clamp(k,blo.z,bhi.z)));
                     }
                 });
               }
            }
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(m_areafrac[0]); mfi.isValid(); ++mfi)
        {
            const Box& ccbx = amrex::enclosedCells((m_areafrac[0])[mfi].box());
            const Box& bbx = amrex::enclosedCells(bounding_box);
            AMREX_D_TERM(auto const& apx = m_areafrac[0].array(mfi);,
                         auto const& apy = m_areafrac[1].array(mfi);,
                         auto const& apz = m_areafrac[2].array(mfi););
            for (const auto& iv : pshifts)
            {
                if (!bbx.contains(ccbx+iv)) {
                    const auto blo = amrex::lbound(bbx);
                    const auto bhi = amrex::ubound(bbx);

                    AMREX_D_TERM(const Box& xbx = amrex::surroundingNodes(ccbx+iv,0);,
                                 const Box& ybx = amrex::surroundingNodes(ccbx+iv,1);,
                                 const Box& zbx = amrex::surroundingNodes(ccbx+iv,2););
                    amrex::ParallelFor(AMREX_D_DECL(xbx,ybx,zbx),
                      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                          if (i < blo.x || i > bhi.x || j < blo.y ||
                              i > bhi.y || k < blo.z || k > bhi.z) {
                              apx(i,j,k) = apx(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                            amrex::Clamp(j,blo.y,bhi.y),
                                                            amrex::Clamp(k,blo.z,bhi.z)));
                          }
                      }
#if (AMREX_SPACEDIM >= 2)
                    , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                          if (i < blo.x || i > bhi.x || j < blo.y ||
                              i > bhi.y || k < blo.z || k > bhi.z) {
                              apy(i,j,k) = apy(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                            amrex::Clamp(j,blo.y,bhi.y),
                                                            amrex::Clamp(k,blo.z,bhi.z)));
                          }
                      }
#if (AMREX_SPACEDIM == 3)
                    , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                          if (i < blo.x || i > bhi.x || j < blo.y ||
                              i > bhi.y || k < blo.z || k > bhi.z) {
                              apz(i,j,k) = apz(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                            amrex::Clamp(j,blo.y,bhi.y),
                                                            amrex::Clamp(k,blo.z,bhi.z)));
                          }
                      }
#endif
#endif
                    );
                }
            }
        }
    }

    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            auto& edgecent = m_edgecent[idim];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                Box const& bbx = amrex::convert(bounding_box, edgecent.ixType());
                for (MFIter mfi(edgecent); mfi.isValid(); ++mfi)
                {
                    auto const& fab = edgecent.array(mfi);
                    const Box& bx = mfi.fabbox();
                    for (const auto& iv : pshifts)
                    {
                        if (!bbx.contains(bx+iv)) {
                            const auto blo = amrex::lbound(bbx);
                            const auto bhi = amrex::ubound(bbx);

                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                            {
                                if (i < blo.x || i > bhi.x || j < blo.y ||
                                    i > bhi.y || k < blo.z || k > bhi.z) {
                                 fab(i,j,k) = fab(AMREX_D_DECL(amrex::Clamp(i,blo.x,bhi.x),
                                                               amrex::Clamp(j,blo.y,bhi.y),
                                                               amrex::Clamp(k,blo.z,bhi.z)));
                                }
                            });
                        }
                    }
                }
            }
        }
    }
}

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
        const Geometry& geom, ChkptFileLevel& fineLevel)
: GShopLevel<ChkptFile>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}}
