#include <AMReX_EB2_Level_chkpt_file.H>
#include <AMReX_EB2_C.H>

#include <AMReX_MultiFabUtil.H>

namespace amrex::EB2 {

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, ChkptFile const& chkpt_file,
                                Geometry const& geom, int max_grid_size, int ngrow, bool extend_domain_face)
    : GShopLevel<ChkptFile>(is, geom)
{
    BL_PROFILE("EB2::ChkptFileLevel()-fine");

    define_fine_chkpt_file(chkpt_file, geom, max_grid_size, ngrow, extend_domain_face);
}

void
ChkptFileLevel::define_fine_chkpt_file (ChkptFile const& chkpt_file,
                                        Geometry const& geom, int max_grid_size,
                                        int ngrow, bool extend_domain_face)
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

    const int ng = GFab::ng;
    chkpt_file.read_from_chkpt_file(m_grids, m_covered_grids,
                                    m_dmap, m_volfrac, m_centroid, m_bndryarea,
                                    m_bndrycent, m_bndrynorm, m_areafrac, m_facecent,
                                    m_edgecent, m_levelset, ng, geom, m_ngrow,
                                    extend_domain_face, max_grid_size);


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

void
ChkptFileLevel::finalize_cell_flags ()
{

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        EBCellFlagFab cellflagtmp;
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
}

ChkptFileLevel::ChkptFileLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                                const Geometry& geom, ChkptFileLevel& fineLevel)
: GShopLevel<ChkptFile>(is, ilev, max_grid_size, ngrow, geom, fineLevel)
{}

}
