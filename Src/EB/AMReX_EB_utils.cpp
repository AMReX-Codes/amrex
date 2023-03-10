#include <AMReX_MultiFab.H>
#include <AMReX_EB_utils.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_REAL.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

namespace amrex {

#if (AMREX_SPACEDIM > 1)
    //
    // Do small cell redistribution on one FAB
    //
    void apply_eb_redistribution ( const Box& bx,
                                   MultiFab& div_mf,
                                   MultiFab& divc_mf,
                                   const MultiFab& weights,
                                   MFIter* mfi,
                                   const int icomp,
                                   const int ncomp,
                                   const EBCellFlagFab& flags_fab,
                                   const MultiFab* volfrac,
                                   Box& /*domain*/,
                                   const Geometry & geom)
    {
        //
        // Check that grid is uniform
        //
        const Real* dx = geom.CellSize();

#if (AMREX_SPACEDIM == 2)
        if (! amrex::almostEqual(dx[0], dx[1]))
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#elif (AMREX_SPACEDIM == 3)
        if( ! amrex::almostEqual(dx[0],dx[1]) ||
            ! amrex::almostEqual(dx[0],dx[2]) ||
            ! amrex::almostEqual(dx[1],dx[2]) )
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#endif

        //
        // Get array4 from arguments
        //
        Array4<Real> const& div  = div_mf.array(*mfi);
        Array4<Real> const& divc = divc_mf.array(*mfi);
        auto const&         wt   = weights.array(*mfi);
        auto const&        flags = flags_fab.array();
        auto const&        vfrac = volfrac->array(*mfi);

        apply_flux_redistribution ( bx, div, divc, wt, icomp, ncomp, flags, vfrac, geom);
    }

    //
    // Do small cell redistribution on one FAB with the Array4's already passed in
    //
    void apply_flux_redistribution ( const Box& bx,
                                     Array4<Real      > const& div,
                                     Array4<Real const> const& divc,
                                     Array4<Real const> const& wt,
                                     const int icomp,
                                     const int ncomp,
                                     Array4<EBCellFlag const> const& flags,
                                     Array4<Real const>    const& vfrac,
                                     const Geometry & geom)
    {
        //
        // Check that grid is uniform
        //
        const Real* dx = geom.CellSize();

#if (AMREX_SPACEDIM == 2)
        if (! amrex::almostEqual(dx[0], dx[1]))
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#elif (AMREX_SPACEDIM == 3)
        if( ! amrex::almostEqual(dx[0],dx[1]) ||
            ! amrex::almostEqual(dx[0],dx[2]) ||
            ! amrex::almostEqual(dx[1],dx[2]) )
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#endif

        const Box dbox = geom.growPeriodicDomain(2);

        //
        // Get array from arguments
        //
        // Array4<Real> const& div  = div_mf.array(*mfi);
        // Array4<Real> const& divc = divc_mf.array(*mfi);
        // auto const&         wt   = weights.array(*mfi);
        // auto const&        flags = flags_fab.array();
        // auto const&        vfrac = volfrac->array(*mfi);

        const Box& grown1_bx = amrex::grow(bx,1);
        const Box& grown2_bx = amrex::grow(bx,2);

        //
        // Working arrays
        //
        FArrayBox  delm_fab(grown1_bx,ncomp);
        FArrayBox  optmp_fab(grown2_bx,ncomp);
        FArrayBox  mask_fab(grown2_bx);

        Array4<Real> const& optmp = optmp_fab.array();
        Array4<Real> const& mask  = mask_fab.array();
        Array4<Real> const& delm  = delm_fab.array();

        //
        // Array "mask" is used to sever the link to ghost cells when the BCs
        // are not periodic
        // It is set to 1 when a cell can be used in computations, 0 otherwise
        //
        AMREX_FOR_3D(grown2_bx, i, j, k,
        {
            mask(i,j,k) = (dbox.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
        });

        //
        // Init to zero tmp array
        //
        AMREX_FOR_4D(grown2_bx, ncomp, i, j, k, n,
        {
            optmp(i,j,k,n) = 0;
        });

        //
        // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
        //
        AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
        {
        if(flags(i,j,k).isSingleValued())
        {
            Real divnc(0.0);
            Real vtot(0.0);
            Real wted_frac(0.0);
            int  ks = (AMREX_SPACEDIM == 3) ? -1 : 0;
            int  ke = (AMREX_SPACEDIM == 3) ?  1 : 0;

            for (int kk(ks); kk <= ke; ++kk) {
                for (int jj(-1); jj <= 1; ++jj) {
                    for (int ii(-1); ii <= 1; ++ii) {
                        if( (ii != 0 || jj != 0 || kk != 0) &&
                            flags(i,j,k).isConnected(ii,jj,kk) &&
                            dbox.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                        {
                            wted_frac = vfrac(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                            vtot   += wted_frac;
                            divnc  += wted_frac * divc(i+ii,j+jj,k+kk,n);
                        }
                    }
                }
            }
            divnc /=  (vtot + 1.e-80);

            // We need to multiply by mask to make sure optmp is zero for cells
            // outside the domain for non-cyclic BCs
            optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n)) * mask(i,j,k);
            delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);

        }
        else
        {
            delm(i,j,k,n) = 0;
        }
        });


        //
        // Step 3: redistribute excess/loss of mass
        //
        AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
        {
        if(flags(i,j,k).isSingleValued())
        {
            Real wtot(0.0);
            int  ks = (AMREX_SPACEDIM == 3) ? -1 : 0;
            int  ke = (AMREX_SPACEDIM == 3) ?  1 : 0;

            for (int kk(ks); kk <= ke; ++kk) {
              for (int jj(-1); jj <= 1; ++jj) {
                for (int ii(-1); ii <= 1; ++ii) {

                        if( (ii != 0 || jj != 0 || kk != 0) &&
                            (flags(i,j,k).isConnected(ii,jj,kk)) )
                        {
                            wtot += wt(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                        }

            }}}

            wtot = 1.0/(wtot + 1.e-80);

            for (int kk(ks); kk <= ke; ++kk) {
              for (int jj(-1); jj <= 1; ++jj) {
                for (int ii(-1); ii <= 1; ++ii) {

                        if( (ii != 0 || jj != 0 || kk != 0) &&
                            (flags(i,j,k).isConnected(ii,jj,kk)) &&
                            bx.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))) )
                        {
                            Gpu::Atomic::AddNoRet(&optmp(i+ii,j+jj,k+kk,n),
                                             delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk));
                        }
                }}}

        }
        });

        //
        // Resume the correct sign, AKA return the negative
        //
        AMREX_FOR_4D(bx, ncomp, i, j, k, n,
        {
            div(i,j,k,icomp+n) = divc(i,j,k,n) + optmp(i,j,k,n);
        });

        Gpu::streamSynchronize();
    }

    //
    // Do small cell redistribution on a MultiFab -- with a weighting function
    //
    void single_level_weighted_redistribute (MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                             int div_comp, int ncomp, const Geometry& geom)
    {
        Box domain(geom.Domain());

#ifdef AMREX_USE_FLOAT
        Real covered_val = Real(1.e20);
#else
        Real covered_val = 1.e40;
#endif

        int nghost = 2;
        AMREX_ASSERT(div_tmp_in.nGrow() >= nghost);

        EB_set_covered(div_tmp_in, 0, ncomp, div_tmp_in.nGrow(), covered_val);

        div_tmp_in.FillBoundary(geom.periodicity());

        // Here we take care of both the regular and covered cases ... all we do below is the cut cell cases
        MultiFab::Copy(div_out, div_tmp_in, 0, div_comp, ncomp, 0);

        // Get EB geometric info
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(div_out.Factory());
        const MultiFab* volfrac   = &(ebfactory. getVolFrac());

        for (MFIter mfi(div_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            const Box& bx = mfi.tilebox ();

            // this is to check efficiently if this tile contains any eb stuff
            const auto&  div_fab = static_cast<EBFArrayBox const&>(div_out[mfi]);
            const EBCellFlagFab&  flags = div_fab.getEBCellFlagFab();

            if ( !(flags.getType(amrex::grow(bx,     0)) == FabType::covered) &&
                 !(flags.getType(amrex::grow(bx,nghost)) == FabType::regular) )
            {
                // Compute div(tau) with EB algorithm
                apply_eb_redistribution(bx, div_out, div_tmp_in, weights, &mfi,
                                               div_comp, ncomp, flags, volfrac, domain,
                                               geom);

            }
        }
    }

    //
    // Do small cell redistribution on a MultiFab -- without a weighting function
    //
    void single_level_redistribute (MultiFab& div_tmp_in, MultiFab& div_out,
                                    int div_comp, int ncomp, const Geometry& geom)
    {
        // We create a weighting array to use inside the redistribution array
        MultiFab weights(div_out.boxArray(), div_out.DistributionMap(), 1, div_tmp_in.nGrow());
        weights.setVal(1.0);

        single_level_weighted_redistribute (div_tmp_in, div_out, weights, div_comp, ncomp, geom);
    }
#endif

void FillSignedDistance (MultiFab& mf, bool fluid_has_positive_sign)
{
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(&(mf.Factory()));
    if (factory) {
        FillSignedDistance(mf, *(factory->getEBLevel()), *factory, 1, fluid_has_positive_sign);
    } else {
        mf.setVal(std::numeric_limits<Real>::max());
    }
}

namespace detail
{
// Purpose: Given a collision between particle and EB surface, and
// given that a neighbour cell owns the EB surface, a collision between
// the particle and the EDGE of the EB facet might occur. This
// function returns the coordinates of the closest point on the edge of
// an EB facet. This function does not check of collisions.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
RealVect
facets_nearest_pt (IntVect const& ind_pt, IntVect const& ind_loop, RealVect const& r_vec,
                   RealVect const& eb_normal, RealVect const& eb_p0,
                   GpuArray<Real,AMREX_SPACEDIM> const& dx)
{
    // Enumerate the possible EB facet edges invovlved.
    int n_facets = 0;
    IntVect ind_facets {AMREX_D_DECL(0, 0, 0)};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ( ind_pt[d] != ind_loop[d] ) {
            ind_facets[n_facets++] = d;
        }
    }

    // scalar characterizing EB facet position
    Real eb_h = eb_normal.dotProduct(eb_p0);
    Real min_dist = std::numeric_limits<Real>::max();
    RealVect c_vec;

    // iterate over EB facet edges and find whichever has the closest nearest point
    for (int i_facet=0; i_facet<n_facets; ++i_facet)
    {
        int tmp_facet = ind_facets[i_facet];

        // determine the normal of the cell's facet (cube faces)
        RealVect facet_normal {AMREX_D_DECL(0._rt, 0._rt, 0._rt)};
        facet_normal[tmp_facet] = 1.; // whether facing inwards or outwards is not important here

        // skip cases where cell faces coincide with the eb facets
        if (AMREX_D_TERM(std::abs(eb_normal[0]) == std::abs(facet_normal[0]),
                      && std::abs(eb_normal[1]) == std::abs(facet_normal[1]),
                      && std::abs(eb_normal[2]) == std::abs(facet_normal[2])))
        { continue; }

        int ind_cell = ind_loop[tmp_facet];
        int ind_nb   = ind_pt[tmp_facet];

        // determine position of the cell's facet
        Real f_c;
        if (ind_cell < ind_nb) {
            f_c = static_cast<Real>( ind_cell + 1 ) * dx[tmp_facet];
        } else {
            f_c = static_cast<Real>( ind_cell     ) * dx[tmp_facet];
        }

        RealVect facet_p0{AMREX_D_DECL((static_cast<Real>(ind_loop[0]) + 0.5_rt) * dx[0],
                                       (static_cast<Real>(ind_loop[1]) + 0.5_rt) * dx[1],
                                       (static_cast<Real>(ind_loop[2]) + 0.5_rt) * dx[2])};
        facet_p0[tmp_facet] = f_c;

        // scalar characterizing cell facet position
        Real facet_h = facet_normal.dotProduct(facet_p0);

        // compute EB facet edge by finding the intercept between EB surface (first plane)
        // and the cell's facet (second plane)
        //
        //  Purpose: Calculates the line (represented by a position and a
        //  direction vector) given by the intersection of two planes (defined
        //  by two normal (n1, n2) and two positions (h1 = n1.p1, h2 = n2.p2).
        //
        //  When one plane is the EB surface, and the other is a face of the
        //  cell. Then this line represents the edge of the EB facet.
        //
        Real c_dp = eb_normal.dotProduct(facet_normal);
        Real c_norm = 1._rt - c_dp*c_dp;
        //
        Real c1 = ( eb_h - facet_h * c_dp ) / c_norm;
        Real c2 = ( facet_h - eb_h * c_dp ) / c_norm;
        //
        RealVect edge_p0{AMREX_D_DECL(c1*eb_normal[0] + c2*facet_normal[0],
                                      c1*eb_normal[1] + c2*facet_normal[1],
                                      c1*eb_normal[2] + c2*facet_normal[2])};
#if (AMREX_SPACEDIM == 3)
        RealVect edge_v = eb_normal.crossProduct(facet_normal);

        // this solution is a line representing the closest EB edge, now compute the point
        // on the line which minimizes the distance to the particle
        //
        // Purpose: Given an a line an a point, this finds the point
        // one the line which minimizes the cartesian distance. It also finds
        // the corresponding distance along the line corresponding to this point
        //
        RealVect c = edge_p0 - r_vec;
        Real lambda_tmp = - edge_v.dotProduct(c) / edge_v.dotProduct(edge_v);
        RealVect c_vec_tmp{AMREX_D_DECL(edge_p0[0] + lambda_tmp*edge_v[0],
                                        edge_p0[1] + lambda_tmp*edge_v[1],
                                        edge_p0[2] + lambda_tmp*edge_v[2])};

        // IMPORTANT: this point might be outside the cell
        //  -> in that case, it will be one of the cell's corners
        //
        // if closest point is outside cell, determine the furthest we can go along the
        // EB edge line whilst staying within the cell.
        //
        // Purpose: Given a line which passes through a box in three dimensions
        // (it can pass through the edges). Let lambda be a real value
        // representing the coordinate along the line. This finds
        // the min/max values of lambda, in order for the point described by
        // lambda to be contained within the box.
        //
        Real cx_lo = -std::numeric_limits<Real>::max();
        Real cy_lo = -std::numeric_limits<Real>::max();
        Real cz_lo = -std::numeric_limits<Real>::max();
        Real cx_hi = std::numeric_limits<Real>::max();
        Real cy_hi = std::numeric_limits<Real>::max();
        Real cz_hi = std::numeric_limits<Real>::max();
        Real eps = std::numeric_limits<Real>::epsilon();
        // if the line runs parallel to any of these dimensions (which is true for
        // EB edges), then skip -> the min/max functions at the end will skip them
        // due to the +/-huge(c...) defaults (above).
        if ( std::abs(edge_v[0]) > eps ) {
            cx_lo = -( edge_p0[0] - static_cast<Real>( ind_loop[0]     ) * dx[0] ) / edge_v[0];
            cx_hi = -( edge_p0[0] - static_cast<Real>( ind_loop[0] + 1 ) * dx[0] ) / edge_v[0];
            if ( edge_v[0] < 0._rt ) amrex::Swap(cx_lo, cx_hi);
        }
        //
        if ( std::abs(edge_v[1]) > eps ) {
            cy_lo = -( edge_p0[1] - static_cast<Real>( ind_loop[1]     ) * dx[1] ) / edge_v[1];
            cy_hi = -( edge_p0[1] - static_cast<Real>( ind_loop[1] + 1 ) * dx[1] ) / edge_v[1];
            if ( edge_v[1] < 0._rt ) amrex::Swap(cy_lo, cy_hi);
        }
        //
        if ( std::abs(edge_v[2]) > eps ) {
            cz_lo = -( edge_p0[2] - static_cast<Real>( ind_loop[2]     ) * dx[2] ) / edge_v[2];
            cz_hi = -( edge_p0[2] - static_cast<Real>( ind_loop[2] + 1 ) * dx[2] ) / edge_v[2];
            if ( edge_v[2] < 0._rt ) amrex::Swap(cz_lo, cz_hi);
        }
        //
        Real lambda_min = amrex::max(cx_lo, cy_lo, cz_lo);
        Real lambda_max = amrex::min(cx_hi, cy_hi, cz_hi);

        if (lambda_tmp < lambda_min) {
            lambda_tmp = lambda_min;
        } else if ( lambda_tmp > lambda_max) {
            lambda_tmp = lambda_max;
        }

        RealVect rc_vec;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            c_vec_tmp[d] = edge_p0[d] + lambda_tmp*edge_v[d];
            rc_vec[d] = c_vec_tmp[d] - r_vec[d];
        }

        // determine new distance to particle
        Real min_dist_tmp = rc_vec.dotProduct(rc_vec);

        // minimize distance
        if (min_dist_tmp < min_dist) {
            min_dist = min_dist_tmp;
            c_vec = c_vec_tmp;
        }
#else
        RealVect c_vec_tmp = edge_p0;
        RealVect rc_vec = c_vec_tmp - r_vec;

        // determine new distance to particle
        Real min_dist_tmp = rc_vec.dotProduct(rc_vec);

        // minimize distance
        if (min_dist_tmp < min_dist) {
            min_dist = min_dist_tmp;
            c_vec = c_vec_tmp;
        }
#endif
    }

    return c_vec;
}
}

void FillSignedDistance (MultiFab& mf, EB2::Level const& ls_lev,
                         EBFArrayBoxFactory const& eb_factory, int refratio,
                         bool fluid_has_positive_sign)
{
    AMREX_ALWAYS_ASSERT(mf.is_nodal());

    ls_lev.fillLevelSet(mf, ls_lev.Geom()); // This is the implicit function, not the SDF.

    const auto& bndrycent = eb_factory.getBndryCent();
    const auto& areafrac = eb_factory.getAreaFrac();
    const auto& flags = eb_factory.getMultiEBCellFlagFab();
    const int eb_pad = bndrycent.nGrow();

    const auto dx_ls = ls_lev.Geom().CellSizeArray();
    const auto dx_eb = eb_factory.Geom().CellSizeArray();
    Real dx_eb_max = amrex::max(AMREX_D_DECL(dx_eb[0],dx_eb[1],dx_eb[2]));
    Real ls_roof = amrex::min(AMREX_D_DECL(dx_eb[0],dx_eb[1],dx_eb[2])) * static_cast<Real>(flags.nGrow()+1);

    Real fluid_sign = fluid_has_positive_sign ? 1._rt : -1._rt;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Box const& gbx = mfi.fabbox();
        Array4<Real> const& fab = mf.array(mfi);

        if (bndrycent.ok(mfi))
        {
            const auto& flag = flags.const_array(mfi);

            Box eb_search = mfi.validbox();
            eb_search.coarsen(refratio).enclosedCells().grow(eb_pad);

            const auto nallcells = static_cast<int>(eb_search.numPts());

            Gpu::DeviceVector<int> is_cut(nallcells);
            int* p_is_cut = is_cut.data();

            Gpu::DeviceVector<int> cutcell_offset(nallcells);
            int* p_cutcell_offset = cutcell_offset.data();

            int ncutcells = Scan::PrefixSum<int>
                (nallcells,
                 [=] AMREX_GPU_DEVICE (int icell) -> int
                 {
                     GpuArray<int,3> ijk = eb_search.atOffset3d(icell);
                     int is_cut_cell = flag(ijk[0],ijk[1],ijk[2]).isSingleValued();
                     p_is_cut[icell] = is_cut_cell;
                     return is_cut_cell;
                 },
                 [=] AMREX_GPU_DEVICE (int icell, int const& x)
                 {
                     p_cutcell_offset[icell] = x;
                 },
                 Scan::Type::exclusive, Scan::retSum);

            if (ncutcells > 0) {
                Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM*2> > facets(ncutcells);
                auto *p_facets = facets.data();
                Array4<Real const> const& bcent = bndrycent.const_array(mfi);
                AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                             Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                             Array4<Real const> const& apz = areafrac[2]->const_array(mfi));
                amrex::ParallelFor(eb_search, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    auto icell = eb_search.index(IntVect(AMREX_D_DECL(i,j,k)));
                    if (p_is_cut[icell]) {
                        GpuArray<Real,AMREX_SPACEDIM*2>& fac = p_facets[p_cutcell_offset[icell]];
                        AMREX_D_TERM(fac[0] = (bcent(i,j,k,0)+Real(i)+0.5_rt) * dx_eb[0];,
                                     fac[1] = (bcent(i,j,k,1)+Real(j)+0.5_rt) * dx_eb[1];,
                                     fac[2] = (bcent(i,j,k,2)+Real(k)+0.5_rt) * dx_eb[2]);

                         Real axm = apx(i,  j  , k  );
                         Real axp = apx(i+1,j  , k  );
                         Real aym = apy(i,  j  , k  );
                         Real ayp = apy(i,  j+1, k  );
#if (AMREX_SPACEDIM == 3)
                         Real azm = apz(i,  j  , k  );
                         Real azp = apz(i,  j  , k+1);
                         Real apnorm = std::sqrt((axm-axp)*(axm-axp) +
                                                 (aym-ayp)*(aym-ayp) +
                                                 (azm-azp)*(azm-azp));
#else
                         Real apnorm = std::sqrt((axm-axp)*(axm-axp) +
                                                 (aym-ayp)*(aym-ayp));
#endif
                         Real apnorminv = 1._rt / apnorm;
                         AMREX_D_TERM(Real anrmx = (axp-axm) * apnorminv;,   // pointing to the wall
                                      Real anrmy = (ayp-aym) * apnorminv;,
                                      Real anrmz = (azp-azm) * apnorminv);

                         // pointing to the fluid
                         AMREX_D_TERM(fac[AMREX_SPACEDIM+0] = -anrmx;,
                                      fac[AMREX_SPACEDIM+1] = -anrmy;,
                                      fac[AMREX_SPACEDIM+2] = -anrmz);
                    }
                });

                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    AMREX_D_TERM(Real dxinv = 1._rt/dx_eb[0];,
                                 Real dyinv = 1._rt/dx_eb[1];,
                                 Real dzinv = 1._rt/dx_eb[2]);
                    AMREX_D_TERM(Real x = i*dx_ls[0];,
                                 Real y = j*dx_ls[1];,
                                 Real z = k*dx_ls[2]);
                    Real min_dist2 = std::numeric_limits<Real>::max();
                    int i_nearest = 0;
                    for (int ifac  = 0; ifac < ncutcells; ++ifac) {
                        AMREX_D_TERM(Real cx = p_facets[ifac][0];,
                                     Real cy = p_facets[ifac][1];,
                                     Real cz = p_facets[ifac][2]);
                        Real dist2 = AMREX_D_TERM((x-cx)*(x-cx),+(y-cy)*(y-cy),+(z-cz)*(z-cz));
                        if (dist2 < min_dist2) {
                            i_nearest = ifac;
                            min_dist2 = dist2;
                        }
                    }

                    // Test if pos "projects onto" the nearest EB facet's interior
                    AMREX_D_TERM(Real cx = p_facets[i_nearest][0];,
                                 Real cy = p_facets[i_nearest][1];,
                                 Real cz = p_facets[i_nearest][2]);
                    AMREX_D_TERM(Real nx = p_facets[i_nearest][AMREX_SPACEDIM+0];,
                                 Real ny = p_facets[i_nearest][AMREX_SPACEDIM+1];,
                                 Real nz = p_facets[i_nearest][AMREX_SPACEDIM+2]);
                    Real dist_proj = AMREX_D_TERM((x-cx)*(-nx),+(y-cy)*(-ny),+(z-cz)*(-nz));
                    AMREX_D_TERM(Real eb_min_x = x + nx*dist_proj;,
                                 Real eb_min_y = y + ny*dist_proj;,
                                 Real eb_min_z = z + nz*dist_proj);
                    AMREX_D_TERM(int vi_cx = static_cast<int>(std::floor(cx * dxinv));,
                                 int vi_cy = static_cast<int>(std::floor(cy * dyinv));,
                                 int vi_cz = static_cast<int>(std::floor(cz * dzinv)));
                    AMREX_D_TERM(int vi_x = static_cast<int>(std::floor(eb_min_x * dxinv));,
                                 int vi_y = static_cast<int>(std::floor(eb_min_y * dyinv));,
                                 int vi_z = static_cast<int>(std::floor(eb_min_z * dzinv)));

                    bool min_pt_valid = false;
                    if ((AMREX_D_TERM(vi_cx == vi_x, && vi_cy == vi_y, && vi_cz == vi_z))  ||
                        std::abs(dist_proj) > ls_roof + dx_eb_max)
                    {
                        // If the distance is very big, we can set it to true as well.
                        // Later the signed distance will be assigned the roof value.
                        min_pt_valid = true;
                    } else { // rounding error might give false negatives
#if (AMREX_SPACEDIM == 3)
                        for (int k_shift = -1; k_shift <= 1; ++k_shift) {
#endif
                        for (int j_shift = -1; j_shift <= 1; ++j_shift) {
                        for (int i_shift = -1; i_shift <= 1; ++i_shift) {
                            AMREX_D_TERM(vi_x = static_cast<int>(std::floor((eb_min_x+i_shift*1.e-6_rt*dx_eb[0])*dxinv));,
                                         vi_y = static_cast<int>(std::floor((eb_min_y+j_shift*1.e-6_rt*dx_eb[1])*dyinv));,
                                         vi_z = static_cast<int>(std::floor((eb_min_z+k_shift*1.e-6_rt*dx_eb[2])*dzinv)));
                            if (AMREX_D_TERM(vi_cx == vi_x, && vi_cy == vi_y, && vi_cz == vi_z)) {
                                min_pt_valid = true;
                                goto after_loops;
                            }
                        }}
#if (AMREX_SPACEDIM == 3)
                        }
#endif
                        after_loops:;
                    }

                    // If projects onto nearest EB facet, then return projected distance
                    // Alternatively: find the nearest point on the EB edge
                    Real min_dist;
                    if ( min_pt_valid ) {
                        // this is a signed distance function
                        min_dist = dist_proj;
                    } else {
                        // fallback: find the nearest point on the EB edge
                        // revert the value of vi_x, vi_y and vi_z
                        AMREX_D_TERM(vi_x = static_cast<int>(std::floor(eb_min_x * dxinv));,
                                     vi_y = static_cast<int>(std::floor(eb_min_y * dyinv));,
                                     vi_z = static_cast<int>(std::floor(eb_min_z * dzinv)));
                        auto c_vec = detail::facets_nearest_pt
                            ({AMREX_D_DECL(vi_x,vi_y,vi_z)}, {AMREX_D_DECL(vi_cx, vi_cy, vi_cz)},
                             {AMREX_D_DECL(x,y,z)}, {AMREX_D_DECL(nx,ny,nz)},
                             {AMREX_D_DECL(cx,cy,cz)}, dx_eb);
                        Real min_edge_dist2 = AMREX_D_TERM( (c_vec[0]-x)*(c_vec[0]-x),
                                                           +(c_vec[1]-y)*(c_vec[1]-y),
                                                           +(c_vec[2]-z)*(c_vec[2]-z));
                        min_dist = -std::sqrt(amrex::min(min_dist2, min_edge_dist2));
                    }

                    Real usd = amrex::min(ls_roof,std::abs(min_dist));
                    if (fab(i,j,k) <= 0._rt) {
                        fab(i,j,k) = fluid_sign * usd;
                    } else {
                        fab(i,j,k) = (-fluid_sign) * usd;
                    }
                });
                Gpu::streamSynchronize();
            }
        } else {
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (fab(i,j,k) <= 0._rt) {
                    fab(i,j,k) = fluid_sign * ls_roof;
                } else {
                    fab(i,j,k) = (-fluid_sign) * ls_roof;
                }
            });
        }
    }

    mf.FillBoundary(0,1,ls_lev.Geom().periodicity());
}

} // end namespace
