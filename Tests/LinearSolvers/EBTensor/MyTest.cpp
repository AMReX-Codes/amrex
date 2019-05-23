#include "MyTest.H"
#include "MyTest_K.H"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBTensorOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_VisMF.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    RealBox rb({AMREX_D_DECL(-1.0,-1.0,-1.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom.define(domain);

    {
        bool static first = true;
        if (first) {
            first = false;
            EB2::Build(geom, 0, 100);
        }
    }

    initData();
}

void
MyTest::solve ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = 1.e-11;
    const Real tol_abs = 0.0;

    MLEBTensorOp ebtensorop({geom}, {grids}, {dmap}, info, {factory.get()});

    ebtensorop.setMaxOrder(linop_maxorder);

    Array<LinOpBCType,AMREX_SPACEDIM> v_lo_bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann)};
    Array<LinOpBCType,AMREX_SPACEDIM> v_hi_bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann)};

    std::string geom_type;
    int cylinder_direction;
    {
        ParmParse pp("eb2");
        pp.get("geom_type", geom_type);
        pp.get("cylinder_direction", cylinder_direction);
    }

    if (geom_type == "all_regular") {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            v_lo_bc[idim] = LinOpBCType::Dirichlet;
            v_hi_bc[idim] = LinOpBCType::Dirichlet;
        }
    } else {
        v_lo_bc[cylinder_direction] = LinOpBCType::Dirichlet;
        v_hi_bc[cylinder_direction] = LinOpBCType::Dirichlet;
    }

    ebtensorop.setDomainBC({AMREX_D_DECL(v_lo_bc,v_lo_bc,v_lo_bc)},
                           {AMREX_D_DECL(v_hi_bc,v_hi_bc,v_hi_bc)});

    ebtensorop.setLevelBC(0, &solution);

    const Real a = 0.0; // 1.0e6;
    {
        MultiFab tmp(grids, dmap, 1, 0, MFInfo(), *factory);
        tmp.setVal(a);

        ebtensorop.setACoeffs(0, tmp);

        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(grids, IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, dmap, 1, 0, MFInfo(), *factory);
        }
        amrex::average_cellcenter_to_face(amrex::GetArrOfPtrs(face_bcoef), eta, geom);
        ebtensorop.setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef));
        ebtensorop.setEBShearViscosity(0, eta);
    }

    MLMG mlmg(ebtensorop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    mlmg.setBottomTolerance(1.e-4);

    MultiFab::Saxpy(rhs, a, exact, 0, 0, AMREX_SPACEDIM, 0);

//    solution.setVal(0.0);
    mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);

    MultiFab error(grids, dmap, 1, 0, MFInfo(), *factory);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        amrex::Print() << "\n";
        MultiFab::Copy(error, solution, idim, 0, 1, 0);
        MultiFab::Subtract(error, exact, idim, 0, 1, 0);
        amrex::Print() << "  max-norm error = " << error.norm0() << std::endl;
        const MultiFab& vfrc = factory->getVolFrac();
        MultiFab::Multiply(error, vfrc, 0, 0, 1, 0);
        const auto dx = geom.CellSize();
        error.mult(dx[0]*dx[1]*dx[2]);
        amrex::Print() << "    1-norm error = " << error.norm1() << std::endl;
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
 
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);
}

void
MyTest::initData ()
{
    grids.define(geom.Domain());
    grids.maxSize(max_grid_size);
    dmap.define(grids);

    factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);

    solution.define(grids, dmap, 3, 1, MFInfo(), *factory);
    exact.define(grids, dmap, 3, 1, MFInfo(), *factory);
    rhs.define(grids, dmap, 3, 1, MFInfo(), *factory);
    eta.define(grids, dmap, 1, 1, MFInfo(), *factory);

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLo();

    int cylinder_direction;
    amrex::Real R;
    ParmParse pp("eb2");
    pp.get("cylinder_direction", cylinder_direction);
    pp.get("cylinder_radius", R);
    amrex::Real R2 = R*R;

    for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const Array4<Real> velfab = exact.array(mfi);
        const Array4<Real> rhsfab = rhs.array(mfi);
        const Array4<Real> etafab = eta.array(mfi);
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real x = (i+0.5)*dx[0] + problo[0];
                    Real y = (j+0.5)*dx[1] + problo[1];
                    Real z = (k+0.5)*dx[2] + problo[2];
                    Real u,v,w,urhs,vrhs,wrhs,eta;
                    if (cylinder_direction == 2) {
                        init(x,y,z,R2,u,v,w,urhs,vrhs,wrhs,eta);
                    } else if (cylinder_direction == 0) {
                        init(y,z,x,R2,v,w,u,vrhs,wrhs,urhs,eta);
                    } else {
                        init(z,x,y,R2,w,u,v,wrhs,urhs,vrhs,eta);
                    }
                    velfab(i,j,k,0) = u;
                    velfab(i,j,k,1) = v;
                    velfab(i,j,k,2) = w;
                    rhsfab(i,j,k,0) = urhs;
                    rhsfab(i,j,k,1) = vrhs;
                    rhsfab(i,j,k,2) = wrhs;
                    etafab(i,j,k) = eta;
                    if (x < -1.0 or x > 1.0 or
                        y < -1.0 or y > 1.0 or
                        z < -1.0 or z > 1.0)
                    {
                        x = std::max(-1.0,std::min(1.0,x));
                        y = std::max(-1.0,std::min(1.0,y));
                        z = std::max(-1.0,std::min(1.0,z));
                        if (cylinder_direction == 2) {
                            init(x,y,z,R2,u,v,w,urhs,vrhs,wrhs,eta);
                        } else if (cylinder_direction == 0) {
                            init(y,z,x,R2,v,w,u,vrhs,wrhs,urhs,eta);
                        } else {
                            init(z,x,y,R2,w,u,v,wrhs,urhs,vrhs,eta);
                        }
                        velfab(i,j,k,0) = u;
                        velfab(i,j,k,1) = v;
                        velfab(i,j,k,2) = w;
                    }
                }
            }
        }
    }

    amrex::EB_set_covered(exact, 0.0);
    amrex::EB_set_covered(rhs, 0.0);
    amrex::EB_set_covered(eta, 0.0);

    solution.setVal(0.0);
    MultiFab::Copy(solution, exact, 0, 0, AMREX_SPACEDIM, 1);
}
