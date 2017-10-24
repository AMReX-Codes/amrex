
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <fort.H>
#include <prob_par.H>

using namespace amrex;

namespace prob {

    Real a     = 1.e-3;
    Real b     = 1.0;
    Real sigma = 10.0;
    Real w     = 0.05;

    MLLinOp::BCType bc_type = MLLinOp::BCType::Dirichlet;
}

using namespace prob;

void init_prob_parms ()
{
    ParmParse pp("prob");
    pp.query("a"    , a);
    pp.query("b"    , b);
    pp.query("sigma", sigma);
    pp.query("w"    , w);
    
    std::string bc_type_s;
    pp.query("bc_type", bc_type_s);
    if (bc_type_s == "Dirichlet") {
        bc_type = MLLinOp::BCType::Dirichlet;
    }
    else if (bc_type_s == "Neumann") {
        bc_type = MLLinOp::BCType::Neumann;
    }
    else if (bc_type_s == "Periodic") {
        bc_type = MLLinOp::BCType::Periodic;
    }
    else {
        amrex::Print() << "Don't know this boundary type: " << bc_type_s << "\n";
        amrex::Error("");
    }
}

void init_prob (const Vector<Geometry>& geom, Vector<MultiFab>& alpha, Vector<MultiFab>& beta,
                Vector<MultiFab>& rhs, Vector<MultiFab>& exact)
{
    char bct;
    if (bc_type == MLLinOp::BCType::Dirichlet) {
        bct = 'd';
    } else if (bc_type == MLLinOp::BCType::Neumann) {
        bct = 'n';
    } else {
        bct = 'p';
    }

    const int nlevels = geom.size();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        const Real* problo = geom[ilev].ProbLo();
        const Real* probhi = geom[ilev].ProbHi();
        const Real* dx     = geom[ilev].CellSize();

        for (MFIter mfi(alpha[ilev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            fort_set_coef(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_ANYD(exact[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(alpha[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(beta[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                          dx, problo, probhi,
                          &a, &b, &sigma, &w, &bct);
        }
    }
}
