
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MacProjector.H>
#include "initEB.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int verbose = 1;
        int n_cell = 128;
        int max_grid_size = 32;
        int is_periodic = 0;

        // read parameters
        {
            ParmParse pp;
            pp.query("v", verbose);
            pp.query("verbose", verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("is_periodic", is_periodic);
        }

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(is_periodic,is_periodic,is_periodic)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
            geom.define(domain);
            
            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        initEB(geom);

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, {2,2,2}, EBSupport::full);

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)),
                             dmap, 1, 1, MFInfo(), factory);
            vel[idim].setVal(2.0*(idim+1));
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)),
                              dmap, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }

        MultiFab divu(grids, dmap, 1, 0, MFInfo(), factory);

        AMREX_D_TERM(amrex::VisMF::Write(vel[0], "vx-pre");,
                     amrex::VisMF::Write(vel[1], "vy-pre");,
                     amrex::VisMF::Write(vel[2], "vz-pre"););
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }
        EB_computeDivergence(divu, amrex::GetArrOfConstPtrs(vel), geom);
//        amrex::VisMF::Write(divu, "divu-pre");
        amrex::Print() << "\nmax-norm of divu before projection is " << divu.norm0() << "\n" << std::endl;

        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                             {amrex::GetArrOfConstPtrs(beta)}, // beta 
                             {geom});                          // Geometry

        macproj.setVerbose(verbose);

        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Neumann,
                                          LinOpBCType::Neumann)},
                            {AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Neumann,
                                          LinOpBCType::Neumann)});

        Real reltol = 1.e-12;
        macproj.project(reltol);

//        AMREX_D_TERM(amrex::VisMF::Write(vel[0], "vx");,
//                     amrex::VisMF::Write(vel[1], "vy");,
//                     amrex::VisMF::Write(vel[2], "vz"););

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }
        EB_computeDivergence(divu, amrex::GetArrOfConstPtrs(vel), geom);
//        amrex::VisMF::Write(divu, "divu");
        amrex::Print() << "\nmax-norm of divu after projection is " << divu.norm0() << "\n" << std::endl;
    }

    amrex::Finalize();
}
