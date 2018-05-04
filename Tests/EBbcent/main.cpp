
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WrappedGShop.H>
#include <AMReX_EBLevelGrid.H>
#include <AMReX_PlaneIF.H>

using namespace amrex;

static IntVect debug_cell(AMREX_D_DECL(32, 6, 0));

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int n_cell = 64;
        int max_grid_size = 32;
        std::string geom_type = "ramp_normal_point";

        BoxArray grids;
        DistributionMapping dmap;
        Geometry geom;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
            std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
            Geometry::Setup(&rb, 0, is_periodic.data());
            Box domain(IntVect(0), IntVect(n_cell-1));
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        {
            Vector<Real> pointvec {AMREX_D_DECL( 0.5 ,0.0, 0.0)};
            Vector<Real> normalvec{AMREX_D_DECL(-1.0, 0.0, 0.0)};

            const Real* dx = geom.CellSize();
            pointvec[0] += 0.1*dx[0];

            int inside = 0;
            bool normalInside = (inside == 0);
            RealVect normal{AMREX_D_DECL(normalvec[0],normalvec[1],normalvec[2])};
            RealVect point {AMREX_D_DECL( pointvec[0], pointvec[1], pointvec[2])};
            
            PlaneIF pif(normal, point, normalInside);

            bool verbose = true;
            WrappedGShop gshop(pif, verbose);
            AMReX_EBIS::instance()->define(geom.Domain(), RealVect::Zero, dx[0], gshop,
                                           max_grid_size, 0);
        }

        {
            EBLevelGrid eblg(grids, dmap, geom.Domain(), 1);
            const EBISLayout& ebisl = eblg.getEBISL();
            for (MFIter mfi(grids,dmap); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.contains(debug_cell))
                {
                    const EBISBox& ebisbox = ebisl[mfi];

                    const auto& vofs = ebisbox.getVoFs(debug_cell);
                    for (const auto& vi : vofs)
                    {
                        amrex::AllPrint().SetPrecision(17)
                            << "volFrac = " << ebisbox.volFrac(vi) << "\n";
                    }

                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        const auto& lo_faces = ebisbox.getAllFaces(debug_cell, idim, Side::Lo);
                        for (const auto& fi : lo_faces)
                        {
                            amrex::AllPrint().SetPrecision(17)
                                << "lo face of dim " << idim << ": areaFrac = "
                                << ebisbox.areaFrac(fi) << "\n"
                                << "                : centroid = "
                                << ebisbox.centroid(fi) << "\n";
                        }

                        const auto& hi_faces = ebisbox.getAllFaces(debug_cell, idim, Side::Hi);
                        for (const auto& fi : hi_faces)
                        {
                            amrex::AllPrint().SetPrecision(17)
                                << "hi face of dim " << idim << ": areaFrac = "
                                << ebisbox.areaFrac(fi) << "\n"
                                << "                : centroid = "
                                << ebisbox.centroid(fi) << "\n";
                        }
                    }

                    for (const auto& vi : vofs)
                    {
                        amrex::AllPrint().SetPrecision(17)
                            << "bndryCentroid = " << ebisbox.bndryCentroid(vi) << "\n";
                    }
                    for (const auto& vi : vofs)
                    {
                      const auto& ebmom = ebisbox.getEBMoments(vi);
                      for(MomItSpaceDim momit; momit.ok(); ++ momit)
                      {
                        amrex::AllPrint().SetPrecision(17)
                          << "EBMoment" << momit() << " = " << ebmom[momit()] << "\n";
                      }
                    }

                }
            }
        }
    }

    amrex::Finalize();
}
