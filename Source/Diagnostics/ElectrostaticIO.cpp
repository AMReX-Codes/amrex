#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>

#include <WarpX.H>
#include <WarpX_f.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

using namespace amrex;

void
WarpX::
WritePlotFileES (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                 const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                 const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E)
{
    BL_PROFILE("WarpX::WritePlotFileES()");

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(plotfile_headerversion);

    const std::string& plotfilename = amrex::Concatenate(plot_file,istep[0]);

    amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

    const int nlevels = finestLevel()+1;

    {
	Vector<std::string> varnames;
	Vector<std::unique_ptr<MultiFab> > mf(finest_level+1);

	for (int lev = 0; lev <= finest_level; ++lev) {
            int ncomp = 5;
            const int ngrow = 0;
            mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

            int dcomp = 0;
            amrex::average_node_to_cellcenter(*mf[lev], dcomp, *rho[lev], 0, 1);
            if (lev == 0) {
                varnames.push_back("rho");
            }
            dcomp += 1;

            amrex::average_node_to_cellcenter(*mf[lev], dcomp  , *E[lev][0], 0, 1);
            amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *E[lev][0], 0, 1);
            amrex::average_node_to_cellcenter(*mf[lev], dcomp+2, *E[lev][0], 0, 1);

            if (lev == 0) {
                varnames.push_back("Ex");
                varnames.push_back("Ey");
                varnames.push_back("Ez");
            }
            dcomp += 3;

            amrex::average_node_to_cellcenter(*mf[lev], dcomp, *phi[lev], 0, 1);
            if (lev == 0) {
                varnames.push_back("phi");
            }
            dcomp += 1;
        }

        Vector<std::string> rfs(1,"raw_fields"); // pre-build raw_fields/
        amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                       amrex::GetVecOfConstPtrs(mf),
                                       varnames, Geom(), t_new[0], istep, refRatio(),
                                       "HyperCLaw-V1.1",
                                       "Level_",
                                       "Cell",
                                       rfs);
    }

    {
        const std::string raw_plotfilename = plotfilename + "/raw_fields";
        const int nlevels = finestLevel()+1;
        for (int lev = 0; lev < nlevels; ++lev) {
            const DistributionMapping& dm = DistributionMap(lev);

            MultiFab Ex( E[lev][0]->boxArray(), dm, 1, 0);
            MultiFab Ey( E[lev][1]->boxArray(), dm, 1, 0);
            MultiFab Ez( E[lev][2]->boxArray(), dm, 1, 0);
            MultiFab charge_density(rho[lev]->boxArray(), dm, 1, 0);
            MultiFab potential(phi[lev]->boxArray(), dm, 1, 0);

            MultiFab::Copy(Ex, *E[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(Ey, *E[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(Ez, *E[lev][2], 0, 0, 1, 0);
            MultiFab::Copy(charge_density, *rho[lev], 0, 0, 1, 0);
            MultiFab::Copy(potential, *phi[lev], 0, 0, 1, 0);

            VisMF::Write(Ex, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex"));
            VisMF::Write(Ey, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey"));
            VisMF::Write(Ez, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez"));
            VisMF::Write(charge_density, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "rho"));
            VisMF::Write(potential, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "phi"));
        }
    }

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("weight");

    particle_varnames.push_back("momentum_x");
    particle_varnames.push_back("momentum_y");
    particle_varnames.push_back("momentum_z");

    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");
    particle_varnames.push_back("Ez");

    particle_varnames.push_back("Bx");
    particle_varnames.push_back("By");
    particle_varnames.push_back("Bz");

    Vector<std::string> int_names;
        
    mypc->Checkpoint(plotfilename, particle_varnames, int_names);

    WriteJobInfo(plotfilename);

    WriteWarpXHeader(plotfilename);

    VisMF::SetHeaderVersion(current_version);
}
