#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "IO.H"

using namespace amrex;

void
WritePlotFile (const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
               const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
               const MultiFab& jx, const MultiFab& jy, const MultiFab& jz,
               const Geometry& geom, Real time, int step)
{
    BL_PROFILE("WritePlotFile()");

    BoxArray grids = Ex.boxArray();
    grids.convert(IntVect(0,0,0));

    const DistributionMapping& dmap = Ex.DistributionMap();

    const std::string& plotfilename = amrex::Concatenate("plt", step);

    amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

    Vector<std::string> varnames;

    const int ngrow = 0;
    const int ncomp = 3*3;
    MultiFab mf(grids, dmap, ncomp, ngrow);

    Vector<const MultiFab*> srcmf(BL_SPACEDIM);

    int dcomp = 0;

    srcmf[0] = &jx;
    srcmf[1] = &jy;
    srcmf[2] = &jz;
    amrex::average_edge_to_cellcenter(mf, dcomp, srcmf);
    varnames.push_back("jx");
    varnames.push_back("jy");
    varnames.push_back("jz");
    dcomp += 3;

    srcmf[0] = &Ex;
    srcmf[1] = &Ey;
    srcmf[2] = &Ez;
    amrex::average_edge_to_cellcenter(mf, dcomp, srcmf);
    varnames.push_back("Ex");
    varnames.push_back("Ey");
    varnames.push_back("Ez");
    dcomp += 3;

    srcmf[0] = &Bx;
    srcmf[1] = &By;
    srcmf[2] = &Bz;
    amrex::average_face_to_cellcenter(mf, dcomp, srcmf);
    varnames.push_back("Bx");
    varnames.push_back("By");
    varnames.push_back("Bz");
    dcomp += 3;

    amrex::WriteSingleLevelPlotfile(plotfilename, mf, varnames, geom, time, step);
}
