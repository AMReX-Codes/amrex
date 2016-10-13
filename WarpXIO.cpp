
#include <MultiFabUtil.H>
#include <PlotFileUtil.H>

#include <WarpX.H>

void
WarpX::WritePlotFile (int istep, Real t) const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "    Writing plotfile " << istep << std::endl;
    }
    
    int ncomp = 3*3;
    int ngrow = 0;
    int lev = 0;
    MultiFab mf(ba_arr[lev], ncomp, ngrow, dmap_arr[lev]);

    std::vector<std::string> varnames {"jx", "jy", "jz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"};

    std::vector<MultiFab*> srcmf(3);
    for (int i = 0; i < 3; ++i) {
	srcmf[i] = current[i].get();
    }
    int dcomp = 0;
    BoxLib::average_edge_to_cellcenter(mf, dcomp, srcmf);

    for (int i = 0; i < 3; ++i) {
	srcmf[i] = Efield[i].get();
    }
    dcomp += 3;
    BoxLib::average_edge_to_cellcenter(mf, dcomp, srcmf);

    for (int i = 0; i < 3; ++i) {
	srcmf[i] = Bfield[i].get();
    }
    dcomp += 3;
    BoxLib::average_face_to_cellcenter(mf, dcomp, srcmf);
    
    const std::string& plotfilename = BoxLib::Concatenate("plt",istep);
    BoxLib::WriteSingleLevelPlotfile(plotfilename, mf, varnames, geom_arr[lev], t);

    mypc->Checkpoint(plotfilename, "particle");
}
