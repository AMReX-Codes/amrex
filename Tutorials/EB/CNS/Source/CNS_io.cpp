
#include <CNS.H>
#include <AMReX_EBMultiFabUtil.H>

using namespace amrex;

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (do_reflux && level > 0) {
        flux_reg.define(grids, papa.boxArray(level-1),
                        dmap, papa.DistributionMap(level-1),
                        geom, papa.Geom(level-1),
                        papa.refRatio(level-1), level, NUM_STATE);
    }

    buildMetrics();
}

void 
CNS::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old) 
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
CNS::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    BL_PROFILE("CNS::writePlotFile()");

//    AmrLevel::writePlotFile(dir, os, how);
    
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++) {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++) {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType()) {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
        {
            derive_names.push_back(it->name());
            num_derive++;
	}
    }

    int n_data_items = plot_var_map.size() + num_derive + 1;

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (int i = 0; i < plot_var_map.size(); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for ( std::list<std::string>::iterator it = derive_names.begin();
	      it != derive_names.end(); ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }

        // volfrac
        os << "vfrac\n";

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (int i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string sLevel = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += sLevel;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(FullPath, 0755))
            amrex::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }

        // volfrac threshhold for amrvis
        if (level == parent->finestLevel()) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                os << "1.0e-6\n";
            }
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow, MFInfo(), Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (int i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
#ifdef BL_TESTING
        // to avoid fcompare failure
        if (typ == Cost_Type) {
            plotMF.setVal(0.0, cnt, 1, nGrow);
        }
#endif
	cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::iterator it = derive_names.begin();
	     it != derive_names.end(); ++it)
	{
            auto derive_dat = derive(*it,cur_time,nGrow);
            MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
	    cnt++;
	}
    }

    plotMF.setVal(0.0, cnt, 1, nGrow);

    MultiFab::Copy(plotMF,volFrac(),0,cnt,1,nGrow);

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}
