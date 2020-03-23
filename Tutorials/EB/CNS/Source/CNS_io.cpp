
#include <CNS.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

}

void 
CNS::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old) 
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
CNS::getPlotData(MultiFab& plot_data, std::vector<std::string>& plot_names)
{
    BL_PROFILE("CNS::getPlotData()");

    plot_names.resize(0);
    
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

    //
    // Names of variables -- first state, then derived
    //
    for (int i = 0; i < plot_var_map.size(); i++)
    {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        plot_names.push_back(desc_lst[typ].name(comp));
    }

    for ( std::list<std::string>::iterator it = derive_names.begin();
          it != derive_names.end(); ++it)
    {
        const DeriveRec* rec = derive_lst.get(*it);
        plot_names.push_back(rec->variableName(0));
    }

    // volfrac
    plot_names.push_back("vfrac");


    int       cnt   = 0;
    const int nGrow = 0;
    plot_data.define(grids,dmap,n_data_items,nGrow, MFInfo(), Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (int i = 0; i < plot_var_map.size(); i++)
    {
        int typ  = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plot_data,*this_dat,comp,cnt,1,nGrow);
#ifdef AMREX_TESTING
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
            MultiFab::Copy(plot_data,*derive_dat,0,cnt,1,nGrow);
            cnt++;
        }
    }

    plot_data.setVal(0.0, cnt, 1, nGrow);

    MultiFab::Copy(plot_data,volFrac(),0,cnt,1,nGrow);
}

void CNS::writePlotFilePost(const std::string &dir, std::ostream &os)
{
    if (ParallelDescriptor::IOProcessor() && (level == parent->finestLevel()))
    {
        // volfrac threshhold for amrvis
        for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
            os << "1.0e-6\n";
        }
    }
}

#ifdef AMREX_USE_HDF5
void CNS::checkPointHDF5Post()
{
    // save the inputs file to the checkpoint file

    if (level == parent->finestLevel()) {
      H5 &h5 = parent->getOutputHDF5();

      H5 cns = h5.createGroup("CNS");

      std::ostringstream cfg;
      ParmParse::dumpTable(cfg, true);
      cns.writeString("config", {cfg.str()});
      cns.closeGroup();
    }
}

void CNS::writePlotHDF5Post()
{
    // save the inputs file to the checkpoint file

    if (level == parent->finestLevel()) {
      H5 &h5 = parent->getOutputHDF5();

      H5 cns = h5.createGroup("CNS");

      std::ostringstream cfg;
      ParmParse::dumpTable(cfg, true);
      cns.writeString("config", {cfg.str()});
      cns.closeGroup();
    }
}
#endif
