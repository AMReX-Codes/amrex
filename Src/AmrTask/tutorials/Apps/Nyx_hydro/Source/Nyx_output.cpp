#include <unistd.h>
#include <iomanip>
#include <Nyx.H>
#include <Nyx_F.H>
#include "Nyx_output.H"

#include "AMReX_buildInfo.H"

#ifdef FORCING
#include "Forcing.H"

void mt_write(std::ofstream& output);
#endif

using namespace amrex;

namespace
{
    const std::string dm_chk_particle_file("DM");
    const std::string dm_plt_particle_file("DM");

    const std::string agn_chk_particle_file("AGN");
    const std::string agn_plt_particle_file("AGN");
}

std::string
Nyx::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");
    return the_plot_file_type;
}

std::string
Nyx::retrieveDM () {
    return dm_chk_particle_file;
}

#ifdef AGN
std::string
Nyx::retrieveAGN () {
    return agn_chk_particle_file;
}
#endif

void
Nyx::setPlotVariables ()
{
    AmrLevel::setPlotVariables();

    ParmParse pp("nyx");
    bool plot_X, plot_rank;
    if (pp.query("plot_rank", plot_rank))
    {
        if (plot_rank)
        {
            //
            // Write the processor ID for each grid into the plotfile
            //
            std::string proc_string = "Rank";
            parent->addDerivePlotVar(proc_string);
        }
    }
    if (pp.query("plot_X", plot_X))
    {
        if (plot_X)
        {
            //
            // Get the number of species from the network model.
            //
            fort_get_num_spec(&NumSpec);
            //
            // Get the species names from the network model.
            //
            for (int i = 0; i < NumSpec; i++)
            {
                int len = 20;
                Vector<int> int_spec_names(len);
                //
                // This call return the actual length of each string in "len"
                //
                fort_get_spec_names(int_spec_names.dataPtr(), &i, &len);
                char* spec_name = new char[len+1];

                for (int j = 0; j < len; j++)
                    spec_name[j] = int_spec_names[j];
                spec_name[len] = '\0';

                // @todo: better string ops
                std::string spec_string = "X(";
                spec_string += spec_name;
                spec_string += ')';
                parent->addDerivePlotVar(spec_string);
                delete [] spec_name;
            }
        }
    }
}

void
Nyx::writePlotFile (const std::string& dir,
                    ostream&           os,
                    VisMF::How         how)
{
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
    {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
        {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp))
                && desc_lst[typ].getType() == IndexType::TheCellType())
            {
                plot_var_map.push_back(std::pair<int,int>(typ, comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;

    for (std::list<DeriveRec>::const_iterator it = derive_lst.dlist().begin();
         it != derive_lst.dlist().end(); ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
        {
            if (it->name() == "particle_count" ||
                it->name() == "total_particle_count" ||
                it->name() == "particle_mass_density" ||
                it->name() == "total_density")
            {
                if (Nyx::theDMPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#ifdef AGN
            } else if (it->name() == "agn_particle_count" ||
                       it->name() == "agn_mass_density")
            {
                if (Nyx::theAPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#endif
#ifdef NEUTRINO_PARTICLES
            } else if (it->name() == "neutrino_particle_count" ||
                       it->name() == "neutrino_mass_density")
            {
                if (Nyx::theNPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#endif
            } else if (it->name() == "Rank") {
                derive_names.push_back(it->name());
                num_derive++;
            } else {
                derive_names.push_back(it->name());
                num_derive++;
            }
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0) {
            amrex::Error("Must specify at least one valid data item to plot");
	}

        os << n_data_items << '\n';
        //
        // Names of variables -- first state, then derived
        //
        for (i = 0; i < plot_var_map.size(); i++)
        {
            int typ = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            os << desc_lst[typ].name(comp) << '\n';
        }

        for (std::list<std::string>::iterator it = derive_names.begin();
             it != derive_names.end(); ++it)
        {
            const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	std::string PrettyLine = std::string(78, '=') + "\n";
	std::string OtherLine = std::string(78, '-') + "\n";
	std::string SkipSpace = std::string(8, ' ') + "\n";

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Nyx Job Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "inputs file: " << inputs_name << "\n\n";

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
	jobInfoFile << "\n";
	jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
	  getCPUTime()/3600.0;

	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << "\n";
	}

	jobInfoFile << "\n\n";


	// cosmology information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Cosmology Information\n";
	jobInfoFile << PrettyLine;

        //	Real comoving_OmM, comoving_OmL, comoving_h;
        //	fort_get_omm(&comoving_OmM);
	// Omega lambda is defined algebraically
	Real comoving_OmL = 1. - comoving_OmM;

	// fort_get_hubble(&comoving_h);

	jobInfoFile << "Omega_m (comoving):      " << comoving_OmM << "\n";
	jobInfoFile << "Omega_lambda (comoving): " << comoving_OmL << "\n";
	jobInfoFile << "h (comoving):            " << comoving_h << "\n";

	jobInfoFile << "\n\n";

        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
	jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
	jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
	jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "Nyx    git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "AMReX git hash:  " << githash2 << "\n";
	}

	jobInfoFile << "\n\n";

	// grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (i = 0; i <= f_lev; i++)
          {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (n = 0; n < BL_SPACEDIM; n++)
              {
                jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
                //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
              }
            jobInfoFile << "\n\n";
          }

        jobInfoFile << " Boundary conditions\n";
        Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
        ParmParse pp("nyx");
        pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
        pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


        // these names correspond to the integer flags setup in the
        // Castro_setup.cpp
        const char* names_bc[] =
          { "interior", "inflow", "outflow",
            "symmetry", "slipwall", "noslipwall" };


        jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
        jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
        if (BL_SPACEDIM >= 2) {
          jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
          jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
        }
        if (BL_SPACEDIM == 3) {
          jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
          jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
        }

        jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;

	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if( ! levelDirectoryCreated) {
      amrex::Print() << "IOIOIOIO:CD  Nyx::writePlotFile:  ! ldc:  creating directory:  "
                     << FullPath << '\n';
      if (ParallelDescriptor::IOProcessor()) {
        if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
	}
      }
      //
      // Force other processors to wait until directory is built.
      //
      ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int cnt = 0;
    const int nGrow = 0;
    MultiFab plotMF(grids, dmap, n_data_items, nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
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
            const auto& derive_dat = derive(*it, cur_time, nGrow);
            MultiFab::Copy(plotMF, *derive_dat, 0, cnt, 1, nGrow);
            cnt++;
        }
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF, TheFullPath, how, true);

    //
    // Write the particles and `comoving_a` in a plotfile directory. 
    //
    particle_plot_file(dir);

    // Write out all parameters into the plotfile
    if (write_parameters_in_plotfile) {
	write_parameter_file(dir);
    }

    if(Nyx::theDMPC()) {
      Nyx::theDMPC()->SetLevelDirectoriesCreated(false);
    }
#ifdef AGN
    if(Nyx::theAPC()) {
      Nyx::theAPC()->SetLevelDirectoriesCreated(false);
    }
#endif

}

void
Nyx::writePlotFilePre (const std::string& dir, ostream& os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePre();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePre();
  }
#endif

}


void
Nyx::writePlotFilePost (const std::string& dir, ostream& os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePost();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePost();
  }
#endif
}


void
Nyx::particle_plot_file (const std::string& dir)
{
    if (level == 0)
    {
        if (Nyx::theDMPC())
          {
            Nyx::theDMPC()->WriteNyxPlotFile(dir, dm_plt_particle_file);
          }

#ifdef AGN
        if (Nyx::theAPC())
          {
            Nyx::theAPC()->WriteNyxPlotFile(dir, agn_plt_particle_file);
          }
#endif

#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif

        // Write comoving_a into its own file in the particle directory
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            if (cur_time == 0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
            File.close();
        }

        // Write particle_plotfile_format into its own file in the particle directory
        if (Nyx::theDMPC() && ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/" + dm_plt_particle_file + "/precision";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            File << particle_plotfile_format << '\n';
            File.close();
        }

#ifdef AGN
        // Write particle_plotfile_format into its own file in the particle directory
        if (Nyx::theAPC() && ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/" + agn_plt_particle_file + "/precision";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            File << particle_plotfile_format << '\n';
            File.close();
        }
#endif
    }
}

void
Nyx::particle_check_point (const std::string& dir)
{
  BL_PROFILE("Nyx::particle_check_point");
  if (level == 0)
    {
      if (Nyx::theDMPC())
        {
          Nyx::theDMPC()->NyxCheckpoint(dir, dm_chk_particle_file);
        }
#ifdef AGN
      if (Nyx::theAPC())
        {
          Nyx::theAPC()->NyxCheckpoint(dir, agn_chk_particle_file);
        }
#endif

#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif

        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            if (cur_time == 0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
        }
    }
}

void
Nyx::write_parameter_file (const std::string& dir)
{
    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/the_parameters";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            ParmParse::dumpTable(File,true);
            File.close();
        }
    }
}

void
Nyx::writeMultiFabAsPlotFile(const std::string& pltfile,
                             const MultiFab&    mf,
                             std::string        componentName)
{
    std::ofstream os;
    if (ParallelDescriptor::IOProcessor())
    {
        if( ! amrex::UtilCreateDirectory(pltfile, 0755)) {
          amrex::CreateDirectoryFailed(pltfile);
	}
        std::string HeaderFileName = pltfile + "/Header";
        os.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        // The first thing we write out is the plotfile type.
        os << thePlotFileType() << '\n';
        // Just one component ...
        os << 1 << '\n';
        // ... with name
        os << componentName << '\n';
        // Dimension
        os << BL_SPACEDIM << '\n';
        // Time
        os << "0\n";
        // One level
        os << "0\n";
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        // Only one level -> no refinement ratios
        os << '\n';
        // Geom
        os << parent->Geom(0).Domain() << ' ';
        os << '\n';
        os << parent->levelSteps(0) << ' ';
        os << '\n';
        for (int k = 0; k < BL_SPACEDIM; k++)
            os << parent->Geom(0).CellSize()[k] << ' ';
        os << '\n';
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    std::string Level = "Level_0";
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = pltfile;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
        if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
	}
    }
    //
    // Force other processors to wait until directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        std::string PathNameInHeader = Level;
        PathNameInHeader += BaseName;
        os << PathNameInHeader << '\n';
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(mf, TheFullPath);
    ParallelDescriptor::Barrier();
}

void
Nyx::checkPoint (const std::string& dir,
                 std::ostream&      os,
                 VisMF::How         how,
                 bool               dump_old_default)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
  particle_check_point(dir);
#ifdef FORCING
  forcing_check_point(dir);
#endif

  if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      {
	// store elapsed CPU time
	std::ofstream CPUFile;
	std::string FullPathCPUFile = dir;
	FullPathCPUFile += "/CPUtime";
	CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

	CPUFile << std::setprecision(15) << getCPUTime();
	CPUFile.close();
      }
    }

    if(Nyx::theDMPC()) {
      Nyx::theDMPC()->SetLevelDirectoriesCreated(false);
    }
#ifdef AGN
    if(Nyx::theAPC()) {
      Nyx::theAPC()->SetLevelDirectoriesCreated(false);
    }
#endif

}

void
Nyx::checkPointPre (const std::string& dir,
                    std::ostream&      os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->CheckpointPre();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->CheckpointPre();
  }
#endif

}


void
Nyx::checkPointPost (const std::string& dir,
                 std::ostream&      os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->CheckpointPost();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->CheckpointPost();
  }
#endif
}


#ifdef FORCING
void
Nyx::forcing_check_point (const std::string& dir)
{
    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/forcing";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.setf(std::ios::scientific, std::ios::floatfield);
            File.precision(16);
            forcing->write_Spectrum(File);
            File.close();

            FileName = dir + "/mt";
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            mt_write(File);
        }
    }
}
#endif
