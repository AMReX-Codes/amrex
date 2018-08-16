#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <Restart.H>

using namespace amrex;


extern "C"
{

void writecheckpointfile (int stepno[], int finest_level, Real dt[], Real time[], BoxArray** ba, MultiFab** phi)
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, time, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010

    const std::string& checkpointname = amrex::Concatenate(chk_file,stepno[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				                        std::ofstream::trunc |
				                        std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

       // write out title line
       HeaderFile << "Checkpoint file for AdvectionScalar_F\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of stepno
       for (int i = 0; i < nlevels; ++i) {
           HeaderFile << stepno[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < nlevels; ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of time
       for (int i = 0; i < nlevels; ++i) {
           HeaderFile << time[i] << " ";
       }
       HeaderFile << "\n";

      // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
	   BoxArray& ba1 = *ba[lev];
           ba1.writeOn(HeaderFile);
           HeaderFile << '\n';
       }

   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab& phi1 = *phi[lev];
        VisMF::Write(phi1,amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
    }	


}

void readheaderandboxarraydata (int finest_level[], int stepno[], Real dt[], Real time[], BoxArray** ba, DistributionMapping** dm )
{

    ParmParse pp("amr");
    pp.query("restart",restart_chkfile);

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level[0];
    GotoNextLine(is);

   // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            stepno[i++] = std::stoi(word);
        }
    }

   // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

   // read in array of time
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            time[i++] = std::stod(word);
        }
    }

// read in level 'lev' BoxArray from Header
     for (int lev = 0; lev <= finest_level[0]; ++lev) {
        BoxArray& ba1=*ba[lev];
        ba1.readFrom(is);
        ba[lev]=&ba1;
        GotoNextLine(is);

        //create a distribution mapping
        DistributionMapping dm1 { ba1, ParallelDescriptor::NProcs() };
        *dm[lev]=dm1;
        }

}

void readmultifabdata(int finest_level, MultiFab** phi) 

{
    ParmParse pp("amr");
    pp.query("restart",restart_chkfile);

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(*phi[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}

void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

}
