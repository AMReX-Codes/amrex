#include "FlushFormatPlotfile.H"
#include "WarpX.H"
#include "Interpolate.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

void
FlushFormatPlotfile::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<const amrex::MultiFab*> mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    MultiParticleContainer& mpc, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_rho, bool plot_F) const
{
    auto & warpx = WarpX::GetInstance();
    const auto step = iteration[0];
    const std::string& filename = amrex::Concatenate(prefix, step);
    amrex::Print() << "  Writing plotfile " << filename << "\n";

    Vector<std::string> rfs;
    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::Version_v1);
    if (plot_raw_fields) rfs.emplace_back("raw_fields");
    amrex::WriteMultiLevelPlotfile(filename, nlev,
                                   mf,
                                   varnames, geom,
                                   time, iteration, warpx.refRatio(),
                                   "HyperCLaw-V1.1",
                                   "Level_",
                                   "Cell",
                                   rfs
                                   );

    WriteAllRawFields(plot_raw_fields, nlev, filename, plot_raw_fields_guards,
                      plot_rho, plot_F);

    mpc.WritePlotFile(filename);

    WriteJobInfo(filename);

    WriteWarpXHeader(filename);

    VisMF::SetHeaderVersion(current_version);
}

void
FlushFormatPlotfile::WriteJobInfo(const std::string& dir) const
{

    auto & warpx = WarpX::GetInstance();

    if (ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = dir;

        std::string PrettyLine = std::string(78, '=') + "\n";
//        std::string OtherLine = std::string(78, '-') + "\n";
//        std::string SkipSpace = std::string(8, ' ') + "\n";

        FullPathJobInfoFile += "/warpx_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " WarpX Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

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
        const char* githash3 = buildInfoGetGitHash(3);
        if (strlen(githash1) > 0) {
            jobInfoFile << "WarpX  git describe: " << githash1 << "\n";
        }
        if (strlen(githash2) > 0) {
            jobInfoFile << "AMReX  git describe: " << githash2 << "\n";
        }
        if (strlen(githash3) > 0) {
            jobInfoFile << "PICSAR git describe: " << githash3 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= warpx.finestLevel(); i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << warpx.boxArray(i).size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                jobInfoFile << warpx.Geom(i).Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << " Boundary conditions\n";

        jobInfoFile << "   -x: " << "interior" << "\n";
        jobInfoFile << "   +x: " << "interior" << "\n";
        if (AMREX_SPACEDIM >= 2) {
            jobInfoFile << "   -y: " << "interior" << "\n";
            jobInfoFile << "   +y: " << "interior" << "\n";
        }
        if (AMREX_SPACEDIM == 3) {
            jobInfoFile << "   -z: " << "interior" << "\n";
            jobInfoFile << "   +z: " << "interior" << "\n";
        }

        jobInfoFile << "\n\n";


        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}

void
FlushFormatPlotfile::WriteWarpXHeader(const std::string& name) const
{
    auto & warpx = WarpX::GetInstance();
    if (ParallelDescriptor::IOProcessor())
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(name + "/WarpXHeader");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << "Checkpoint version: 1\n";

        const int nlevels = warpx.finestLevel()+1;
        HeaderFile << nlevels << "\n";

        for (int i = 0; i < warpx.getistep().size(); ++i) {
            HeaderFile << warpx.getistep(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.getnsubsteps().size(); ++i) {
            HeaderFile << warpx.getnsubsteps(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.gett_new().size(); ++i) {
            HeaderFile << warpx.gett_new(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.gett_old().size(); ++i) {
            HeaderFile << warpx.gett_old(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.getdt().size(); ++i) {
            HeaderFile << warpx.getdt(i) << " ";
        }
        HeaderFile << "\n";

        HeaderFile << warpx.getmoving_window_x() << "\n";

        HeaderFile << warpx.getis_synchronized() << "\n";

        // Geometry
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << warpx.Geom(0).ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << warpx.Geom(0).ProbHi(i) << ' ';
        }
        HeaderFile << '\n';

        // BoxArray
        for (int lev = 0; lev < nlevels; ++lev) {
            warpx.boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }

        warpx.GetPartContainer().WriteHeader(HeaderFile);
    }
}

/** \brief Write the data from MultiFab `F` into the file `filename`
 *  as a raw field (i.e. no interpolation to cell centers).
 *  Write guard cells if `plot_guards` is True.
 */
void
WriteRawMF ( const MultiFab& F, const DistributionMapping& dm,
             const std::string& filename,
             const std::string& level_prefix,
             const std::string& field_name,
             const int lev, const bool plot_guards )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);
    if (plot_guards) {
        // Dump original MultiFab F
        VisMF::Write(F, prefix);
    } else {
        // Copy original MultiFab into one that does not have guard cells
        MultiFab tmpF( F.boxArray(), dm, F.nComp(), 0);
        MultiFab::Copy(tmpF, F, 0, 0, F.nComp(), 0);
        VisMF::Write(tmpF, prefix);
    }
}

/** \brief Write a multifab of the same shape as `F` but filled with 0.
 *  (The shape includes guard cells if `plot_guards` is True.)
 *  This is mainly needed because the yt reader requires all levels of the
 *  coarse/fine patch to be written, but WarpX does not have data for
 *  the coarse patch of level 0 (meaningless).
 */
void
WriteZeroRawMF( const MultiFab& F, const DistributionMapping& dm,
                const std::string& filename,
                const std::string& level_prefix,
                const std::string& field_name,
                const int lev, const int ng )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);

    MultiFab tmpF(F.boxArray(), dm, F.nComp(), ng);
    tmpF.setVal(0.);
    VisMF::Write(tmpF, prefix);
}

/** \brief Write the coarse vector multifab `F*_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F*_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseVector( const std::string field_name,
    const MultiFab* Fx_cp,
    const MultiFab* Fy_cp,
    const MultiFab* Fz_cp,
    const MultiFab* Fx_fp,
    const MultiFab* Fy_fp,
    const MultiFab* Fz_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards )
{
    int ng = 0;
    if (plot_guards) ng = Fx_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawMF( *Fx_fp, dm, filename, level_prefix, field_name+"x_cp", lev, ng );
        WriteZeroRawMF( *Fy_fp, dm, filename, level_prefix, field_name+"y_cp", lev, ng );
        WriteZeroRawMF( *Fz_fp, dm, filename, level_prefix, field_name+"z_cp", lev, ng );
    } else {
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = Interpolate::getInterpolatedVector( Fx_cp, Fy_cp, Fz_cp, Fx_fp, Fy_fp, Fz_fp,
                                    dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawMF( *F[0], dm, filename, level_prefix, field_name+"x_cp", lev, plot_guards );
        WriteRawMF( *F[1], dm, filename, level_prefix, field_name+"y_cp", lev, plot_guards );
        WriteRawMF( *F[2], dm, filename, level_prefix, field_name+"z_cp", lev, plot_guards );
    }
}

/** \brief Write the coarse scalar multifab `F_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseScalar( const std::string field_name,
    const MultiFab* F_cp,
    const MultiFab* F_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards,
    const int icomp )
{
    int ng = 0;
    if (plot_guards) ng = F_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawMF( *F_fp, dm, filename, level_prefix, field_name+"_cp", lev, ng );
    } else {
        // Create an alias to the component `icomp` of F_cp
        MultiFab F_comp(*F_cp, amrex::make_alias, icomp, 1);
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = Interpolate::getInterpolatedScalar( F_comp, *F_fp, dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawMF( *F, dm, filename, level_prefix, field_name+"_cp", lev, plot_guards );
    }
}

void
FlushFormatPlotfile::WriteAllRawFields(
    const bool plot_raw_fields, const int nlevels, const std::string& plotfilename,
    const bool plot_raw_fields_guards, const bool plot_rho, bool plot_F) const
{
    if (!plot_raw_fields) return;
    auto & warpx = WarpX::GetInstance();
    for (int lev = 0; lev < nlevels; ++lev)
    {
        const std::unique_ptr<MultiFab> empty_ptr;
        const std::string raw_pltname = plotfilename + "/raw_fields";
        const DistributionMapping& dm = warpx.DistributionMap(lev);

        // Auxiliary patch

        WriteRawMF( warpx.getEfield(lev, 0), dm, raw_pltname, level_prefix, "Ex_aux", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getEfield(lev, 1), dm, raw_pltname, level_prefix, "Ey_aux", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getEfield(lev, 2), dm, raw_pltname, level_prefix, "Ez_aux", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield(lev, 0), dm, raw_pltname, level_prefix, "Bx_aux", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield(lev, 1), dm, raw_pltname, level_prefix, "By_aux", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield(lev, 2), dm, raw_pltname, level_prefix, "Bz_aux", lev, plot_raw_fields_guards);

        // fine patch
        WriteRawMF( warpx.getEfield_fp(lev, 0), dm, raw_pltname, level_prefix, "Ex_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getEfield_fp(lev, 1), dm, raw_pltname, level_prefix, "Ey_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getEfield_fp(lev, 2), dm, raw_pltname, level_prefix, "Ez_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getcurrent_fp(lev, 0), dm, raw_pltname, level_prefix, "jx_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getcurrent_fp(lev, 1), dm, raw_pltname, level_prefix, "jy_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getcurrent_fp(lev, 2), dm, raw_pltname, level_prefix, "jz_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield_fp(lev, 0), dm, raw_pltname, level_prefix, "Bx_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield_fp(lev, 1), dm, raw_pltname, level_prefix, "By_fp", lev, plot_raw_fields_guards);
        WriteRawMF( warpx.getBfield_fp(lev, 2), dm, raw_pltname, level_prefix, "Bz_fp", lev, plot_raw_fields_guards);
        if (plot_F) WriteRawMF( warpx.getF_fp(lev), dm, raw_pltname, level_prefix, "F_fp", lev, plot_raw_fields_guards);
        if (plot_rho) {
            // Use the component 1 of `rho_fp`, i.e. rho_new for time synchronization
            // If nComp > 1, this is the upper half of the list of components.
            MultiFab rho_new(warpx.getF_fp(lev), amrex::make_alias, warpx.getF_fp(lev).nComp()/2, warpx.getF_fp(lev).nComp()/2);
            WriteRawMF( rho_new, dm, raw_pltname, level_prefix, "rho_fp", lev, plot_raw_fields_guards);
        }

        // Coarse path
        if (lev > 0){
            WriteCoarseVector( "E",
                               warpx.get_pointer_Efield_cp(lev, 0), warpx.get_pointer_Efield_cp(lev, 1), warpx.get_pointer_Efield_cp(lev, 2),
                               warpx.get_pointer_Efield_fp(lev, 0), warpx.get_pointer_Efield_fp(lev, 1), warpx.get_pointer_Efield_fp(lev, 2),
                               dm, raw_pltname, level_prefix, lev, plot_raw_fields_guards);
            WriteCoarseVector( "B",
                               warpx.get_pointer_Bfield_cp(lev, 0), warpx.get_pointer_Bfield_cp(lev, 1), warpx.get_pointer_Bfield_cp(lev, 2),
                               warpx.get_pointer_Bfield_fp(lev, 0), warpx.get_pointer_Bfield_fp(lev, 1), warpx.get_pointer_Bfield_fp(lev, 2),
                               dm, raw_pltname, level_prefix, lev, plot_raw_fields_guards);
            WriteCoarseVector( "j",
                               warpx.get_pointer_current_cp(lev, 0), warpx.get_pointer_current_cp(lev, 1), warpx.get_pointer_current_cp(lev, 2),
                               warpx.get_pointer_current_fp(lev, 0), warpx.get_pointer_current_fp(lev, 1), warpx.get_pointer_current_fp(lev, 2),
                               dm, raw_pltname, level_prefix, lev, plot_raw_fields_guards);
        }
        if (plot_F) WriteCoarseScalar(
            "F", warpx.get_pointer_F_cp(lev), warpx.get_pointer_F_fp(lev),
            dm, raw_pltname, level_prefix, lev,
            plot_raw_fields_guards, 0);
        if (plot_rho) WriteCoarseScalar(
            "rho", warpx.get_pointer_rho_cp(lev), warpx.get_pointer_rho_fp(lev),
            dm, raw_pltname, level_prefix, lev,
            plot_raw_fields_guards, 1);
        // Use the component 1 of `rho_cp`, i.e. rho_new for time synchronization
    }
}
