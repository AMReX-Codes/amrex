#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_RealVect.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#include <MultiParticleContainer.H>
#include <PML.H>
#include <BoostedFrameDiagnostic.H>
#include <BilinearFilter.H>
#include <NCIGodfreyFilter.H>

#ifdef WARPX_USE_PSATD
#include <fftw3.h>
#include <SpectralSolver.H>
#endif

#if defined(BL_USE_SENSEI_INSITU)
namespace amrex {
class AmrMeshInSituBridge;
}
#endif

enum struct DtType : int
{
    Full = 0,
    FirstHalf,
    SecondHalf
};

enum struct PatchType : int
{
    fine,
    coarse
};

class WarpX
    : public amrex::AmrCore
{
public:

    friend class PML;

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    static std::string Version ();
    static std::string PicsarVersion ();

    int Verbose () const { return verbose; }

    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    static void shiftMF(amrex::MultiFab& mf, const amrex::Geometry& geom, int num_shift, int dir);

    static void GotoNextLine (std::istream& is);

    // External fields
    static amrex::Vector<amrex::Real> B_external;

    // Algorithms
    static long use_picsar_deposition;
    static long current_deposition_algo;
    static long charge_deposition_algo;
    static long field_gathering_algo;
    static long particle_pusher_algo;
    static int maxwell_fdtd_solver_id;

    // Interpolation order
    static long nox;
    static long noy;
    static long noz;

    static bool use_fdtd_nci_corr;
    static int  l_lower_order_in_v;

    static bool use_filter;
    static bool serialize_ics;

    // Back transformation diagnostic
    static bool do_boosted_frame_diagnostic;
    static std::string lab_data_directory;
    static int  num_snapshots_lab;
    static amrex::Real dt_snapshots_lab;
    static bool do_boosted_frame_fields;
    static bool do_boosted_frame_particles;

    // Boosted frame parameters
    static amrex::Real gamma_boost;
    static amrex::Real beta_boost;
    static amrex::Vector<int> boost_direction;
    static amrex::Real zmax_plasma_to_compute_max_step;
    static int do_compute_max_step_from_zmax;

    static bool do_dynamic_scheduling;
    static bool refine_plasma;

    static int sort_int;

    // buffers
    static int n_field_gather_buffer;
    static int n_current_deposition_buffer;

    // do nodal
    static int do_nodal;

    const amrex::MultiFab& getcurrent (int lev, int direction) {return *current_fp[lev][direction];}
    const amrex::MultiFab& getEfield  (int lev, int direction) {return *Efield_aux[lev][direction];}
    const amrex::MultiFab& getBfield  (int lev, int direction) {return *Bfield_aux[lev][direction];}

    const amrex::MultiFab& getcurrent_cp (int lev, int direction) {return *current_cp[lev][direction];}
    const amrex::MultiFab& getEfield_cp  (int lev, int direction) {return  *Efield_cp[lev][direction];}
    const amrex::MultiFab& getBfield_cp  (int lev, int direction) {return  *Bfield_cp[lev][direction];}

    const amrex::MultiFab& getcurrent_fp (int lev, int direction) {return *current_fp[lev][direction];}
    const amrex::MultiFab& getEfield_fp  (int lev, int direction) {return *Efield_fp[lev][direction];}
    const amrex::MultiFab& getBfield_fp  (int lev, int direction) {return *Bfield_fp[lev][direction];}

    static amrex::MultiFab* getCosts (int lev) {
        if (m_instance) {
            return m_instance->costs[lev].get();
        } else {
            return nullptr;
        }
    }

    static amrex::IntVect filter_npass_each_dir;
    BilinearFilter bilinear_filter;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_exeybz;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_bxbyez;
    
    static int num_mirrors;
    amrex::Vector<amrex::Real> mirror_z;
    amrex::Vector<amrex::Real> mirror_z_width;
    amrex::Vector<int> mirror_z_npoints;
    
    void applyMirrors(amrex::Real time);

    void ComputeDt ();
    // Compute max_step automatically for simulations in a boosted frame.
    void computeMaxStepBoostAccelerator(amrex::Geometry geom);
    int  MoveWindow (bool move_j);
    void UpdatePlasmaInjectionPosition (amrex::Real dt);

    void ResetProbDomain (const amrex::RealBox& rb);

    void EvolveE (         amrex::Real dt);
    void EvolveE (int lev, amrex::Real dt);
    void EvolveB (         amrex::Real dt);
    void EvolveB (int lev, amrex::Real dt);
    void EvolveF (         amrex::Real dt, DtType dt_type);
    void EvolveF (int lev, amrex::Real dt, DtType dt_type);
    void EvolveB (int lev, PatchType patch_type, amrex::Real dt);
    void EvolveE (int lev, PatchType patch_type, amrex::Real dt);
    void EvolveF (int lev, PatchType patch_type, amrex::Real dt, DtType dt_type);

    void DampPML ();
    void DampPML (int lev);
    void DampPML (int lev, PatchType patch_type);

    void PushParticlesandDepose (int lev, amrex::Real cur_time);
    void PushParticlesandDepose (         amrex::Real cur_time);

    // This function does aux(lev) = fp(lev) + I(aux(lev-1)-cp(lev)).
    // Caller must make sure fp and cp have ghost cells filled.
    void UpdateAuxilaryData ();

    // Fill boundary cells including coarse/fine boundaries
    void FillBoundaryB ();
    void FillBoundaryE ();
    void FillBoundaryF ();
    void FillBoundaryE (int lev);
    void FillBoundaryB (int lev);
    void FillBoundaryF (int lev);

    void SyncCurrent ();
    void SyncRho ();


    int getistep (int lev) const {return istep[lev];}
    void setistep (int lev, int ii) {istep[lev] = ii;}
    amrex::Real gett_new (int lev) const {return t_new[lev];}
    void sett_new (int lev, amrex::Real time) {t_new[lev] = time;}
    amrex::Real getdt (int lev) const {return dt[lev];}

    int maxStep () const {return max_step;}
    amrex::Real stopTime () const {return stop_time;}

    int checkInt () const {return check_int;}
    int plotInt () const {return plot_int;}

    void WriteCheckPointFile () const;
    void WritePlotFile () const;
    void UpdateInSitu () const;
    void AverageAndPackFields( amrex::Vector<std::string>& varnames,
        amrex::Vector<amrex::MultiFab>& mf_avg, const int ngrow) const;

    void WritePlotFileES(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                         const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E);

    static std::array<amrex::Real,3> CellSize (int lev);
    static amrex::RealBox getRealBox(const amrex::Box& bx, int lev);
    static std::array<amrex::Real,3> LowerCorner (const amrex::Box& bx, int lev);
    static std::array<amrex::Real,3> UpperCorner (const amrex::Box& bx, int lev);

    static amrex::IntVect RefRatio (int lev);

    static const amrex::iMultiFab* CurrentBufferMasks (int lev);
    static const amrex::iMultiFab* GatherBufferMasks (int lev);

    static amrex::IntVect Bx_nodal_flag;
    static amrex::IntVect By_nodal_flag;
    static amrex::IntVect Bz_nodal_flag;

    static amrex::IntVect Ex_nodal_flag;
    static amrex::IntVect Ey_nodal_flag;
    static amrex::IntVect Ez_nodal_flag;

    static amrex::IntVect jx_nodal_flag;
    static amrex::IntVect jy_nodal_flag;
    static amrex::IntVect jz_nodal_flag;

    static int do_moving_window;
    static int moving_window_dir;

    // slice generation //
    void InitializeSliceMultiFabs ();
    void SliceGenerationForDiagnostics ();
    void WriteSlicePlotFile () const;
    void ClearSliceMultiFabs ();

    // these should be private, but can't due to Cuda limitations
    static void ComputeDivB (amrex::MultiFab& divB, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx);

    static void ComputeDivB (amrex::MultiFab& divB, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx, int ngrow);

    static void ComputeDivE (amrex::MultiFab& divE, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx);

    static void ComputeDivE (amrex::MultiFab& divE, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx, int ngrow);

protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) final;

    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping.  Only used during initialization.  Called
    //! by AmrCoreInitFromScratch.
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) final;

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) final
        { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) final;

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) final;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    ///
    /// Advance the simulation by numsteps steps, electromagnetic case.
    ///
    void EvolveEM(int numsteps);

    void FillBoundaryB (int lev, PatchType patch_type);
    void FillBoundaryE (int lev, PatchType patch_type);
    void FillBoundaryF (int lev, PatchType patch_type);

    void OneStep_nosub (amrex::Real t);
    void OneStep_sub1 (amrex::Real t);

    void RestrictCurrentFromFineToCoarsePatch (int lev);
    void AddCurrentFromFineLevelandSumBoundary (int lev);
    void StoreCurrent (int lev);
    void RestoreCurrent (int lev);
    void ApplyFilterandSumBoundaryJ (int lev, PatchType patch_type);
    void NodalSyncJ (int lev, PatchType patch_type);

    void RestrictRhoFromFineToCoarsePatch (int lev);
    void ApplyFilterandSumBoundaryRho (int lev, PatchType patch_type, int icomp, int ncomp);
    void AddRhoFromFineLevelandSumBoundary (int lev, int icomp, int ncomp);
    void NodalSyncRho (int lev, PatchType patch_type, int icomp, int ncomp);

#ifdef WARPX_DO_ELECTROSTATIC
    ///
    /// Advance the simulation by numsteps steps, electrostatic case.
    ///
    void EvolveES(int numsteps);

    ///
    /// Compute the gravitational potential from rho by solving Poisson's equation.
    /// Both rho and phi are assumed to be node-centered. This method is only used
    /// in electrostatic mode.
    ///
    void computePhi(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const;

    ///
    /// Compute the electric field in each direction by computing the gradient
    /// the potential phi using 2nd order centered differences. Both rho and phi
    /// are assumed to be node-centered. This method is only used in electrostatic mode.
    ///
    void computeE(amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                  const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const;

    //
    // This stuff is needed by the nodal multigrid solver when running in
    // electrostatic mode.
    //
    void zeroOutBoundary(amrex::MultiFab& input_data, amrex::MultiFab& bndry_data,
                         const amrex::FabArray<amrex::BaseFab<int> >& mask) const;

    void sumFineToCrseNodal(const amrex::MultiFab& fine, amrex::MultiFab& crse,
                            const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

    void fixRHSForSolve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rhs,
                        const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) const ;

    void getLevelMasks(amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks,
                       const int nnodes = 1);

    // used to zero out fine level data on points shared with the coarse grid
    // in electrostatic mode
    amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > > masks;

    // used to gather the field from the coarse level in electrostatic mode.
    amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > > gather_masks;
#endif // WARPX_DO_ELECTROSTATIC

    void ReadParameters ();

    void InitFromScratch ();

    void AllocLevelData (int lev, const amrex::BoxArray& new_grids,
                         const amrex::DistributionMapping& new_dmap);
    void InitLevelData (int lev, amrex::Real time);

    void InitFromCheckpoint ();
    void PostRestart ();

    void InitOpenbc ();

    void InitPML ();
    void ComputePMLFactors ();

    void InitFilter ();

    void InitDiagnostics ();

    void InitNCICorrector ();

    void WriteWarpXHeader(const std::string& name) const;
    void WriteJobInfo (const std::string& dir) const;

    std::unique_ptr<amrex::MultiFab> GetCellCenteredData();

    std::array<std::unique_ptr<amrex::MultiFab>, 3> getInterpolatedE(int lev) const;

    std::array<std::unique_ptr<amrex::MultiFab>, 3> getInterpolatedB(int lev) const;

    void SyncCurrent (const std::array<const amrex::MultiFab*,3>& fine,
                      const std::array<      amrex::MultiFab*,3>& crse,
                      int ref_ratio);

    void SyncRho (const amrex::MultiFab& fine, amrex::MultiFab& crse, int ref_ratio);

    void ExchangeWithPmlB (int lev);
    void ExchangeWithPmlE (int lev);
    void ExchangeWithPmlF (int lev);

    void LoadBalance ();

    void BuildBufferMasks ();
    const amrex::iMultiFab* getCurrentBufferMasks (int lev) const {
        return current_buffer_masks[lev].get();
    }
    const amrex::iMultiFab* getGatherBufferMasks (int lev) const {
        return gather_buffer_masks[lev].get();
    }

    void AllocLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                        const amrex::IntVect& ngE, const amrex::IntVect& ngJ,
                        const amrex::IntVect& ngRho, int ngF);

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;

    // Boosted Frame Diagnostics
    std::unique_ptr<BoostedFrameDiagnostic> myBFD;

    //
    // Fields: First array for level, second for direction
    //

    // Full solution
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_aux;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_aux;

    // Fine patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_fp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_fp;

    // store fine patch
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_store;

    // Coarse patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_cp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cp;

    // Copy of the coarse aux
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cax;
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cax;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > current_buffer_masks;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > gather_buffer_masks;

    // If charge/current deposition buffers are used
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_buf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf;

    amrex::Vector<std::array<std::unique_ptr<amrex::iMultiFab>, 3 > > current_fp_owner_masks;
    amrex::Vector<std::array<std::unique_ptr<amrex::iMultiFab>, 3 > > current_cp_owner_masks;

    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > rho_fp_owner_masks;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > rho_cp_owner_masks;

    // div E cleaning
    int do_dive_cleaning = 0;

    // PML
    int do_pml = 1;
    int pml_ncell = 10;
    int pml_delta = 10;
    amrex::Vector<std::unique_ptr<PML> > pml;

    amrex::Real moving_window_x = std::numeric_limits<amrex::Real>::max();
    amrex::Real moving_window_v = std::numeric_limits<amrex::Real>::max();
    amrex::Real current_injection_position = 0;

    // Plasma injection parameters
    int warpx_do_continuous_injection = 0;
    int num_injected_species = -1;
    amrex::Vector<int> injected_plasma_species;
    
    int do_electrostatic = 0;
    int n_buffer = 4;
    amrex::Real const_dt = 0.5e-11;

    int load_balance_int = -1;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > costs;
    int load_balance_with_sfc = 0;
    amrex::Real load_balance_knapsack_factor = 1.24;

    // Other runtime parameters
    int verbose = 1;

    int do_subcycling = 0;

    int max_step   = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    int regrid_int = -1;

    amrex::Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"checkpoints/chk"};
    std::string plot_file {"diags/plotfiles/plt"};
    std::string slice_plot_file {"diags/slice_plotfiles/plt"};
    int check_int = -1;
    int plot_int = -1;

#ifdef WARPX_USE_OPENPMD
    bool dump_plotfiles = false;
    bool dump_openpmd = true;
#else
    bool dump_plotfiles = true;
    bool dump_openpmd = false;
#endif
    bool plot_E_field       = true;
    bool plot_B_field       = true;
    bool plot_J_field       = true;
    bool plot_part_per_cell = true;
    bool plot_part_per_grid = false;
    bool plot_part_per_proc = false;
    bool plot_proc_number   = false;
    bool plot_dive          = false;
    bool plot_divb          = false;
    bool plot_rho           = false;
    bool plot_F             = false;
    bool plot_finepatch     = false;
    bool plot_crsepatch     = false;
    bool plot_raw_fields    = false;
    bool plot_raw_fields_guards = false;
    int plot_coarsening_ratio = 1;

    amrex::VisMF::Header::Version checkpoint_headerversion = amrex::VisMF::Header::NoFabHeader_v1;
//    amrex::VisMF::Header::Version plotfile_headerversion   = amrex::VisMF::Header::NoFabHeader_v1;
    amrex::VisMF::Header::Version plotfile_headerversion  = amrex::VisMF::Header::Version_v1;
    amrex::VisMF::Header::Version slice_plotfile_headerversion  = amrex::VisMF::Header::Version_v1;

    bool use_single_read = true;
    bool use_single_write = true;
    int mffile_nstreams = 4;
    int field_io_nfiles = 1024;
    int particle_io_nfiles = 1024;

    amrex::RealVect fine_tag_lo;
    amrex::RealVect fine_tag_hi;

    bool is_synchronized = true;

    //Slice Parameters
    int slice_max_grid_size;
    int slice_plot_int = -1;
    amrex::RealBox slice_realbox;
    amrex::IntVect slice_cr_ratio;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_slice;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_slice;

#ifdef WARPX_USE_PSATD
    // Store fields in real space on the dual grid (i.e. the grid for the FFT push of the fields)
    // This includes data for the FFT guard cells (between FFT groups)
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp_fft;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_fp_fft;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_fft;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_fp_fft;

    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cp_fft;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cp_fft;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_cp_fft;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_cp_fft;

public:
    // FFTData is a stuct containing a 22 pointers to auxiliary arrays
    // 1-11: padded arrays in real space (required by FFTW); 12-22: arrays in spectral space
    struct FFTData
    {
        static constexpr int N = 22;
        void* data[N] = { nullptr };

        ~FFTData () {
            for (int i = 0; i < N; ++i) { // The memory is allocated with fftw_alloc.
                fftw_free(data[i]);
                data[i] = nullptr;
            }
        }

        FFTData () = default;

        FFTData (FFTData && rhs) noexcept {
            for (int i = 0; i < N; ++i) {
                data[i] = rhs.data[i];
                rhs.data[i] = nullptr;
            }
        }

        FFTData (FFTData const&) = delete;
        void operator= (FFTData const&) = delete;
        void operator= (FFTData&&) = delete;
    };

private:

    amrex::Vector<std::unique_ptr<amrex::LayoutData<FFTData> > > dataptr_fp_fft;
    amrex::Vector<std::unique_ptr<amrex::LayoutData<FFTData> > > dataptr_cp_fft;

    amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_fp;
    amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_cp;

  // Domain decomposition containing the valid part of the dual grids (i.e. without FFT guard cells)
    amrex::Vector<amrex::BoxArray> ba_valid_fp_fft;
    amrex::Vector<amrex::BoxArray> ba_valid_cp_fft;

    amrex::Vector<amrex::Box> domain_fp_fft;  //  "global" domain for the group this process belongs to
    amrex::Vector<amrex::Box> domain_cp_fft;

    amrex::Vector<MPI_Comm> comm_fft;
    amrex::Vector<int>      color_fft;

    bool fft_hybrid_mpi_decomposition = false;
    int ngroups_fft = 4;
    int fftw_plan_measure = 1;
    int nox_fft = 16;
    int noy_fft = 16;
    int noz_fft = 16;

    void AllocLevelDataFFT (int lev);
    void InitLevelDataFFT (int lev, amrex::Real time);

    void InitFFTComm (int lev);
    void FFTDomainDecomposition (int lev, amrex::BoxArray& ba_fft, amrex::DistributionMapping& dm_fft,
                                 amrex::BoxArray& ba_valid, amrex::Box& domain_fft,
                                 const amrex::Box& domain);
    void InitFFTDataPlan (int lev);
    void FreeFFT (int lev);

    void EvolvePSATD (int numsteps);
    void PushPSATD (amrex::Real dt);
    void PushPSATD (int lev, amrex::Real dt);

#endif

#if defined(BL_USE_SENSEI_INSITU)
    amrex::AmrMeshInSituBridge *insitu_bridge;
#endif
    int insitu_int;
    int insitu_start;
    std::string insitu_config;
    int insitu_pin_mesh;
};

#endif
