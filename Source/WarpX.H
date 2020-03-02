/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 *                     Axel Huebl, Burlen Loring, David Grote
 *                     Glenn Richardson, Junmin Gu, Luca Fedeli
 *                     Mathieu Lobet, Maxence Thevenet, Remi Lehe
 *                     Revathi Jambunathan, Weiqun Zhang, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_H_
#define WARPX_H_

#include "Evolve/WarpXDtType.H"
#include "Particles/MultiParticleContainer.H"
#include "BoundaryConditions/PML.H"
#include "Diagnostics/BackTransformedDiagnostic.H"
#include "Filter/BilinearFilter.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#ifdef WARPX_USE_PSATD_HYBRID
#   include "FieldSolver/PicsarHybridSpectralSolver/PicsarHybridFFTData.H"
#endif

#include "Parallelization/GuardCellManager.H"

#ifdef WARPX_USE_OPENPMD
#   include "Diagnostics/WarpXOpenPMD.H"
#endif
#include "Parser/WarpXParserWrapper.H"

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_RealVect.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#ifdef _OPENMP
#   include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <array>


#if defined(BL_USE_SENSEI_INSITU)
namespace amrex {
class AmrMeshInSituBridge;
}
#endif

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

    static std::string Version (); //! Version of WarpX executable
    static std::string PicsarVersion (); //! Version of PICSAR dependency

    int Verbose () const { return verbose; }

    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    static void shiftMF (amrex::MultiFab& mf, const amrex::Geometry& geom,
                         int num_shift, int dir, amrex::IntVect ng_extra,
                         amrex::Real external_field=0.0, bool useparser = false,
                         ParserWrapper<3> *field_parser=nullptr);

    static void GotoNextLine (std::istream& is);

    //! Author of an input file / simulation setup
    static std::string authors;

    // Initial field on the grid.
    static amrex::Vector<amrex::Real> E_external_grid;
    static amrex::Vector<amrex::Real> B_external_grid;

    // Initialization Type for External E and B on grid
    static std::string B_ext_grid_s;
    static std::string E_ext_grid_s;

    // Parser for B_external on the grid
    static std::string str_Bx_ext_grid_function;
    static std::string str_By_ext_grid_function;
    static std::string str_Bz_ext_grid_function;
    // Parser for E_external on the grid
    static std::string str_Ex_ext_grid_function;
    static std::string str_Ey_ext_grid_function;
    static std::string str_Ez_ext_grid_function;

    // ParserWrapper for B_external on the grid
    std::unique_ptr<ParserWrapper<3> > Bxfield_parser;
    std::unique_ptr<ParserWrapper<3> > Byfield_parser;
    std::unique_ptr<ParserWrapper<3> > Bzfield_parser;
    // ParserWrapper for E_external on the grid
    std::unique_ptr<ParserWrapper<3> > Exfield_parser;
    std::unique_ptr<ParserWrapper<3> > Eyfield_parser;
    std::unique_ptr<ParserWrapper<3> > Ezfield_parser;

    // Algorithms
    static long current_deposition_algo;
    static long charge_deposition_algo;
    static long field_gathering_algo;
    static long particle_pusher_algo;
    static int maxwell_fdtd_solver_id;
    static long load_balance_costs_update_algo;

    // div E cleaning
    static int do_dive_cleaning;

    // Interpolation order
    static long nox;
    static long noy;
    static long noz;

    // Number of modes for the RZ multimode version
    static long n_rz_azimuthal_modes;
    static long ncomps;

    static bool use_fdtd_nci_corr;
    static int  l_lower_order_in_v;

    static bool use_filter;
    static bool serialize_ics;

    // Back transformation diagnostic
    static bool do_back_transformed_diagnostics;
    static std::string lab_data_directory;
    static int  num_snapshots_lab;
    static amrex::Real dt_snapshots_lab;
    static bool do_back_transformed_fields;
    static bool do_back_transformed_particles;

    // Boosted frame parameters
    static amrex::Real gamma_boost;
    static amrex::Real beta_boost;
    static amrex::Vector<int> boost_direction;
    static amrex::Real zmax_plasma_to_compute_max_step;
    static int do_compute_max_step_from_zmax;

    static bool do_dynamic_scheduling;
    static bool refine_plasma;

    static int sort_int;
    static amrex::IntVect sort_bin_size;

    static int do_subcycling;

    static bool do_device_synchronize_before_profile;
    static bool safe_guard_cells;

    // buffers
    static int n_field_gather_buffer;       //! in number of cells from the edge (identical for each dimension)
    static int n_current_deposition_buffer; //! in number of cells from the edge (identical for each dimension)

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

    /** get low-high-low-high-... vector for each direction indicating if mother grid PMLs are enabled */
    std::vector<bool> getPMLdirections() const;

    static amrex::MultiFab* getCosts (int lev) {
        if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            if (m_instance) {
            return m_instance->costs[lev].get();
            } else {
                return nullptr;
            }
        } else {
            return nullptr;
        }
    }

    static amrex::Vector<amrex::Real>* getCostsHeuristic (int lev) {
        if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic) {
            if (m_instance) {
                return m_instance->costs_heuristic[lev].get();
            } else {
                return nullptr;
            }
        } else {
            return nullptr;
        }
    }

    static amrex::IntVect filter_npass_each_dir;
    BilinearFilter bilinear_filter;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_exeybz;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_bxbyez;

    amrex::Real time_of_last_gal_shift  = 0;
    amrex::Array<amrex::Real,3> v_galilean = {{0}};

    static int num_mirrors;
    amrex::Vector<amrex::Real> mirror_z;
    amrex::Vector<amrex::Real> mirror_z_width;
    amrex::Vector<int> mirror_z_npoints;

    /// object with all reduced diagnotics, similar to MultiParticleContainer for species.
    MultiReducedDiags* reduced_diags;

    void applyMirrors(amrex::Real time);

    void ComputeDt ();
    // Compute max_step automatically for simulations in a boosted frame.
    void computeMaxStepBoostAccelerator(const amrex::Geometry& geom);
    int  MoveWindow (bool move_j);

    /**
     * \brief
     * This function shifts the boundary of the grid by 'v_galilean*dt'.
     * In doding so, only positions attributes are changed while fields remain unchanged.
     */
    void ShiftGalileanBoundary ();
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

    /** \brief apply QED correction on electric field
     * \param dt vector of time steps (for all levels)
     */
    void Hybrid_QED_Push (         amrex::Vector<amrex::Real> dt);

    /** \brief apply QED correction on electric field for level lev
     * \param lev mesh refinement level
     * \param dt time step
     */
    void Hybrid_QED_Push (int lev, amrex::Real dt);

    /** \brief apply QED correction on electric field for level lev and patch type patch_type
     * \param lev mesh refinement level
     * \param dt patch_type which MR patch: PatchType::fine or PatchType::coarse
     * \param dt time step
     */
    void Hybrid_QED_Push (int lev, PatchType patch_type, amrex::Real dt);

    static amrex::Real quantum_xi_c2;

    /** \brief perform load balance; compute and communicate new `amrex::DistributionMapping`
     */
    void LoadBalance ();
    void LoadBalanceTimers ();
    void LoadBalanceHeuristic ();

#ifdef WARPX_DIM_RZ
    void ApplyInverseVolumeScalingToCurrentDensity(amrex::MultiFab* Jx,
                                                   amrex::MultiFab* Jy,
                                                   amrex::MultiFab* Jz,
                                                   int lev);

    void ApplyInverseVolumeScalingToChargeDensity(amrex::MultiFab* Rho,
                                                  int lev);
#endif

    void DampPML ();
    void DampPML (int lev);
    void DampPML (int lev, PatchType patch_type);

    void DampJPML ();
    void DampJPML (int lev);
    void DampJPML (int lev, PatchType patch_type);

    void CopyJPML ();

    PML* GetPML (int lev);

    void PushParticlesandDepose (int lev, amrex::Real cur_time, DtType a_dt_type=DtType::Full);
    void PushParticlesandDepose (         amrex::Real cur_time);

    // This function does aux(lev) = fp(lev) + I(aux(lev-1)-cp(lev)).
    // Caller must make sure fp and cp have ghost cells filled.
    void UpdateAuxilaryData ();
    void UpdateAuxilaryDataStagToNodal ();
    void UpdateAuxilaryDataSameType ();

    // Fill boundary cells including coarse/fine boundaries
    void FillBoundaryB   (amrex::IntVect ng, amrex::IntVect ng_extra_fine=amrex::IntVect::TheZeroVector());
    void FillBoundaryE   (amrex::IntVect ng, amrex::IntVect ng_extra_fine=amrex::IntVect::TheZeroVector());
    void FillBoundaryF   (amrex::IntVect ng);
    void FillBoundaryAux (amrex::IntVect ng);
    void FillBoundaryE   (int lev, amrex::IntVect ng, amrex::IntVect ng_extra_fine=amrex::IntVect::TheZeroVector());
    void FillBoundaryB   (int lev, amrex::IntVect ng, amrex::IntVect ng_extra_fine=amrex::IntVect::TheZeroVector());
    void FillBoundaryF   (int lev, amrex::IntVect ng);
    void FillBoundaryAux (int lev, amrex::IntVect ng);

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
    int openpmdInt () const {return openpmd_int;}

    void WriteCheckPointFile () const;
    void WriteOpenPMDFile () const;
    void WritePlotFile () const;
    void UpdateInSitu () const;
    void AverageAndPackFields( amrex::Vector<std::string>& varnames,
        amrex::Vector<amrex::MultiFab>& mf_avg, const int ngrow) const;
    void prepareFields( int const step, amrex::Vector<std::string>& varnames,
        amrex::Vector<amrex::MultiFab>& mf_avg,
        amrex::Vector<const amrex::MultiFab*>& output_mf,
        amrex::Vector<amrex::Geometry>& output_geom ) const;

    void WritePlotFileES(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                         const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E);

    static std::array<amrex::Real,3> CellSize (int lev);
    static amrex::RealBox getRealBox(const amrex::Box& bx, int lev);
    static std::array<amrex::Real,3> LowerCorner (const amrex::Box& bx,
                                                  std::array<amrex::Real,3> galilean_shift, int lev);
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
    static amrex::Real moving_window_v;

    // slice generation //
    void InitializeSliceMultiFabs ();
    void SliceGenerationForDiagnostics ();
    void WriteSlicePlotFile () const;
    void ClearSliceMultiFabs ();
    static int num_slice_snapshots_lab;
    static amrex::Real dt_slice_snapshots_lab;
    static amrex::Real particle_slice_width_lab;
    amrex::RealBox getSliceRealBox() const {return slice_realbox;}

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

    const amrex::IntVect getngE() const { return guard_cells.ng_alloc_EB; };
    const amrex::IntVect getngF() const { return guard_cells.ng_alloc_F; };

    void ComputeSpaceChargeField (bool const reset_fields);
    void AddSpaceChargeField (WarpXParticleContainer& pc);
    void computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                     std::array<amrex::Real, 3> const beta = {{0,0,0}},
                     amrex::Real const required_precision=1.e-11 ) const;
    void computeE (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                   const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<amrex::Real, 3> const beta = {{0,0,0}} ) const;
    void computeB (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& B,
                   const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<amrex::Real, 3> const beta = {{0,0,0}} ) const;

    /**
     * \brief
     * This function initializes the E and B fields on each level
     * using the parser and the user-defined function for the external fields.
     * The subroutine will parse the x_/y_z_external_grid_function and
     * then, the B or E multifab is initialized based on the (x,y,z) position
     * on the staggered yee-grid or cell-centered grid.
     */
    void InitializeExternalFieldsOnGridUsingParser (
         amrex::MultiFab *mfx, amrex::MultiFab *mfy, amrex::MultiFab *mfz,
         ParserWrapper<3> *xfield_parser, ParserWrapper<3> *yfield_parser,
         ParserWrapper<3> *zfield_parser, amrex::IntVect x_nodal_flag,
         amrex::IntVect y_nodal_flag, amrex::IntVect z_nodal_flag,
         const int lev);

    /** \brief add particle and cell contributions to multifabs in `costs_heuristic`
     */
    void ComputeCostsHeuristic ();

protected:

    /**
     * \brief
     *  This function initializes E, B, rho, and F, at all the levels
     *  of the multifab. rho and F are initialized with 0.
     *  The E and B fields are initialized using user-defined inputs.
     *  The initialization type is set using "B_ext_grid_init_style"
     *  and "E_ext_grid_init_style". The initialization style is set to "default"
     *  if not explicitly defined by the user, and the E and B fields are
     *  initialized with E_external_grid and B_external_grid, respectively, each with
     *  a default value of 0.
     *  If the initialization type for the E and B field is "constant",
     *  then, the E and B fields at all the levels are initialized with
     *  user-defined values for E_external_grid and B_external_grid.
     *  If the initialization type for B-field is set to
     *  "parse_B_ext_grid_function", then, the parser is used to read
     *  Bx_external_grid_function(x,y,z), By_external_grid_function(x,y,z),
     *  and Bz_external_grid_function(x,y,z).
     *  Similarly, if the E-field initialization type is set to
     *  "parse_E_ext_grid_function", then, the parser is used to read
     *  Ex_external_grid_function(x,y,z), Ey_external_grid_function(x,y,z),
     *  and Ex_external_grid_function(x,y,z). The parser for the E and B
     *  initialization assumes that the function has three independent
     *  variables, at max, namely, x, y, z. However, any number of constants
     *  can be used in the function used to define the E and B fields on the grid.
     */
    void InitLevelData (int lev, amrex::Real time);

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

    void FillBoundaryB (int lev, PatchType patch_type, amrex::IntVect ng);
    void FillBoundaryE (int lev, PatchType patch_type, amrex::IntVect ng);
    void FillBoundaryF (int lev, PatchType patch_type, amrex::IntVect ng);

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

    //
    // This stuff is needed by the nodal multigrid solver when running in
    // electrostatic mode.
    //
    void zeroOutBoundary(amrex::MultiFab& input_data, amrex::MultiFab& bndry_data,
                         const amrex::FabArray<amrex::BaseFab<int> >& mask) const;

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

    void InitFromCheckpoint ();
    void PostRestart ();

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

    void ExchangeWithPmlB (int lev);
    void ExchangeWithPmlE (int lev);
    void ExchangeWithPmlF (int lev);

    void BuildBufferMasks ();
    void BuildBufferMasksInBox ( const amrex::Box tbx, amrex::IArrayBox &buffer_mask,
                                 const amrex::IArrayBox &guard_mask, const int ng );
    const amrex::iMultiFab* getCurrentBufferMasks (int lev) const {
        return current_buffer_masks[lev].get();
    }
    const amrex::iMultiFab* getGatherBufferMasks (int lev) const {
        return gather_buffer_masks[lev].get();
    }

    void AllocLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                        const amrex::IntVect& ngE, const amrex::IntVect& ngJ,
                        const amrex::IntVect& ngRho, const amrex::IntVect& ngF,
                        const amrex::IntVect& ngextra, const bool aux_is_nodal);

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;

    // Boosted Frame Diagnostics
    std::unique_ptr<BackTransformedDiagnostic> myBFD;

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

    // PML
    int do_pml = 1;
    int pml_ncell = 10;
    int pml_delta = 10;
    int pml_has_particles = 0;
    int do_pml_j_damping = 0;
    int do_pml_in_domain = 0;
    amrex::IntVect do_pml_Lo = amrex::IntVect::TheUnitVector();
    amrex::IntVect do_pml_Hi = amrex::IntVect::TheUnitVector();
    amrex::Vector<std::unique_ptr<PML> > pml;

    amrex::Real moving_window_x = std::numeric_limits<amrex::Real>::max();
    amrex::Real current_injection_position = 0;

    // Plasma injection parameters
    int warpx_do_continuous_injection = 0;
    int num_injected_species = -1;
    amrex::Vector<int> injected_plasma_species;

    int do_electrostatic = 0;
    int n_buffer = 4;
    amrex::Real const_dt = 0.5e-11;

    // Load balancing
    /** Interval to perform load balance; `load_balance_int=5`, e.g., will load
     * balance during steps 5, 10, 15, etc. */
    int load_balance_int = -1;
    /** Collection of multifabs to keep track of weights, based on in-code timers,
     * for use in load balancing routines. */
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > costs;
    /** Collection of multifabs to keep track of weights, based on number of cells
     * and number of particles per box; for use in load balancing routines. */
    amrex::Vector<std::unique_ptr<amrex::Vector<amrex::Real> > > costs_heuristic;
    /** Load balance with 'space filling curve' strategy. */
    int load_balance_with_sfc = 0;
    /** Controls the maximum number of boxes that can be assigned to a rank during
     * load balance via the 'knapsack' strategy; e.g., if there are 4 boxes per rank,
     * `load_balance_knapsack_factor=2` limits the maximum number of boxes that can
     * be assigned to a rank to 8. */
    amrex::Real load_balance_knapsack_factor = 1.24;
    /** Weight factor for cells in updating `costs_heuristic`. */
    amrex::Real costs_heuristic_cells_wt = 0.1;
    /** Weight factor for particles in updating `costs_heuristic`. */
    amrex::Real costs_heuristic_particles_wt = 0.9;

    // Override sync every ? steps
    int override_sync_int = 1;

    // Other runtime parameters
    int verbose = 1;

    bool use_hybrid_QED = 0;

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

    std::string openpmd_backend {"default"};
    int openpmd_int = -1;
    bool openpmd_tspf = true; //!< one file per timestep (or one file for all steps)
#ifdef WARPX_USE_OPENPMD
    WarpXOpenPMDPlot* m_OpenPMDPlotWriter = nullptr;
#endif

    bool plot_rho           = false;
    bool plot_costs         = true;
    bool plot_finepatch     = false;
    bool plot_crsepatch     = false;
    bool plot_raw_fields    = false;
    bool plot_raw_fields_guards = false;
    amrex::Vector<std::string> fields_to_plot;
    int plot_coarsening_ratio = 1;

    amrex::VisMF::Header::Version checkpoint_headerversion = amrex::VisMF::Header::NoFabHeader_v1;
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

    guardCellManager guard_cells;

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

#ifdef WARPX_USE_PSATD_HYBRID
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
#endif

    bool fft_hybrid_mpi_decomposition = false;
    int nox_fft = 16;
    int noy_fft = 16;
    int noz_fft = 16;

#ifdef WARPX_USE_PSATD
private:
    void EvolvePSATD (int numsteps);
    void PushPSATD (amrex::Real dt);
    void PushPSATD_localFFT (int lev, amrex::Real dt);

    int ngroups_fft = 4;
    int fftw_plan_measure = 1;

    amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_fp;
    amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_cp;
#endif
    amrex::Vector<std::unique_ptr<FiniteDifferenceSolver>> m_fdtd_solver_fp;
    amrex::Vector<std::unique_ptr<FiniteDifferenceSolver>> m_fdtd_solver_cp;

#ifdef WARPX_USE_PSATD_HYBRID
private:
    amrex::Vector<std::unique_ptr<amrex::LayoutData<FFTData> > > dataptr_fp_fft;
    amrex::Vector<std::unique_ptr<amrex::LayoutData<FFTData> > > dataptr_cp_fft;

    // Domain decomposition containing the valid part of the dual grids (i.e. without FFT guard cells)
    amrex::Vector<amrex::BoxArray> ba_valid_fp_fft;
    amrex::Vector<amrex::BoxArray> ba_valid_cp_fft;

    amrex::Vector<amrex::Box> domain_fp_fft;  //  "global" domain for the group this process belongs to
    amrex::Vector<amrex::Box> domain_cp_fft;

    amrex::Vector<MPI_Comm> comm_fft;
    amrex::Vector<int>      color_fft;

    void AllocLevelDataFFT (int lev);
    void InitLevelDataFFT (int lev, amrex::Real time);

    void InitFFTComm (int lev);
    void FFTDomainDecomposition (int lev, amrex::BoxArray& ba_fft, amrex::DistributionMapping& dm_fft,
                                 amrex::BoxArray& ba_valid, amrex::Box& domain_fft,
                                 const amrex::Box& domain);
    void InitFFTDataPlan (int lev);
    void FreeFFT (int lev);
    void PushPSATD_hybridFFT (int lev, amrex::Real dt);
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
