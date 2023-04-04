#include <AMReX_HypreIJIface.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

namespace amrex {

namespace {

/** Helper object to parse HYPRE inputs and call API functions
 */
struct HypreOptParse
{
    //! Input file parser instance for the given namespace
    amrex::ParmParse pp;

    //! Hypre solver/preconditioner whose options are being set
    HYPRE_Solver solver;

    HypreOptParse (const std::string& prefix, HYPRE_Solver sinp)
        : pp(prefix), solver(sinp)
    {}

    template <typename F>
    void operator() (const std::string& key, F&& func)
    {
        if (pp.contains(key.c_str())) {
            int val;
            pp.query(key.c_str(), val);
            func(solver, val);
        }
    }

    template <typename F, typename T>
    void operator() (const std::string& key, F&& func, T default_val)
    {
        T val = default_val;
        pp.queryAdd(key.c_str(), val);
        func(solver, val);
    }

    template <typename F, typename T>
    void operator() (const std::string& key, F&& func, T default_val, int index)
    {
        T val = default_val;
        pp.queryAdd(key.c_str(), val);
        func(solver, val, index);
    }

    template <typename T, typename F>
    void set (const std::string& key, F&& func)
    {
        if (pp.contains(key.c_str())) {
            T val;
            pp.query(key.c_str(), val);
            func(solver, val);
        }
    }
};

} // namespace

HypreIJIface::HypreIJIface (
    MPI_Comm comm, HypreIntType ilower, HypreIntType iupper, int verbose)
    : m_comm(comm), m_ilower(ilower), m_iupper(iupper), m_verbose(verbose)
{
    HYPRE_IJMatrixCreate(
        m_comm, m_ilower, m_iupper, m_ilower, m_iupper, &m_mat);
    HYPRE_IJMatrixSetObjectType(m_mat, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(m_mat);

    HYPRE_IJVectorCreate(m_comm, m_ilower, m_iupper, &m_rhs);
    HYPRE_IJVectorSetObjectType(m_rhs, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(m_rhs);

    HYPRE_IJVectorCreate(m_comm, m_ilower, m_iupper, &m_sln);
    HYPRE_IJVectorSetObjectType(m_sln, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(m_sln);
}

HypreIJIface::~HypreIJIface ()
{
    HYPRE_IJMatrixDestroy(m_mat);
    m_mat = nullptr;
    HYPRE_IJVectorDestroy(m_rhs);
    m_rhs = nullptr;
    HYPRE_IJVectorDestroy(m_sln);
    m_sln = nullptr;

    if (m_solver != nullptr) {
        m_solverDestroyPtr(m_solver);
        m_solver = nullptr;
    }

    if (m_precond != nullptr) {
        m_precondDestroyPtr(m_precond);
        m_precond = nullptr;
    }
}

void HypreIJIface::run_hypre_setup ()
{
    if (m_need_setup || m_recompute_preconditioner) {
        BL_PROFILE("HypreIJIface::run_hypre_setup()");
        if (m_has_preconditioner)
            m_solverPrecondPtr(
                m_solver, m_precondSolvePtr, m_precondSetupPtr, m_precond);

        m_solverSetupPtr(m_solver, m_parA, m_parRhs, m_parSln);
        m_need_setup = false;
    }
}

void HypreIJIface::run_hypre_solve ()
{
    BL_PROFILE("HypreIJIface::run_hypre_solve()");
    m_solverSolvePtr(m_solver, m_parA, m_parRhs, m_parSln);
}

void HypreIJIface::solve (
    HypreRealType rel_tol, HypreRealType abs_tol, HypreIntType max_iter)
{
    // Assuming that Matrix/rhs etc. has been assembled by calling code

    HYPRE_IJMatrixGetObject(m_mat, (void**)&m_parA);
    HYPRE_IJVectorGetObject(m_rhs, (void**)&m_parRhs);
    HYPRE_IJVectorGetObject(m_sln, (void**)&m_parSln);

    if (m_write_files) {
        const std::string matfile = amrex::Concatenate(
            m_file_prefix + "_A", static_cast<int>(m_write_counter)) + ".out";
        const std::string rhsfile = amrex::Concatenate(
            m_file_prefix + "_b", static_cast<int>(m_write_counter)) + ".out";
        HYPRE_IJMatrixPrint(m_mat, matfile.c_str());
        HYPRE_IJVectorPrint(m_rhs, rhsfile.c_str());
    }

    m_solverSetTolPtr(m_solver, rel_tol);
    m_solverSetMaxIterPtr(m_solver, max_iter);
    if ((abs_tol > 0.0) && (m_solverSetAbsTolPtr != nullptr))
        m_solverSetAbsTolPtr(m_solver, abs_tol);

    // setup
    run_hypre_setup();

    // solve
    run_hypre_solve();

    // diagnostics
    m_solverNumItersPtr(m_solver, &m_num_iterations);
    m_solverFinalResidualNormPtr(m_solver, &m_final_res_norm);

    if (m_write_files) {
        const std::string slnfile = amrex::Concatenate(
            m_file_prefix + "_x", static_cast<int>(m_write_counter)) + ".out";
        HYPRE_IJVectorPrint(m_sln, slnfile.c_str());

        // Increment counter if the user has requested output of multiple solves
        if (!m_overwrite_files) ++m_write_counter;
    }

    if (m_verbose > 1)
        amrex::Print() << "HYPRE " << m_solver_name
                       << ": Num. iterations = " << m_num_iterations
                       << "; Relative residual = " << m_final_res_norm
                       << std::endl;
}

void HypreIJIface::parse_inputs (const std::string& prefix)
{
    amrex::ParmParse pp(prefix);

    pp.queryAdd("hypre_solver", m_solver_name);
    pp.queryAdd("hypre_preconditioner", m_preconditioner_name);
    pp.queryAdd("recompute_preconditioner", m_recompute_preconditioner);
    pp.queryAdd("write_matrix_files", m_write_files);
    pp.queryAdd("overwrite_existing_matrix_files", m_overwrite_files);
    pp.queryAdd("adjust_singular_matrix", m_adjust_singular_matrix);

    if (m_verbose > 2)
        amrex::Print() << "HYPRE: solver = " << m_solver_name
                       << "; preconditioner = " << m_preconditioner_name
                       << std::endl;

    if (m_preconditioner_name == "none") {
        m_has_preconditioner = false;
    } else {
        m_has_preconditioner = true;
        init_preconditioner(prefix, m_preconditioner_name);
    }

    init_solver(prefix, m_solver_name);
}

void HypreIJIface::init_preconditioner (
    const std::string& prefix, const std::string& name)
{
    if (name == "BoomerAMG") {
        boomeramg_precond_configure(prefix);
    } else if (name == "euclid") {
        euclid_precond_configure(prefix);
    } else {
        amrex::Abort("Invalid HYPRE preconditioner specified: " + name);
    }
}

void HypreIJIface::init_solver (
    const std::string& prefix, const std::string& name)
{
    if (name == "BoomerAMG") {
        boomeramg_solver_configure(prefix);
    } else if (name == "GMRES") {
        gmres_solver_configure(prefix);
    } else if (name == "COGMRES") {
        cogmres_solver_configure(prefix);
    } else if (name == "LGMRES") {
        lgmres_solver_configure(prefix);
    } else if (name == "FlexGMRES") {
        flex_gmres_solver_configure(prefix);
    } else if (name == "BiCGSTAB") {
        bicgstab_solver_configure(prefix);
    } else if (name == "PCG") {
        pcg_solver_configure(prefix);
    } else if (name == "Hybrid") {
        hybrid_solver_configure(prefix);
    }else {
        amrex::Abort("Invalid HYPRE solver specified: " + name);
    }
}

void HypreIJIface::boomeramg_precond_configure (const std::string& prefix)
{
    if (m_verbose > 2)
        amrex::Print() << "Creating BoomerAMG preconditioner" << std::endl;
    HYPRE_BoomerAMGCreate(&m_precond);

    // Setup the pointers
    m_precondDestroyPtr = &HYPRE_BoomerAMGDestroy;
    m_precondSetupPtr = &HYPRE_BoomerAMGSetup;
    m_precondSolvePtr = &HYPRE_BoomerAMGSolve;

    // Parse options
    HypreOptParse hpp(prefix, m_precond);
    hpp("bamg_verbose", HYPRE_BoomerAMGSetPrintLevel);
    hpp("bamg_logging", HYPRE_BoomerAMGSetLogging);

    hpp("bamg_max_iterations", HYPRE_BoomerAMGSetMaxIter, 1);
    hpp("bamg_precond_tolerance", HYPRE_BoomerAMGSetTol, 0.0);
    hpp("bamg_coarsen_type", HYPRE_BoomerAMGSetCoarsenType, 6);
    hpp("bamg_cycle_type", HYPRE_BoomerAMGSetCycleType, 1);
    hpp("bamg_relax_order", HYPRE_BoomerAMGSetRelaxOrder, 1);

    if (hpp.pp.contains("bamg_down_relax_type") && hpp.pp.contains("bamg_up_relax_type") && hpp.pp.contains("bamg_coarse_relax_type")) {
        hpp("bamg_down_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 1);
        hpp("bamg_up_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 2);
        hpp("bamg_coarse_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 3);
    } else {
        hpp("bamg_relax_type", HYPRE_BoomerAMGSetRelaxType, 6);
    }

    if (hpp.pp.contains("bamg_num_down_sweeps") && hpp.pp.contains("bamg_num_up_sweeps") && hpp.pp.contains("bamg_num_coarse_sweeps")) {
        hpp("bamg_num_down_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 2, 1);
        hpp("bamg_num_up_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 2, 2);
        hpp("bamg_num_coarse_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 1, 3);
    } else {
        hpp("bamg_num_sweeps", HYPRE_BoomerAMGSetNumSweeps, 2);
    }

    hpp("bamg_max_levels", HYPRE_BoomerAMGSetMaxLevels, 20);
    hpp("bamg_strong_threshold", HYPRE_BoomerAMGSetStrongThreshold,
        (AMREX_SPACEDIM == 3) ? 0.57 : 0.25);
    hpp("bamg_interp_type", HYPRE_BoomerAMGSetInterpType, 0);

    hpp("bamg_variant", HYPRE_BoomerAMGSetVariant);
    hpp("bamg_keep_transpose", HYPRE_BoomerAMGSetKeepTranspose);
    hpp("bamg_min_coarse_size", HYPRE_BoomerAMGSetMinCoarseSize);
    hpp("bamg_max_coarse_size", HYPRE_BoomerAMGSetMaxCoarseSize);
    hpp("bamg_pmax_elmts", HYPRE_BoomerAMGSetPMaxElmts);
    hpp("bamg_agg_num_levels", HYPRE_BoomerAMGSetAggNumLevels);
    hpp("bamg_agg_interp_type", HYPRE_BoomerAMGSetAggInterpType);
    hpp("bamg_agg_pmax_elmts", HYPRE_BoomerAMGSetAggPMaxElmts);
    hpp("bamg_trunc_factor", HYPRE_BoomerAMGSetTruncFactor, 0.1);
    hpp("bamg_set_restriction", HYPRE_BoomerAMGSetRestriction, 0);

    if (hpp.pp.contains("bamg_non_galerkin_tol")) {
        hpp("bamg_non_galerkin_tol", HYPRE_BoomerAMGSetNonGalerkinTol);

        if (hpp.pp.contains("bamg_non_galerkin_level_tols")) {
            std::vector<int> levels;
            std::vector<amrex::Real> tols;
            hpp.pp.getarr("bamg_non_galerkin_level_levels", levels);
            hpp.pp.getarr("bamg_non_galerkin_level_tols", tols);

            if (levels.size() != tols.size())
                amrex::Abort(
                    "HypreIJIface: Invalid sizes for non-Galerkin level "
                    "tolerances");

            for (size_t i = 0; i < levels.size(); ++i)
                HYPRE_BoomerAMGSetLevelNonGalerkinTol(
                    m_precond, tols[i], levels[i]);
        }
    }

    if (hpp.pp.contains("bamg_smooth_type")) {
        int smooth_type;
        hpp.pp.get("bamg_smooth_type", smooth_type);

        hpp("bamg_smooth_type", HYPRE_BoomerAMGSetSmoothType);

#if defined(HYPRE_RELEASE_NUMBER) && (HYPRE_RELEASE_NUMBER >= 22100)
        // Process ILU smoother parameters
        if (smooth_type == 5) { // ParILUK
            hpp("bamg_smooth_num_sweeps", HYPRE_BoomerAMGSetSmoothNumSweeps);
            hpp("bamg_smooth_num_levels", HYPRE_BoomerAMGSetSmoothNumLevels);
            hpp("bamg_ilu_type", HYPRE_BoomerAMGSetILUType);
            hpp("bamg_ilu_level", HYPRE_BoomerAMGSetILULevel);
            hpp("bamg_ilu_max_iter", HYPRE_BoomerAMGSetILUMaxIter);
        }
        else if (smooth_type == 7) { // Pilut
            hpp("bamg_smooth_num_sweeps", HYPRE_BoomerAMGSetSmoothNumSweeps);
            hpp("bamg_smooth_num_levels", HYPRE_BoomerAMGSetSmoothNumLevels);
            hpp("bamg_ilu_max_iter", HYPRE_BoomerAMGSetILUMaxIter);
            hpp("bamg_ilu_max_row_nnz", HYPRE_BoomerAMGSetILUMaxRowNnz);
            hpp("bamg_ilu_drop_tol", HYPRE_BoomerAMGSetILUDroptol, 1.e-10);
        }
#endif

        // Process Euclid smoother parameters
        if (smooth_type == 9) {
            if (hpp.pp.contains("bamg_euclid_file")) {
                std::string euclid_file;
                hpp.pp.get("bamg_euclid_file", euclid_file);
                HYPRE_BoomerAMGSetEuclidFile(
                    m_precond, const_cast<char*>(euclid_file.c_str()));
            }
            hpp("bamg_smooth_num_levels", HYPRE_BoomerAMGSetSmoothNumLevels);
            hpp("bamg_smooth_num_sweeps", HYPRE_BoomerAMGSetSmoothNumSweeps);
        }
    }
}

void HypreIJIface::euclid_precond_configure (const std::string& prefix)
{
    HYPRE_EuclidCreate(m_comm, &m_precond);

    // Setup the pointers
    m_precondDestroyPtr = &HYPRE_EuclidDestroy;
    m_precondSetupPtr = &HYPRE_EuclidSetup;
    m_precondSolvePtr = &HYPRE_EuclidSolve;

    HypreOptParse hpp(prefix, m_precond);
    // for 3D problems set to 1, for 2D set to 4-8
    hpp("euclid_level", HYPRE_EuclidSetLevel, (AMREX_SPACEDIM == 3) ? 1 : 4);
    // 0 = PILU; 1 = Block Jacobi
    hpp("euclid_use_block_jacobi", HYPRE_EuclidSetBJ, 0);
    // Flag indicating whether to write out euclid stats
    hpp("euclid_stats", HYPRE_EuclidSetStats, 0);
    // Flag indicating whether to print out memory diagnostic
    hpp("euclid_mem", HYPRE_EuclidSetMem, 0);
}

void HypreIJIface::boomeramg_solver_configure (const std::string& prefix)
{
    if (m_has_preconditioner) {
        amrex::Warning(
            "HYPRE: Cannot use preconditioner with BoomerAMG solver");
        m_has_preconditioner = false;
    }

    HYPRE_BoomerAMGCreate(&m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_BoomerAMGDestroy;
    m_solverSetupPtr = &HYPRE_BoomerAMGSetup;
    m_solverPrecondPtr = nullptr;
    m_solverSolvePtr = &HYPRE_BoomerAMGSolve;

    m_solverSetTolPtr = &HYPRE_BoomerAMGSetTol;
    m_solverSetAbsTolPtr = nullptr;
    m_solverSetMaxIterPtr = &HYPRE_BoomerAMGSetMaxIter;
    m_solverNumItersPtr = &HYPRE_BoomerAMGGetNumIterations;
    m_solverFinalResidualNormPtr = &HYPRE_BoomerAMGGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_BoomerAMGSetPrintLevel);
    hpp("logging", HYPRE_BoomerAMGSetLogging);
    hpp("bamg_relax_order", HYPRE_BoomerAMGSetRelaxOrder, 1);

    if (hpp.pp.contains("bamg_down_relax_type") && hpp.pp.contains("bamg_up_relax_type") && hpp.pp.contains("bamg_coarse_relax_type")) {
        hpp("bamg_down_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 1);
        hpp("bamg_up_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 2);
        hpp("bamg_coarse_relax_type", HYPRE_BoomerAMGSetCycleRelaxType, 11, 3);
    } else {
        hpp("bamg_relax_type", HYPRE_BoomerAMGSetRelaxType, 6);
    }

    if (hpp.pp.contains("bamg_num_down_sweeps") && hpp.pp.contains("bamg_num_up_sweeps") && hpp.pp.contains("bamg_num_coarse_sweeps")) {
        hpp("bamg_num_down_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 2, 1);
        hpp("bamg_num_up_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 2, 2);
        hpp("bamg_num_coarse_sweeps", HYPRE_BoomerAMGSetCycleNumSweeps, 1, 3);
    } else {
        hpp("bamg_num_sweeps", HYPRE_BoomerAMGSetNumSweeps, 2);
    }

    hpp("bamg_strong_threshold", HYPRE_BoomerAMGSetStrongThreshold,
        (AMREX_SPACEDIM == 3) ? 0.57 : 0.25);
    hpp("bamg_coarsen_type", HYPRE_BoomerAMGSetCoarsenType);
    hpp("bamg_cycle_type", HYPRE_BoomerAMGSetCycleType);
    hpp("bamg_max_levels", HYPRE_BoomerAMGSetMaxLevels);

    bool use_old_default = true;
    hpp.pp.queryAdd("bamg_use_old_default", use_old_default);
    if (use_old_default)
        HYPRE_BoomerAMGSetOldDefault(m_solver);
}

void HypreIJIface::gmres_solver_configure (const std::string& prefix)
{
    if (m_verbose > 2)
        amrex::Print() << "Creating GMRES solver" << std::endl;
    HYPRE_ParCSRGMRESCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRGMRESDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRGMRESSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRGMRESSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRGMRESSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRGMRESSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRGMRESSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRGMRESSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRGMRESGetNumIterations;
    m_solverFinalResidualNormPtr =
        &HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRGMRESSetPrintLevel);
    hpp("logging", HYPRE_ParCSRGMRESSetLogging);

    hpp("num_krylov", HYPRE_ParCSRGMRESSetKDim, 50);
    hpp("max_iterations", HYPRE_ParCSRGMRESSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRGMRESSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRGMRESSetAbsoluteTol);
}

void HypreIJIface::cogmres_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRCOGMRESCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRCOGMRESDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRCOGMRESSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRCOGMRESSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRCOGMRESSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRCOGMRESSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRCOGMRESSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRCOGMRESSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRCOGMRESGetNumIterations;
    m_solverFinalResidualNormPtr =
        &HYPRE_ParCSRCOGMRESGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRCOGMRESSetPrintLevel);
    hpp("logging", HYPRE_ParCSRCOGMRESSetLogging);

    hpp("num_krylov", HYPRE_ParCSRCOGMRESSetKDim, 50);
    hpp("max_iterations", HYPRE_ParCSRCOGMRESSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRCOGMRESSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRCOGMRESSetAbsoluteTol);
}

void HypreIJIface::lgmres_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRLGMRESCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRLGMRESDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRLGMRESSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRLGMRESSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRLGMRESSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRLGMRESSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRLGMRESSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRLGMRESSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRLGMRESGetNumIterations;
    m_solverFinalResidualNormPtr =
        &HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRLGMRESSetPrintLevel);
    hpp("logging", HYPRE_ParCSRLGMRESSetLogging);

    hpp("num_krylov", HYPRE_ParCSRLGMRESSetKDim, 50);
    hpp("max_iterations", HYPRE_ParCSRLGMRESSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRLGMRESSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRLGMRESSetAbsoluteTol);
}

void HypreIJIface::flex_gmres_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRFlexGMRESCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRFlexGMRESDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRFlexGMRESSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRFlexGMRESSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRFlexGMRESSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRFlexGMRESSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRFlexGMRESSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRFlexGMRESSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRFlexGMRESGetNumIterations;
    m_solverFinalResidualNormPtr =
        &HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRFlexGMRESSetPrintLevel);
    hpp("logging", HYPRE_ParCSRFlexGMRESSetLogging);

    hpp("num_krylov", HYPRE_ParCSRFlexGMRESSetKDim, 50);
    hpp("max_iterations", HYPRE_ParCSRFlexGMRESSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRFlexGMRESSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRFlexGMRESSetAbsoluteTol);
}

void HypreIJIface::bicgstab_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRBiCGSTABCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRBiCGSTABDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRBiCGSTABSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRBiCGSTABSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRBiCGSTABSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRBiCGSTABSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRBiCGSTABSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRBiCGSTABSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRBiCGSTABGetNumIterations;
    m_solverFinalResidualNormPtr =
        &HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRBiCGSTABSetPrintLevel);
    hpp("logging", HYPRE_ParCSRBiCGSTABSetLogging);

    hpp("max_iterations", HYPRE_ParCSRBiCGSTABSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRBiCGSTABSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRBiCGSTABSetAbsoluteTol);
}

void HypreIJIface::pcg_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRPCGCreate(m_comm, &m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRPCGDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRPCGSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRPCGSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRPCGSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRPCGSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRPCGSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRPCGSetMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRPCGGetNumIterations;
    m_solverFinalResidualNormPtr = &HYPRE_ParCSRPCGGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRPCGSetPrintLevel);
    hpp("logging", HYPRE_ParCSRPCGSetLogging);

    hpp("max_iterations", HYPRE_ParCSRPCGSetMaxIter, 200);
    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRPCGSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRPCGSetAbsoluteTol);
}

void HypreIJIface::hybrid_solver_configure (const std::string& prefix)
{
    HYPRE_ParCSRHybridCreate(&m_solver);

    // Setup pointers
    m_solverDestroyPtr = &HYPRE_ParCSRHybridDestroy;
    m_solverSetupPtr = &HYPRE_ParCSRHybridSetup;
    m_solverPrecondPtr = &HYPRE_ParCSRHybridSetPrecond;
    m_solverSolvePtr = &HYPRE_ParCSRHybridSolve;

    m_solverSetTolPtr = &HYPRE_ParCSRHybridSetTol;
    m_solverSetAbsTolPtr = &HYPRE_ParCSRHybridSetAbsoluteTol;
    m_solverSetMaxIterPtr = &HYPRE_ParCSRHybridSetDSCGMaxIter;
    m_solverNumItersPtr = &HYPRE_ParCSRHybridGetNumIterations;
    m_solverFinalResidualNormPtr = &HYPRE_ParCSRHybridGetFinalRelativeResidualNorm;

    // Parse options
    HypreOptParse hpp(prefix, m_solver);
    hpp("verbose", HYPRE_ParCSRHybridSetPrintLevel);
    hpp("logging", HYPRE_ParCSRHybridSetLogging);

    hpp.set<amrex::Real>("rtol", HYPRE_ParCSRHybridSetTol);
    hpp.set<amrex::Real>("atol", HYPRE_ParCSRHybridSetAbsoluteTol);

    hpp("num_krylov", HYPRE_ParCSRHybridSetKDim, 50);
    hpp("hybrid_dscg_max_iter", HYPRE_ParCSRHybridSetDSCGMaxIter);
    hpp("hybrid_pcg_max_iter", HYPRE_ParCSRHybridSetPCGMaxIter);
    hpp("hybrid_setup_type", HYPRE_ParCSRHybridSetSetupType);
    hpp("hybrid_solver_type", HYPRE_ParCSRHybridSetSolverType);
    hpp.set<HypreRealType>(
        "hybrid_set_strong_threshold", HYPRE_ParCSRHybridSetStrongThreshold);
    hpp("hybrid_recompute_residual", HYPRE_ParCSRHybridSetRecomputeResidual);
    hpp("hybrid_recompute_residual_period", HYPRE_ParCSRHybridSetRecomputeResidualP);
}

} // namespace amrex
