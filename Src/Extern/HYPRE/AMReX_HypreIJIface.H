#ifndef AMREX_HYPREIJIFACE_H
#define AMREX_HYPREIJIFACE_H
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_parcsr_mv.h>

#include <memory>

namespace amrex {

class HypreIJIface
{
public:
    using HypreIntType = HYPRE_Int;
    using HypreRealType = HYPRE_Real;

    HypreIJIface(MPI_Comm comm, HypreIntType ilower, HypreIntType iupper,
                 int verbose);

    ~HypreIJIface();

    HypreIJIface (HypreIJIface const&) = delete;
    HypreIJIface (HypreIJIface &&) = delete;
    HypreIJIface& operator= (HypreIJIface const&) = delete;
    HypreIJIface& operator= (HypreIJIface &&) = delete;

    void parse_inputs(const std::string& prefix = "hypre");

    void solve(HypreRealType rel_tol, HypreRealType abs_tol, HypreIntType max_iter);

    //! IJ matrix instance
    [[nodiscard]] HYPRE_IJMatrix A() { return m_mat; }

    //! Right hand side IJ vector instance
    [[nodiscard]] HYPRE_IJVector b() { return m_rhs; }

    //! Solution IJ vector instance
    [[nodiscard]] HYPRE_IJVector x() { return m_sln; }

    //! Number of iterations taken by the solver to reach the desired tolerance
    [[nodiscard]] HypreIntType getNumIters() const { return m_num_iterations; }

    //! Final residual norm after a linear solve
    [[nodiscard]] HypreRealType getFinalResidualNorm() const { return m_final_res_norm; }

    [[nodiscard]] bool adjustSingularMatrix() const { return m_adjust_singular_matrix; }

private:
    void init_preconditioner(const std::string& prefix, const std::string& name);
    void init_solver(const std::string& prefix, const std::string& name);

    //! run_setup
    void run_hypre_setup();
    void run_hypre_solve();

    // Preconditioners
    void boomeramg_precond_configure(const std::string& prefix);
    void euclid_precond_configure(const std::string& prefix);

    // Solvers
    void boomeramg_solver_configure(const std::string& prefix);
    void gmres_solver_configure(const std::string& prefix);
    void cogmres_solver_configure(const std::string& prefix);
    void lgmres_solver_configure(const std::string& prefix);
    void flex_gmres_solver_configure(const std::string& prefix);
    void bicgstab_solver_configure(const std::string& prefix);
    void pcg_solver_configure(const std::string& prefix);
    void hybrid_solver_configure(const std::string& prefix);


    MPI_Comm m_comm{MPI_COMM_NULL};

    HYPRE_IJMatrix m_mat{nullptr};
    HYPRE_IJVector m_rhs{nullptr};
    HYPRE_IJVector m_sln{nullptr};

    HYPRE_ParCSRMatrix m_parA{nullptr};
    HYPRE_ParVector m_parRhs{nullptr};
    HYPRE_ParVector m_parSln{nullptr};

    HYPRE_Solver m_solver{nullptr};
    HYPRE_Solver m_precond{nullptr};

    HypreIntType (*m_solverDestroyPtr)(HYPRE_Solver){nullptr};
    HypreIntType (*m_solverSetupPtr)(
        HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector){nullptr};
    HypreIntType (*m_solverSolvePtr)(
        HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector){nullptr};
    HypreIntType (*m_solverPrecondPtr)(
        HYPRE_Solver,
        HYPRE_PtrToParSolverFcn,
        HYPRE_PtrToParSolverFcn,
        HYPRE_Solver){nullptr};

    HypreIntType (*m_precondDestroyPtr)(HYPRE_Solver){nullptr};
    HypreIntType (*m_precondSetupPtr)(
        HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector){nullptr};
    HypreIntType (*m_precondSolvePtr)(
        HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector){nullptr};

    HypreIntType (*m_solverSetTolPtr)(HYPRE_Solver, amrex::Real){nullptr};
    HypreIntType (*m_solverSetAbsTolPtr)(HYPRE_Solver, amrex::Real){nullptr};
    HypreIntType (*m_solverSetMaxIterPtr)(HYPRE_Solver, HypreIntType){nullptr};
    HypreIntType (*m_solverNumItersPtr)(HYPRE_Solver, HypreIntType*){nullptr};
    HypreIntType (*m_solverFinalResidualNormPtr)(HYPRE_Solver, amrex::Real*){nullptr};

    HypreIntType m_ilower{0};
    HypreIntType m_iupper{0};

    HypreRealType m_final_res_norm;
    HypreIntType m_num_iterations;

    std::string m_solver_name{"BoomerAMG"};
    std::string m_preconditioner_name{"none"};
    std::string m_file_prefix{"IJ"};

    //! Verbosity of the HYPRE solvers
    int m_verbose{0};

    unsigned int m_write_counter{0};

    //! Flag indicating whether a preconditioner has been set
    bool m_has_preconditioner{false};

    //! Flag indicating whether the solver/preconditioner has been setup
    bool m_need_setup{true};

    //! Flag indicating whether user has requested recomputation of preconditioner
    bool m_recompute_preconditioner{true};

    //! Should singular matrix be adjusted to fix solution on a node/cell?
    bool m_adjust_singular_matrix{false};

    //! Flag indicating whether to dump matrix files
    bool m_write_files{false};

    //! Flag indicating whether the files are overwritten on subsequent writes
    bool m_overwrite_files{true};
};

} // namespace amrex

#endif /* AMREX_HYPREIJIFACE_H */
