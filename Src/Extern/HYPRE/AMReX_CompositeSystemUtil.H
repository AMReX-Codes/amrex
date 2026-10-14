#ifndef AMREX_COMPOSITE_SYSTEM_UTIL_H_
#define AMREX_COMPOSITE_SYSTEM_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_AlgMG.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SpMatrix.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>

/**
 * \file AMReX_CompositeSystemUtil.H
 *
 * Helpers for the composite (all AMR levels) solvers such as
 * HypreMLABecLap: an owned-row CSR built from (row, column, value)
 * triplets, the AlgMG backend options, and Matrix Market output of the
 * assembled systems for analysis outside AMReX.
 */

namespace amrex::composite {

//! (row, column, value) of a matrix entry with global indices
struct Triplet { Long row; Long col; Real val; };

/**
 * \brief Rows owned by this rank in CSR form, columns sorted within each
 * row and duplicates summed.
 */
struct OwnedCSR
{
    Long row_begin = 0;          //!< first global row of this rank
    Long nrows = 0;              //!< number of rows owned by this rank
    Vector<Long> row_offset;     //!< size nrows+1
    Vector<Long> col_index;      //!< global column indices
    Vector<Real> values;

    [[nodiscard]] Long nnz () const { return values.size(); }

    //! Build from the triplets of this rank's rows (rows outside
    //! [row_begin, row_begin+nrows) are an error).
    void define (Long a_row_begin, Long a_nrows, Vector<Triplet> trips)
    {
        row_begin = a_row_begin;
        nrows = a_nrows;
        std::sort(trips.begin(), trips.end(), [] (Triplet const& a, Triplet const& b) {
            return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        });
        row_offset.assign(nrows+1, 0);
        col_index.clear();
        values.clear();
        col_index.reserve(trips.size());
        values.reserve(trips.size());
        Long prev_row = -1, prev_col = -1;
        for (auto const& t : trips) {
            AMREX_ALWAYS_ASSERT(t.row >= row_begin && t.row < row_begin+nrows);
            if (t.row == prev_row && t.col == prev_col) {
                values.back() += t.val;
            } else {
                col_index.push_back(t.col);
                values.push_back(t.val);
                ++row_offset[t.row-row_begin+1];
                prev_row = t.row;
                prev_col = t.col;
            }
        }
        for (Long r = 0; r < nrows; ++r) { row_offset[r+1] += row_offset[r]; }
    }

    //! Define an SpMatrix (device memory in GPU builds) with this CSR
    void toSpMatrix (SpMatrix<Real>& mat, Vector<Long> const& rank_offset) const
    {
        SpMatrix<Real>::csr_type csr;
        csr.mat.resize(values.size());
        csr.col_index.resize(col_index.size());
        csr.row_offset.resize(row_offset.size());
        csr.nnz = nnz();
        Gpu::copyAsync(Gpu::hostToDevice, values.begin(), values.end(), csr.mat.begin());
        Gpu::copyAsync(Gpu::hostToDevice, col_index.begin(), col_index.end(), csr.col_index.begin());
        Gpu::copyAsync(Gpu::hostToDevice, row_offset.begin(), row_offset.end(), csr.row_offset.begin());
        Gpu::streamSynchronize();
        mat.define(AlgPartition(rank_offset), std::move(csr), CsrSorted{true});
    }
};

//! Dumps written so far per dump prefix, over all solver objects of the run
inline std::map<std::string,int>& dump_counters ()
{
    static std::map<std::string,int> counter;
    return counter;
}

//! Whether a dump with this prefix is still pending (fewer than \p max written)
inline bool dump_pending (std::string const& prefix, int max)
{
    return max > 0 && dump_counters()[prefix] < max;
}

//! Index of the next dump for this prefix, or -1 once \p max have been written
inline int next_dump_index (std::string const& prefix, int max)
{
    if (!dump_pending(prefix, max)) { return -1; }
    return dump_counters()[prefix]++;
}

/**
 * \brief Gather the chunks of a distributed array on rank 0 (ranks in
 * order); other ranks return an empty vector.  The total must fit in int.
 */
template <typename T>
Vector<T> gather_on_root (Vector<T> const& v, MPI_Comm comm)
{
    const int nprocs = ParallelDescriptor::NProcs(comm);
    if (nprocs == 1) { return v; }
#ifdef AMREX_USE_MPI
    int nloc = static_cast<int>(v.size());
    Vector<int> counts(nprocs), displs(nprocs, 0);
    MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    Vector<T> g;
    if (ParallelDescriptor::MyProc(comm) == 0) {
        int ntot = 0;
        for (int p = 0; p < nprocs; ++p) { displs[p] = ntot; ntot += counts[p]; }
        g.resize(ntot);
    }
    auto const mt = ParallelDescriptor::Mpi_typemap<T>::type();
    MPI_Gatherv(v.data(), nloc, mt, g.data(), counts.data(), displs.data(), mt, 0, comm);
    return g;
#else
    return v;
#endif
}

/**
 * \brief Whether a dump of \p total values goes to one file written by rank
 * 0 (always on one rank).  Above \p gather_max, or above what MPI's int
 * counts can address, every rank writes its own part to `<file>.<rank>`.
 */
inline bool dump_single_file (Long total, Long gather_max, MPI_Comm comm)
{
    return (ParallelDescriptor::NProcs(comm) == 1) ||
        (total <= std::min(gather_max, Long(std::numeric_limits<int>::max())));
}

//! Open a dump file for writing, aborting on failure
inline std::ofstream open_dump (std::string const& name)
{
    std::ofstream ofs(name);
    if (!ofs) { amrex::Abort("Cannot open dump file " + name); }
    ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);
    return ofs;
}

//! Abort if writing a dump file failed
inline void check_dump (std::ofstream const& ofs, std::string const& name)
{
    if (!ofs) { amrex::Abort("Failed to write dump file " + name); }
}

/**
 * \brief Write the row map of a composite system: one line per row with the
 * 1-based row, AMR level, index (i j k) and the unknown's scaling.  \p info
 * holds (level,i,j,k) of this rank's rows, which start at the global row
 * \p row_begin.  Gathered on rank 0 or one file per rank, see
 * dump_single_file.
 */
inline void write_row_map (std::string const& file, std::string const& title,
                           Vector<std::array<Long,4>> const& info, Vector<Real> const& scaling,
                           Long row_begin, MPI_Comm comm, Long gather_max)
{
    const int nprocs = ParallelDescriptor::NProcs(comm);
    const int myproc = ParallelDescriptor::MyProc(comm);
    Long nglobal = info.size();
    ParallelAllReduce::Sum(nglobal, comm);

    Vector<Long> flat(info.size()*4);
    for (Long r = 0; r < info.size(); ++r) {
        for (int c = 0; c < 4; ++c) { flat[r*4+c] = info[r][c]; }
    }
    auto write = [&] (std::string const& name, std::string const& note, Long first_row,
                      Vector<Long> const& fl, Vector<Real> const& sc)
    {
        auto ofs = open_dump(name);
        ofs << "# " << title << "\n" << note << "# row level i j k scaling\n";
        for (Long r = 0; r < sc.size(); ++r) {
            ofs << first_row+r+1 << ' ' << fl[r*4] << ' ' << fl[r*4+1] << ' '
                << fl[r*4+2] << ' ' << fl[r*4+3] << ' ' << sc[r] << '\n';
        }
        check_dump(ofs, name);
    };
    if (dump_single_file(4*nglobal, gather_max, comm)) {
        auto gflat = gather_on_root(flat, comm);
        auto gscl = gather_on_root(scaling, comm);
        if (myproc == 0) { write(file, "", 0, gflat, gscl); }
    } else {
        write(file + "." + std::to_string(myproc),
              "# part " + std::to_string(myproc) + " of " + std::to_string(nprocs) +
              ": rows owned by this rank\n", row_begin, flat, scaling);
    }
}

/**
 * \brief Configure an AlgMG solver from `<prefix>.algmg.*`: verbose,
 * krylov (none|bicgstab|gmres|pcg; default \p krylov), smoother
 * (jacobi|l1_jacobi|chebyshev|l1_gauss_seidel), interp_type
 * (direct|mm_ext|mm_ext_i), bottom_solver (jacobi|bicgstab|gmres|direct),
 * strong_threshold, max_levels, max_coarse_size, agg_num_levels,
 * pmax_elmts, trunc_factor, pre_smooth, post_smooth, cheby_degree.
 */
inline void algmg_parse_options (AlgMG<Real>& mg, std::string const& prefix,
                                 AlgMGKrylovSolver krylov)
{
    ParmParse pp(prefix + ".algmg");
    int iv = 0;
    Real rv = 0;
    Long lv = 0;
    if (pp.query("verbose", iv)) { mg.setVerbose(iv); }
    pp.query_enum_case_insensitive("krylov", krylov);
    mg.setKrylovSolver(krylov);
    AlgMGSmoother sm{};
    if (pp.query_enum_case_insensitive("smoother", sm)) { mg.setSmoother(sm); }
    AlgMGInterpType it{};
    if (pp.query_enum_case_insensitive("interp_type", it)) { mg.setInterpType(it); }
    AlgMGBottomSolver bs{};
    if (pp.query_enum_case_insensitive("bottom_solver", bs)) { mg.setBottomSolver(bs); }
    if (pp.query("strong_threshold", rv)) { mg.setStrongThreshold(rv); }
    if (pp.query("max_levels", iv)) { mg.setMaxLevels(iv); }
    if (pp.query("max_coarse_size", lv)) { mg.setMaxCoarseSize(lv); }
    if (pp.query("agg_num_levels", iv)) { mg.setAggressiveNumLevels(iv); }
    if (pp.query("pmax_elmts", iv)) { mg.setPMaxElmts(iv); }
    if (pp.query("trunc_factor", rv)) { mg.setTruncFactor(rv); }
    if (pp.query("pre_smooth", iv)) { mg.setPreSmooth(iv); }
    if (pp.query("post_smooth", iv)) { mg.setPostSmooth(iv); }
    if (pp.query("cheby_degree", iv)) { mg.setChebyshevDegree(iv); }
}

/**
 * \brief Write a distributed matrix in Matrix Market coordinate format
 * (1-based indices, "general"), gathered on rank 0 or one file per rank
 * (see dump_single_file).  A per-rank file holds this rank's rows with the
 * global dimensions; concatenating the entry lines gives the matrix.
 */
inline void write_matrix_market (std::string const& file, OwnedCSR const& csr,
                                 Long nglobal, MPI_Comm comm, Long gather_max)
{
    const int nprocs = ParallelDescriptor::NProcs(comm);
    const int myproc = ParallelDescriptor::MyProc(comm);
    Long nnz_global = csr.nnz();
    ParallelAllReduce::Sum(nnz_global, comm);

    Vector<Long> rows(csr.nnz());
    for (Long r = 0; r < csr.nrows; ++r) {
        for (Long n = csr.row_offset[r]; n < csr.row_offset[r+1]; ++n) { rows[n] = csr.row_begin + r; }
    }
    auto write = [&] (std::string const& name, std::string const& note, Vector<Long> const& rr,
                      Vector<Long> const& cc, Vector<Real> const& vv)
    {
        auto ofs = open_dump(name);
        ofs << "%%MatrixMarket matrix coordinate real general\n" << note;
        ofs << nglobal << ' ' << nglobal << ' ' << vv.size() << '\n';
        for (Long n = 0; n < vv.size(); ++n) {
            ofs << rr[n]+1 << ' ' << cc[n]+1 << ' ' << vv[n] << '\n';
        }
        check_dump(ofs, name);
    };
    if (dump_single_file(nnz_global, gather_max, comm)) {
        auto grows = gather_on_root(rows, comm);
        auto gcols = gather_on_root(csr.col_index, comm);
        auto gvals = gather_on_root(csr.values, comm);
        if (myproc == 0) { write(file, "", grows, gcols, gvals); }
    } else {
        write(file + "." + std::to_string(myproc),
              "% part " + std::to_string(myproc) + " of " + std::to_string(nprocs) +
              ": rows owned by this rank\n", rows, csr.col_index, csr.values);
    }
}

/**
 * \brief Write a distributed vector (this rank's rows start at the global
 * row \p row_begin) as a Matrix Market dense array, gathered on rank 0 or
 * one file per rank (see dump_single_file); a per-rank file gives its
 * global row range in a comment.
 */
inline void write_vector_market (std::string const& file, Vector<Real> const& v,
                                 Long nglobal, Long row_begin, MPI_Comm comm, Long gather_max)
{
    const int nprocs = ParallelDescriptor::NProcs(comm);
    const int myproc = ParallelDescriptor::MyProc(comm);
    auto write = [&] (std::string const& name, std::string const& note, Vector<Real> const& x)
    {
        auto ofs = open_dump(name);
        ofs << "%%MatrixMarket matrix array real general\n" << note;
        ofs << x.size() << " 1\n";
        for (auto const& xi : x) { ofs << xi << '\n'; }
        check_dump(ofs, name);
    };
    if (dump_single_file(nglobal, gather_max, comm)) {
        auto g = gather_on_root(v, comm);
        if (myproc == 0) { write(file, "", g); }
    } else {
        write(file + "." + std::to_string(myproc),
              "% part " + std::to_string(myproc) + " of " + std::to_string(nprocs) +
              ": rows " + std::to_string(row_begin+1) + " to " +
              std::to_string(row_begin+v.size()) + " of " + std::to_string(nglobal) + "\n", v);
    }
}

}

#endif
