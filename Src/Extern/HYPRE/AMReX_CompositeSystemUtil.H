#ifndef AMREX_COMPOSITE_SYSTEM_UTIL_H_
#define AMREX_COMPOSITE_SYSTEM_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_AlgMG.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SpMatrix.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <map>
#include <fstream>
#include <iomanip>
#include <limits>
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

    [[nodiscard]] Long nnz () const { return static_cast<Long>(values.size()); }

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
        Gpu::DeviceVector<Real> d_val(values.size());
        Gpu::DeviceVector<Long> d_col(col_index.size());
        Gpu::DeviceVector<Long> d_off(row_offset.size());
        Gpu::copyAsync(Gpu::hostToDevice, values.begin(), values.end(), d_val.begin());
        Gpu::copyAsync(Gpu::hostToDevice, col_index.begin(), col_index.end(), d_col.begin());
        Gpu::copyAsync(Gpu::hostToDevice, row_offset.begin(), row_offset.end(), d_off.begin());
        Gpu::streamSynchronize();
        mat.define(AlgPartition(rank_offset), d_val.data(), d_col.data(), nnz(),
                   d_off.data(), CsrSorted{true}, CsrValid{true});
    }
};

/**
 * \brief Send each triplet to the rank owning its row (rank p owns rows
 * [rank_offset[p], rank_offset[p+1])) and return this rank's triplets.
 */
inline Vector<Triplet> exchange_triplets (Vector<Triplet>&& trips, Vector<Long> const& rank_offset,
                                          MPI_Comm comm)
{
    int nprocs = 1, myproc = 0;
#ifdef AMREX_USE_MPI
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myproc);
#else
    amrex::ignore_unused(comm);
#endif
    if (nprocs == 1) { return std::move(trips); }
#ifdef AMREX_USE_MPI
    auto owner = [&] (Long row) {
        auto it = std::upper_bound(rank_offset.begin(), rank_offset.end(), row);
        return static_cast<int>(it - rank_offset.begin()) - 1;
    };
    Vector<Vector<Triplet>> sendbuf(nprocs);
    Vector<Triplet> local;
    for (auto const& t : trips) {
        const int p = owner(t.row);
        if (p == myproc) { local.push_back(t); } else { sendbuf[p].push_back(t); }
    }
    Vector<int> sendcnt(nprocs), recvcnt(nprocs), sdispl(nprocs), rdispl(nprocs);
    for (int p = 0; p < nprocs; ++p) {
        sendcnt[p] = static_cast<int>(sendbuf[p].size() * Long(sizeof(Triplet)));
    }
    MPI_Alltoall(sendcnt.data(), 1, MPI_INT, recvcnt.data(), 1, MPI_INT, comm);
    int stot = 0, rtot = 0;
    for (int p = 0; p < nprocs; ++p) {
        sdispl[p] = stot; stot += sendcnt[p];
        rdispl[p] = rtot; rtot += recvcnt[p];
    }
    Vector<char> sbuf(stot), rbuf(rtot);
    for (int p = 0; p < nprocs; ++p) {
        if (sendcnt[p] > 0) { std::memcpy(sbuf.data()+sdispl[p], sendbuf[p].data(), sendcnt[p]); }
    }
    MPI_Alltoallv(sbuf.data(), sendcnt.data(), sdispl.data(), MPI_CHAR,
                  rbuf.data(), recvcnt.data(), rdispl.data(), MPI_CHAR, comm);
    const Long nrecv = rtot / Long(sizeof(Triplet));
    const Long nlocal = local.size();
    local.resize(nlocal + nrecv);
    if (nrecv > 0) { std::memcpy(local.data()+nlocal, rbuf.data(), rtot); }
    return local;
#endif
}

/**
 * \brief Index of the next dump for this file prefix, counted over all
 * solver objects of the run (objects are often rebuilt between solves), or
 * -1 once \p max dumps have been written.
 */
inline int next_dump_index (std::string const& prefix, int max)
{
    static std::map<std::string,int> counter;
    int& n = counter[prefix];
    if (n >= max) { return -1; }
    return n++;
}

/**
 * \brief Write the row map of a composite system: one line per row with the
 * 1-based row, AMR level, cell or node index (i j k), and the unknown's
 * scaling (x = scaling * x').  \p info holds (level,i,j,k) of this rank's
 * rows in row order; ranks own consecutive row ranges in rank order.
 */
inline void write_row_map (std::string const& file, std::string const& title,
                           Vector<std::array<Long,4>> const& info, Vector<Real> const& scaling,
                           MPI_Comm comm)
{
    int nprocs = 1, myproc = 0;
#ifdef AMREX_USE_MPI
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myproc);
#else
    amrex::ignore_unused(comm);
#endif
    Vector<Long> flat(info.size()*4);
    for (Long r = 0; r < info.size(); ++r) {
        for (int c = 0; c < 4; ++c) { flat[r*4+c] = info[r][c]; }
    }
    Vector<Long> gflat;
    Vector<Real> gscl;
#ifdef AMREX_USE_MPI
    int n4 = static_cast<int>(flat.size()), n1 = static_cast<int>(scaling.size());
    Vector<int> c4(nprocs), d4(nprocs, 0), c1(nprocs), d1(nprocs, 0);
    MPI_Gather(&n4, 1, MPI_INT, c4.data(), 1, MPI_INT, 0, comm);
    MPI_Gather(&n1, 1, MPI_INT, c1.data(), 1, MPI_INT, 0, comm);
    if (myproc == 0) {
        int t4 = 0, t1 = 0;
        for (int p = 0; p < nprocs; ++p) { d4[p] = t4; t4 += c4[p]; d1[p] = t1; t1 += c1[p]; }
        gflat.resize(t4);
        gscl.resize(t1);
    }
    auto const lt = ParallelDescriptor::Mpi_typemap<Long>::type();
    auto const rt = ParallelDescriptor::Mpi_typemap<Real>::type();
    MPI_Gatherv(flat.data(), n4, lt, gflat.data(), c4.data(), d4.data(), lt, 0, comm);
    MPI_Gatherv(scaling.data(), n1, rt, gscl.data(), c1.data(), d1.data(), rt, 0, comm);
#else
    gflat = flat;
    gscl = scaling;
#endif
    if (myproc == 0) {
        std::ofstream ofs(file);
        ofs << "# " << title << "\n"
            << "# row level i j k scaling   (row is 1-based; x = scaling * x')\n";
        ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);
        for (Long r = 0; r < gscl.size(); ++r) {
            ofs << r+1 << ' ' << gflat[r*4] << ' ' << gflat[r*4+1] << ' '
                << gflat[r*4+2] << ' ' << gflat[r*4+3] << ' ' << gscl[r] << '\n';
        }
    }
}

/**
 * \brief Configure an AlgMG solver from `<prefix>.algmg.*`: verbose,
 * krylov (none|bicgstab|gmres|pcg), smoother
 * (jacobi|l1_jacobi|chebyshev|l1_gauss_seidel), interp_type
 * (direct|mm_ext|mm_ext_i), bottom_solver (jacobi|bicgstab|gmres|direct),
 * strong_threshold, max_levels, max_coarse_size, agg_num_levels,
 * pmax_elmts, trunc_factor, pre_smooth, post_smooth, cheby_degree.
 */
inline void algmg_parse_options (AlgMG<Real>& mg, std::string const& prefix)
{
    ParmParse pp(prefix + ".algmg");
    int iv = 0;
    Real rv = 0;
    Long lv = 0;
    if (pp.query("verbose", iv)) { mg.setVerbose(iv); }
    AlgMGKrylovSolver krylov{};
    if (pp.query_enum_case_insensitive("krylov", krylov)) { mg.setKrylovSolver(krylov); }
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
 * (1-based indices, "general").
 *
 * If the global number of entries is at most \p gather_max, rank 0 writes
 * the whole matrix to \p file. Otherwise every rank writes its rows to
 * `file.<rank>`, each a complete Matrix Market file with the global
 * dimensions and only that rank's entries (concatenate the entry lines
 * to obtain the whole matrix).
 */
inline void write_matrix_market (std::string const& file, OwnedCSR const& csr,
                                 Long nglobal, MPI_Comm comm, Long gather_max)
{
    int nprocs = 1, myproc = 0;
#ifdef AMREX_USE_MPI
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myproc);
#else
    amrex::ignore_unused(comm);
#endif
    Long nnz_global = csr.nnz();
#ifdef AMREX_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nnz_global, 1, ParallelDescriptor::Mpi_typemap<Long>::type(),
                  MPI_SUM, comm);
#endif
    auto write_entries = [&] (std::ofstream& ofs, Vector<Long> const& rows,
                              Vector<Long> const& cols, Vector<Real> const& vals)
    {
        for (Long n = 0; n < vals.size(); ++n) {
            ofs << rows[n]+1 << ' ' << cols[n]+1 << ' ' << vals[n] << '\n';
        }
    };

    Vector<Long> rows(csr.nnz()), cols(csr.col_index);
    Vector<Real> vals(csr.values);
    for (Long r = 0; r < csr.nrows; ++r) {
        for (Long n = csr.row_offset[r]; n < csr.row_offset[r+1]; ++n) { rows[n] = csr.row_begin + r; }
    }

    auto header = [&] (std::ofstream& ofs, Long nnz, bool part) {
        ofs << "%%MatrixMarket matrix coordinate real general\n";
        if (part) { ofs << "% part " << myproc << " of " << nprocs << " (rows owned by this rank)\n"; }
        ofs << nglobal << ' ' << nglobal << ' ' << nnz << '\n';
        ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);
    };

    if (nnz_global > gather_max || nprocs == 1) {
        if (nprocs == 1) {
            std::ofstream ofs(file);
            header(ofs, csr.nnz(), false);
            write_entries(ofs, rows, cols, vals);
        } else {
            std::ofstream ofs(file + "." + std::to_string(myproc));
            header(ofs, csr.nnz(), true);
            write_entries(ofs, rows, cols, vals);
        }
        return;
    }

#ifdef AMREX_USE_MPI
    // Gather to rank 0
    int nloc = static_cast<int>(csr.nnz());
    Vector<int> counts(nprocs), displs(nprocs, 0);
    MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    int ntot = 0;
    if (myproc == 0) {
        for (int p = 0; p < nprocs; ++p) { displs[p] = ntot; ntot += counts[p]; }
    }
    Vector<Long> grows(ntot), gcols(ntot);
    Vector<Real> gvals(ntot);
    auto const lt = ParallelDescriptor::Mpi_typemap<Long>::type();
    auto const rt = ParallelDescriptor::Mpi_typemap<Real>::type();
    MPI_Gatherv(rows.data(), nloc, lt, grows.data(), counts.data(), displs.data(), lt, 0, comm);
    MPI_Gatherv(cols.data(), nloc, lt, gcols.data(), counts.data(), displs.data(), lt, 0, comm);
    MPI_Gatherv(vals.data(), nloc, rt, gvals.data(), counts.data(), displs.data(), rt, 0, comm);
    if (myproc == 0) {
        std::ofstream ofs(file);
        header(ofs, ntot, false);
        write_entries(ofs, grows, gcols, gvals);
    }
#endif
}

/**
 * \brief Write a distributed vector (rows row_begin .. row_begin+size-1 on
 * this rank, ranks in row order) as a Matrix Market dense array, gathered
 * on rank 0.
 */
inline void write_vector_market (std::string const& file, Vector<Real> const& v,
                                 Long nglobal, MPI_Comm comm)
{
    int nprocs = 1, myproc = 0;
#ifdef AMREX_USE_MPI
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myproc);
#else
    amrex::ignore_unused(comm);
#endif
    Vector<Real> g;
#ifdef AMREX_USE_MPI
    int nloc = static_cast<int>(v.size());
    Vector<int> counts(nprocs), displs(nprocs, 0);
    MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    if (myproc == 0) {
        int ntot = 0;
        for (int p = 0; p < nprocs; ++p) { displs[p] = ntot; ntot += counts[p]; }
        g.resize(ntot);
    }
    MPI_Gatherv(v.data(), nloc, ParallelDescriptor::Mpi_typemap<Real>::type(), g.data(),
                counts.data(), displs.data(), ParallelDescriptor::Mpi_typemap<Real>::type(), 0, comm);
#else
    g = v;
#endif
    if (myproc == 0) {
        std::ofstream ofs(file);
        ofs << "%%MatrixMarket matrix array real general\n";
        ofs << nglobal << " 1\n";
        ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);
        for (auto const& x : g) { ofs << x << '\n'; }
    }
}

}

#endif
