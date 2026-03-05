
#include <AMReX_SpMatUtil.H>
#include <AMReX.H>
#include <AMReX_Random.H>

using namespace amrex;

int main (int argc, char *argv[])
{
    amrex::Initialize(argc, argv);

    std::array<int,2> shapes[]
        {{100, 100}, {75, 64}, {64, 75}, {16, 120}, {120, 16}, {6,8}, {8,6}, {7,7}};

    constexpr int nnz = 5;

    auto value = [] AMREX_GPU_DEVICE (Long row, Long col) {
        return double(100000000 + row*10000 + col);
    };

    for (auto const& nn : shapes)
    {
        // We are not doing auto [nrows, ncols] : shapes, because capturing
        // structured bindings in lambda is C++20.
        auto nrows = nn[0];
        auto ncols = nn[1];

        amrex::Print() << "  Shape: " << nrows << " x " << ncols << "\n";

        AlgPartition row_part(nrows);
        AlgPartition col_part(ncols);

        SpMatrix<double> mat(row_part, nnz);

        auto row_begin = mat.globalRowBegin();

        auto *matVals = mat.data();
        auto *matCols = mat.columnIndex();

        ParallelForRNG(mat.numLocalRows(),
                       [=] AMREX_GPU_DEVICE (Long lrow, RandomEngine const& eng)
        {
            auto row = lrow + row_begin; // global row index
            auto col = (row-1 + ncols) % ncols;
            matVals[lrow*nnz  ] = value(row,col);
            matCols[lrow*nnz  ] = col;

            col = row % ncols;
            matVals[lrow*nnz+1] = value(row,col);
            matCols[lrow*nnz+1] = col;

            col = (row+1) % ncols;
            matVals[lrow*nnz+2] = value(row,col);
            matCols[lrow*nnz+2] = col;

            while (col == matCols[lrow*nnz  ] ||
                   col == matCols[lrow*nnz+1] ||
                   col == matCols[lrow*nnz+2]) {
                col = int(amrex::Random_int(ncols,eng));
            }
            matVals[lrow*nnz+3] = value(row,col);
            matCols[lrow*nnz+3] = col;


            while (col == matCols[lrow*nnz  ] ||
                   col == matCols[lrow*nnz+1] ||
                   col == matCols[lrow*nnz+2] ||
                   col == matCols[lrow*nnz+3]) {
                col = int(amrex::Random_int(ncols,eng));
            }
            matVals[lrow*nnz+4] = value(row,col);
            matCols[lrow*nnz+4] = col;
        });

        mat.sortCSR();

#ifdef AMREX_DEBUG
        ParallelDescriptor::Barrier();
        mat.printToFile("mat-"+std::to_string(nrows)+"x"+std::to_string(ncols));
#endif

        auto matt = amrex::transpose(mat, col_part);
#ifdef AMREX_DEBUG
        ParallelDescriptor::Barrier();
        matt.printToFile("mat-"+std::to_string(nrows)+"x"+std::to_string(ncols)+"-t");
#endif

        {
            auto const& csr = matt.const_parcsr();
            auto err = Reduce::Sum<int>(matt.numLocalRows(),
                                        [=] AMREX_GPU_DEVICE (Long i)
            {
                int r = 0;
                for (auto idx = csr.csr0.row_offset[i]; idx < csr.csr0.row_offset[i+1]; ++idx) {
                    auto v1 = csr.csr0.mat[idx];
                    auto v2 = value(csr.csr0.col_index[idx]+csr.col_begin, i+csr.row_begin);
                    if (! amrex::almostEqual(v1,v2)) {
                        ++r;
                    }
                }
                if (csr.row_map && csr.row_map[i] >= 0) {
                    auto ii = csr.row_map[i];
                    for (auto idx = csr.csr1.row_offset[ii]; idx < csr.csr1.row_offset[ii+1]; ++idx) {
                        auto v1 = csr.csr1.mat[idx];
                        auto v2 = value(csr.col_map[csr.csr1.col_index[idx]], i+csr.row_begin);
                        if (! amrex::almostEqual(v1,v2)) {
                            ++r;
                        }
                    }
                }
                return r;
            });

            AMREX_ALWAYS_ASSERT(err == 0);

            Long total_nnz = csr.csr0.nnz + csr.csr1.nnz;
            ParallelDescriptor::ReduceLongSum(total_nnz);
            AMREX_ALWAYS_ASSERT(total_nnz == nrows*nnz);
        }

        auto mattt = amrex::transpose(matt, row_part);
#ifdef AMREX_DEBUG
        ParallelDescriptor::Barrier();
        mattt.printToFile("mat-"+std::to_string(nrows)+"x"+std::to_string(ncols)+"-tt");
#endif

        {
            auto const& csr = mattt.const_parcsr();
            auto err = Reduce::Sum<int>(mattt.numLocalRows(),
                                        [=] AMREX_GPU_DEVICE (Long i)
            {
                int r = 0;
                for (auto idx = csr.csr0.row_offset[i]; idx < csr.csr0.row_offset[i+1]; ++idx) {
                    auto v1 = csr.csr0.mat[idx];
                    auto v2 = value(i+csr.row_begin, csr.csr0.col_index[idx]+csr.col_begin);
                    if (! amrex::almostEqual(v1,v2)) {
                        ++r;
                    }
                }
                if (csr.row_map && csr.row_map[i] >= 0) {
                    auto ii = csr.row_map[i];
                    for (auto idx = csr.csr1.row_offset[ii]; idx < csr.csr1.row_offset[ii+1]; ++idx) {
                        auto v1 = csr.csr1.mat[idx];
                        auto v2 = value(i+csr.row_begin, csr.col_map[csr.csr1.col_index[idx]]);
                        if (! amrex::almostEqual(v1,v2)) {
                            ++r;
                        }
                    }
                }
                return r;
            });

            AMREX_ALWAYS_ASSERT(err == 0);

            Long total_nnz = csr.csr0.nnz + csr.csr1.nnz;
            ParallelDescriptor::ReduceLongSum(total_nnz);
            AMREX_ALWAYS_ASSERT(total_nnz == nrows*nnz);
        }
    }

    amrex::Finalize();
}
