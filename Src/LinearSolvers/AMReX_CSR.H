#ifndef AMREX_CSR_H_
#define AMREX_CSR_H_
#include <AMReX_Config.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_INT.H>

#include <type_traits>

namespace amrex {

template <typename T>
struct CsrView {
    using U = std::conditional_t<std::is_const_v<T>, Long const, Long>;
    T* AMREX_RESTRICT mat = nullptr;
    U* AMREX_RESTRICT col_index = nullptr;
    U* AMREX_RESTRICT row_offset = nullptr;
    Long nnz = 0;
    Long nrows = 0;
};

template <typename T, template <typename> class V>
struct CSR {
    V<T> mat;
    V<Long> col_index;
    V<Long> row_offset;
    Long nnz = 0;

    [[nodiscard]] Long nrows () const {
        return row_offset.empty() ? Long(0) : Long(row_offset.size())-1;
    }

    void resize (Long num_rows, Long num_non_zeros) {
        mat.resize(num_non_zeros);
        col_index.resize(num_non_zeros);
        row_offset.resize(num_rows+1);
        nnz = num_non_zeros;
    }

    CsrView<T> view () {
        return CsrView<T>{mat.data(), col_index.data(), row_offset.data(),
                          nnz, Long(row_offset.size())-1};
    }

    [[nodiscard]] CsrView<T const> view () const {
        return CsrView<T>{mat.data(), col_index.data(), row_offset.data(),
                          nnz, Long(row_offset.size())-1};
    }

    [[nodiscard]] CsrView<T const> const_view () const {
        return CsrView<T const>{mat.data(), col_index.data(), row_offset.data(),
                                nnz, Long(row_offset.size())-1};
    }
};

template <typename C, typename T, template<typename> class AD, template<typename> class AS,
          std::enable_if_t<std::is_same_v<C,Gpu::HostToDevice> ||
                           std::is_same_v<C,Gpu::DeviceToHost> ||
                           std::is_same_v<C,Gpu::DeviceToDevice>, int> = 0>
void duplicateCSR (C c, CSR<T,AD>& dst, CSR<T,AS> const& src)
{
    dst.mat.resize(src.mat.size());
    dst.col_index.resize(src.col_index.size());
    dst.row_offset.resize(src.row_offset.size());
    Gpu::copyAsync(c,
                   src.mat.begin(),
                   src.mat.end(),
                   dst.mat.begin());
    Gpu::copyAsync(c,
                   src.col_index.begin(),
                   src.col_index.end(),
                   dst.col_index.begin());
    Gpu::copyAsync(c,
                   src.row_offset.begin(),
                   src.row_offset.end(),
                   dst.row_offset.begin());
    dst.nnz = src.nnz;
}

}

#endif
