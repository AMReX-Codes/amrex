#ifndef AMREX_SMALL_MATRIX_H_
#define AMREX_SMALL_MATRIX_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ConstexprFor.H>

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <tuple>
#include <type_traits>

namespace amrex {

    enum struct Order { C, F, RowMajor=C, ColumnMajor=F };

    /**
     * \brief Matrix class with compile-time size
     *
     * Note that column vectors and row vectors are special cases of a
     * Matrix.
     *
     * \tparam T Matrix element data type.
     * \tparam NRows Number of rows.
     * \tparam NCols Number of columns.
     * \tparam ORDER Memory layout order. Order::F (i.e., column-major) by default.
     * \tparam StartIndex Starting index. Either 0 or 1.
     */
    template <class T, int NRows, int NCols, Order ORDER = Order::F, int StartIndex = 0>
    struct SmallMatrix
    {
        using value_type = T;
        using reference_type = T&;
        static constexpr int row_size = NRows;
        static constexpr int column_size = NCols;
        static constexpr Order ordering = ORDER;
        static constexpr int starting_index = StartIndex;

        /**
         * \brief Default constructor
         *
         * The data are uninitialized by default. If you want to initialize
         * to zero, you can do `SmallMatrix<T,NRows,NCols> M{};`.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr SmallMatrix () = default;

        /**
         * \brief Constructs column- or row-vector
         *
         * The data are initialized with the given variadic arguments. If
         * the number of argument is less than the size of the vector, the
         * rest of the vector is initialized to zero.
         */
        template <typename... Ts, int MM=NRows, int NN=NCols,
                  std::enable_if_t<MM==1 || NN==1, int> = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr explicit SmallMatrix (Ts... vs)
            : m_mat{vs...}
        {
            static_assert(sizeof...(vs) <= std::max(NRows,NCols));
        }

        /**
         * \brief Constructs SmallMatrix with nested std::initializer_list
         *
         * The initializer list is assumed to be in row-major order, even when
         * the ordering for the SmallMatrix object is colum-major. Below is
         * an example of constructing a matrix with 2 rows and 3 columns.
         \verbatim
             SmallMatrix<double,2,3> M{{11., 12., 13.},
                                       {21., 22., 23.}};
         \endverbatim
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        explicit SmallMatrix (std::initializer_list<std::initializer_list<T>> const& init)
        {
            AMREX_ASSERT(NRows == init.size());
            int i = StartIndex;
            for (auto const& row : init) {
                AMREX_ASSERT(NCols == row.size());
                int j = StartIndex;
                for (auto const& x : row) {
                    (*this)(i,j) = x;
                    ++j;
                }
                ++i;
            }
        }

        //! Returns a const reference to the element at row i and column j.
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
                --j;
            }
            AMREX_ASSERT(i < NRows && j < NCols);
            if constexpr (ORDER == Order::F) {
                return m_mat[i+j*NRows];
            } else {
                return m_mat[j+i*NCols];
            }
        }

        //! Returns a reference to the element at row i and column j.
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
                --j;
            }
            AMREX_ASSERT(i < NRows && j < NCols);
            if constexpr (ORDER == Order::F) {
                return m_mat[i+j*NRows];
            } else {
                return m_mat[j+i*NCols];
            }
        }

        //! Returns a const reference to element i of a vector
        template <int MM=NRows, int NN=NCols, std::enable_if_t<(MM==1 || NN==1), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i) const noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
            }
            AMREX_ASSERT(i < NRows*NCols);
            return m_mat[i];
        }

        //! Returns a reference to element i of a vector
        template <int MM=NRows, int NN=NCols, std::enable_if_t<(MM==1 || NN==1), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i) noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
            }
            AMREX_ASSERT(i < NRows*NCols);
            return m_mat[i];
        }

        //! Returns a const reference to element i of a vector
        template <int MM=NRows, int NN=NCols, std::enable_if_t<(MM==1 || NN==1), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator[] (int i) const noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
            }
            AMREX_ASSERT(i < NRows*NCols);
            return m_mat[i];
        }

        //! Returns a reference to element i of a vector
        template <int MM=NRows, int NN=NCols, std::enable_if_t<(MM==1 || NN==1), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator[] (int i) noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            if constexpr (StartIndex == 1) {
                --i;
            }
            AMREX_ASSERT(i < NRows*NCols);
            return m_mat[i];
        }

        /**
         * Returns a \c const pointer address to the first element of the
         * SmallMatrix object, as if the object is treated as one-dimensional.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin () const noexcept { return m_mat; }

        /**
         * Returns a \c const pointer address right after the last element of the
         * SmallMatrix object, as if the object is treated as one-dimensional.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end () const noexcept { return m_mat + NRows*NCols; }

        /**
         * Returns a pointer address to the first element of the
         * SmallMatrix object, as if the object is treated as one-dimensional.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin () noexcept { return m_mat; }

        /**
         * Returns a pointer address right after the last element of the
         * SmallMatrix object, as if the object is treated as one-dimensional.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end () noexcept { return m_mat + NRows*NCols; }

        //! Set all elements in the matrix to the given value
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>&
        setVal (T val)
        {
            for (auto& x : m_mat) { x = val; }
            return *this;
        }

        //! Returns an identity matrix
        template <int MM=NRows, int NN=NCols, std::enable_if_t<MM==NN, int> = 0>
        static constexpr
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        Identity () noexcept {
            static_assert(StartIndex == 0 || StartIndex == 1);
            SmallMatrix<T,NRows,NCols,ORDER,StartIndex> I{};
            constexpr_for<StartIndex,NRows+StartIndex>(
                [&] (int i) { I(i,i) = T(1); });
            return I;
        }

        //! Returns a matrix initialized with zeros
        static constexpr
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        Zero () noexcept {
            SmallMatrix<T,NRows,NCols,ORDER,StartIndex> Z{};
            return Z;
        }

        //! Returns transposed matrix
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NCols,NRows,ORDER,StartIndex>
        transpose () const
        {
            SmallMatrix<T,NCols,NRows,ORDER,StartIndex> r;
            for (int j = StartIndex; j < NRows+StartIndex; ++j) {
                for (int i = StartIndex; i < NCols+StartIndex; ++i) {
                    r(i,j) = (*this)(j,i);
                }
            }
            return r;
        }

        //! Transposes a square matrix in-place.
        template <int MM=NRows, int NN=NCols, std::enable_if_t<MM==NN,int> = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>&
        transposeInPlace ()
        {
            static_assert(StartIndex == 0 || StartIndex == 1);
            for (int j = 1+StartIndex; j < NCols+StartIndex; ++j) {
                for (int i = StartIndex; i < j; ++i) {
                    amrex::Swap((*this)(i,j), (*this)(j,i));
                }
            }
            return *this;
        }

        //! Returns the product of all elements in the matrix
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T product () const
        {
            T p = 1;
            for (auto const& x : m_mat) {
                p *= x;
            }
            return p;
        }

        //! Returns the sum of all elements in the matrix
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T sum () const
        {
            T s = 0;
            for (auto const& x : m_mat) {
                s += x;
            }
            return s;
        }

        //! Returns the trace of a square matrix
        template <int MM=NRows, int NN=NCols, std::enable_if_t<MM==NN,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T trace () const
        {
            T t = 0;
            constexpr_for<StartIndex,MM+StartIndex>([&] (int i) { t += (*this)(i,i); });
            return t;
        }

        //! Operator += performing matrix addition as in (*this) += rhs
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>&
        operator += (SmallMatrix<T,NRows,NCols,ORDER,StartIndex> const& rhs)
        {
            for (int n = 0; n < NRows*NCols; ++n) {
                m_mat[n] += rhs.m_mat[n];
            }
            return *this;
        }

        //! Binary operator + returning the result of maxtrix addition, lhs+rhs
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        friend SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        operator+ (SmallMatrix<T,NRows,NCols,ORDER,StartIndex>        lhs,
                   SmallMatrix<T,NRows,NCols,ORDER,StartIndex> const& rhs)
        {
            lhs += rhs;
            return lhs;
        }

        //! Operator -= performing matrix subtraction as in (*this) -= rhs
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>&
        operator -= (SmallMatrix<T,NRows,NCols,ORDER,StartIndex> const& rhs)
        {
            for (int n = 0; n < NRows*NCols; ++n) {
                m_mat[n] -= rhs.m_mat[n];
            }
            return *this;
        }

        //! Binary operator - returning the result of maxtrix subtraction, lhs-rhs
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        friend SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        operator- (SmallMatrix<T,NRows,NCols,ORDER,StartIndex>        lhs,
                   SmallMatrix<T,NRows,NCols,ORDER,StartIndex> const& rhs)
        {
            lhs -= rhs;
            return lhs;
        }

        //! Unary minus operator
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        operator- () const
        {
            return (*this) * T(-1);
        }

        //! Operator *= that scales this matrix in place by a scalar.
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SmallMatrix<T,NRows,NCols,ORDER,StartIndex>&
        operator *= (T a)
        {
            for (auto& x : m_mat) {
                x *= a;
            }
            return *this;
        }

        //! Returns the product of a matrix and a scalar
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        friend SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        operator* (SmallMatrix<T,NRows,NCols,ORDER,StartIndex> m, T a)
        {
            m *= a;
            return m;
        }

        //! Returns the product of a scalar and a matrix
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        friend SmallMatrix<T,NRows,NCols,ORDER,StartIndex>
        operator* (T a, SmallMatrix<T,NRows,NCols,ORDER,StartIndex> m)
        {
            m *= a;
            return m;
        }

        //! Returns matrix product of two matrices
        template <class U, int N1, int N2, int N3, Order Ord, int SI>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        friend SmallMatrix<U,N1,N3,Ord,SI>
        operator* (SmallMatrix<U,N1,N2,Ord,SI> const& lhs,
                   SmallMatrix<U,N2,N3,Ord,SI> const& rhs);

        //! Returns the dot product of two vectors
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T dot (SmallMatrix<T,NRows,NCols,ORDER,StartIndex> const& rhs) const
        {
            T r = 0;
            for (int n = 0; n < NRows*NCols; ++n) {
                r += m_mat[n] * rhs.m_mat[n];
            }
            return r;
        }

        template <int N, std::enable_if_t<(N<NRows*NCols),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T const& get () const { return m_mat[N]; }

        template <int N, std::enable_if_t<(N<NRows*NCols),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T& get () { return m_mat[N]; }

    private:
        T m_mat[NRows*NCols];
    };

    template <class U, int N1, int N2, int N3, Order Ord, int SI>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SmallMatrix<U,N1,N3,Ord,SI>
    operator* (SmallMatrix<U,N1,N2,Ord,SI> const& lhs,
               SmallMatrix<U,N2,N3,Ord,SI> const& rhs)
    {
        static_assert(SI == 0 || SI == 1);
        SmallMatrix<U,N1,N3,Ord,SI> r;
        if constexpr (Ord == Order::F) {
            for (int j = SI; j < N3+SI; ++j) {
                constexpr_for<SI,N1+SI>([&] (int i) { r(i,j) = U(0); });
                for (int k = SI; k < N2+SI; ++k) {
                    auto b = rhs(k,j);
                    constexpr_for<SI,N1+SI>([&] (int i)
                    {
                        r(i,j) += lhs(i,k) * b;
                    });
                }
            }
        } else {
            for (int i = SI; i < N1+SI; ++i) {
                constexpr_for<SI,N3+SI>([&] (int j) { r(i,j) = U(0); });
                for (int k = SI; k < N2+SI; ++k) {
                    auto a = lhs(i,k);
                    constexpr_for<SI,N3+SI>([&] (int j)
                    {
                        r(i,j) += a * rhs(k,j);
                    });
                }
            }
        }
        return r;
    }

    template <class T, int NRows, int NCols, Order ORDER, int SI>
    std::ostream& operator<< (std::ostream& os,
                              SmallMatrix<T,NRows,NCols,ORDER,SI> const& mat)
    {
        for (int i = SI; i < NRows+SI; ++i) {
            os << mat(i,SI);
            for (int j = 1+SI; j < NCols+SI; ++j) {
                os << " " << mat(i,j);
            }
            os << "\n";
        }
        return os;
    }

    template <class T, int N, int StartIndex = 0>
    using SmallVector = SmallMatrix<T,N,1,Order::F,StartIndex>;

    template <class T, int N, int StartIndex = 0>
    using SmallRowVector = SmallMatrix<T,1,N,Order::F,StartIndex>;
}

template <class T, int NRows, int NCols, amrex::Order ORDER, int StartIndex>
struct std::tuple_size<amrex::SmallMatrix<T,NRows,NCols,ORDER,StartIndex> >
    : std::integral_constant<std::size_t,NRows*NCols> {};

template <std::size_t N, class T, int NRows, int NCols, amrex::Order ORDER, int StartIndex>
struct std::tuple_element<N, amrex::SmallMatrix<T,NRows,NCols,ORDER,StartIndex> >
{
    using type = T;
};

#endif

/*
 * Notes on why SmallMatrix matrix{} is zero initialized.
 *
 * SmallMatrix is not an aggregate, because it has a user declared default
 * constructor. The rule is that, for `SmallMatrix matrix{}` with an empty
 * brace-enclosed initializer list, value-initialization is performed. The
 * effects of value-initialization of SmallMatrix (which has a user-declared
 * but not user-provided default constructor) are that the matrix object is
 * first zero-initialized and then the object's default constructor is
 * applied. Since the default constructor does nothing, the final result is
 * the object is zero-initialized.
 *
 * Why is SmallMatrix's default constructor user-declared not user-provided?
 * It's because we first declare it with `SmallMatrix () = default`.
 *
 * Reference:
 *   https://en.cppreference.com/w/cpp/language/list_initialization
 *   https://en.cppreference.com/w/cpp/language/value_initialization
 *   https://en.cppreference.com/w/cpp/language/zero_initialization
 */
