#ifndef AMREX_PARSER_H_
#define AMREX_PARSER_H_

/**
 * \file AMReX_Parser.H
 * \brief Runtime expression parser front-end and compiled executor helpers.
 */

#include <AMReX_Arena.H>
#include <AMReX_Array.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Parser_Exe.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <map>
#include <memory>
#include <string>
#include <set>
#include <utility>

namespace amrex {

/**
 * \struct amrex::ParserExecutor
 * \brief Callable wrapper around a compiled parser expression with \p N variables.
 *
 * Instances are returned by Parser::compile() / Parser::compileHost() and may be
 * invoked on host or device. Device execution is available only when AMREX_USE_GPU
 * is enabled and a device executor buffer has been created.
 *
 * \tparam N Number of variables referenced in the expression.
 */
template <int N>
struct ParserExecutor
{
    /**
     * \brief Evaluate the expression with no arguments.
     *
     * Enabled only for constant expressions (\p N == 0).
     *
     * \return Expression value as double precision.
     */
    template <int M=N, std::enable_if_t<M==0,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() () const noexcept
    {
        AMREX_IF_ON_DEVICE((return parser_exe_eval(m_device_executor, nullptr);))
        AMREX_IF_ON_HOST((return parser_exe_eval(m_host_executor, nullptr);))
    }

    /**
     * \brief Evaluate the expression with \p N double-compatible arguments.
     *
     * \param var Scalar arguments supplied in parser variable order.
     * \return Expression value promoted to double precision.
     */
    template <typename... Ts>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<sizeof...(Ts) == N && !amrex::Same<float,Ts...>::value, double>
    operator() (Ts... var) const noexcept
    {
        amrex::GpuArray<double,N> l_var{var...};
        AMREX_IF_ON_DEVICE((return parser_exe_eval(m_device_executor, l_var.data());))
        AMREX_IF_ON_HOST((return parser_exe_eval(m_host_executor, l_var.data());))
    }

    /**
     * \brief Evaluate the expression with \p N float arguments.
     *
     * Returns a float to avoid unnecessary promotion when all inputs are float.
     *
     * \param var Scalar arguments supplied in parser variable order.
     */
    template <typename... Ts>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<sizeof...(Ts) == N &&  amrex::Same<float,Ts...>::value, float>
    operator() (Ts... var) const noexcept
    {
        amrex::GpuArray<double,N> l_var{var...};
        AMREX_IF_ON_DEVICE((return static_cast<float>(parser_exe_eval(m_device_executor, l_var.data()));))
        AMREX_IF_ON_HOST((return static_cast<float>(parser_exe_eval(m_host_executor, l_var.data()));))
    }

    /**
     * \brief Evaluate the expression using an explicit array of inputs.
     * \param var Array storing the \p N variables in parser order.
     * \return Expression value promoted to double precision.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (GpuArray<double,N> const& var) const noexcept
    {
        AMREX_IF_ON_DEVICE((return parser_exe_eval(m_device_executor, var.data());))
        AMREX_IF_ON_HOST((return parser_exe_eval(m_host_executor, var.data());))
    }

    /**
     * \brief True when a compiled executor buffer is available.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    explicit operator bool () const {
        AMREX_IF_ON_DEVICE((return m_device_executor != nullptr;))
        AMREX_IF_ON_HOST((return m_host_executor != nullptr;))
    }

    char* m_host_executor = nullptr;  ///< Pointer to host-visible bytecode.
#ifdef AMREX_USE_GPU
    char* m_device_executor = nullptr; ///< Pointer to device-visible bytecode (if copied).
#endif
};

/// \ingroup amrex_utilities
/**
 * \brief Front-end for parsing scalar expressions into GPU/CPU executors.
 *
 * A Parser instance owns the expression string, parser state, registered variables,
 * and any user-defined intrinsic functions. After calling define(), registerVariables(),
 * and (optionally) registerUserFn*, the expression can be compiled into a ParserExecutor
 * that evaluates it efficiently on host or device.
 */
class Parser
{
public:
    /**
     * \brief Construct a parser by immediately parsing \p func_body.
     * \param func_body Expression text; pass an empty string to defer parsing.
     */
    Parser (std::string const& func_body);
    //! Default-construct; call define() before compilation.
    Parser () = default;
    /**
     * \brief Parse and own a new expression, replacing any previous state.
     * \param func_body Expression text to parse.
     */
    void define (std::string const& func_body);

    //! True when an expression has been parsed successfully.
    explicit operator bool () const;

    /**
     * \brief Bind a named constant to \p c.
     * \param name Identifier referenced inside the expression.
     * \param c    Constant value substituted when compiling.
     */
    void setConstant (std::string const& name, double c);

    /**
     * \brief Register the ordered list of variable names used by the expression.
     *
     * The ordering determines how arguments are passed to ParserExecutor::operator().
     *
     * \param vars Variable names supplied in parser argument order.
     */
    void registerVariables (Vector<std::string> const& vars);

    /**
     * \brief Register a user-defined unary function available to expressions.
     * \param name Function identifier used inside expressions.
     * \param fh   Host function pointer (double(double)).
     * \param fd   Device function pointer (double(double)); may be nullptr on CPU-only builds.
     */
    void registerUserFn1 (std::string const& name, ParserUserFn1 fh, ParserUserFn1 fd);
    /**
     * \brief Register a user-defined binary function.
     * \param name Function identifier.
     * \param fh   Host implementation taking two doubles.
     * \param fd   Device implementation taking two doubles; may be nullptr on CPU builds.
     */
    void registerUserFn2 (std::string const& name, ParserUserFn2 fh, ParserUserFn2 fd);
    /**
     * \brief Register a user-defined ternary function.
     * \param name Function identifier.
     * \param fh   Host implementation taking three doubles.
     * \param fd   Device implementation taking three doubles; may be nullptr on CPU builds.
     */
    void registerUserFn3 (std::string const& name, ParserUserFn3 fh, ParserUserFn3 fd);
    /**
     * \brief Register a user-defined function with four arguments.
     * \param name Function identifier.
     * \param fh   Host implementation taking four doubles.
     * \param fd   Device implementation taking four doubles; may be nullptr on CPU builds.
     */
    void registerUserFn4 (std::string const& name, ParserUserFn4 fh, ParserUserFn4 fd);

    //! Print a human-readable representation of the parser tree to stdout.
    void print () const;
    //! Dump the compiled bytecode to stdout for debugging.
    void printExe () const;

    /**
     * \brief Return the maximum tree depth of the parsed expression.
     * \return Depth measured in parse-tree levels (0 if undefined).
     */
    [[nodiscard]] int depth () const;
    /**
     * \brief Return the maximum stack usage recorded during compilation.
     * \return Stack depth in parser instructions (0 if undefined).
     */
    [[nodiscard]] int maxStackSize () const;

    /**
     * \brief Return the sanitized expression string currently managed by the Parser.
     * \return Copy of the parsed expression, or an empty string if undefined.
     */
    [[nodiscard]] std::string expr () const;

    /**
     * \brief Return the set of symbols (variables and functions) referenced by the expression.
     * \return Alphabetically sorted symbol names.
     */
    [[nodiscard]] std::set<std::string> symbols () const;

    /**
     * \brief Access the map of registered user-defined functions.
     * \return Map from function name to argument count.
     */
    [[nodiscard]] std::map<std::string,int> const& userFunctions () const;

    /**
     * \brief Compile the current expression into a host/device executor.
     * \tparam N Number of registered variables.
     * \return ParserExecutor callable on both host and GPU (if available).
     */
    template <int N> [[nodiscard]] ParserExecutor<N> compile () const;

    /**
     * \brief Compile the current expression into a host-only executor.
     * \tparam N Number of registered variables.
     * \return ParserExecutor that executes only on the host.
     */
    template <int N> [[nodiscard]] ParserExecutor<N> compileHost () const;

private:

    template <int argc, typename F>
    void register_user_fn (std::string const& name, F fh, F fd);

    /// \cond DOXYGEN_IGNORE
    struct Data {
        std::string m_expression;
        struct amrex_parser* m_parser = nullptr;
        int m_nvars = 0;
        bool m_use_arena = true;
        char* m_host_executor = nullptr;
#ifdef AMREX_USE_GPU
        char* m_device_executor = nullptr;
#endif
        int m_max_stack_size = 0;
        int m_exe_size = 0;
        Vector<char const*> m_locals;
        Data () = default;
        ~Data ();
        Data (Data const&) = delete;
        Data (Data &&) = delete;
        Data& operator= (Data const&) = delete;
        Data& operator= (Data &&) = delete;
    };
    /// \endcond

    std::shared_ptr<Data> m_data;
    std::map<std::string,int> m_ufs;
    std::map<std::string,std::pair<void*,void*>> m_uf_ptrs;
    Vector<std::string> m_vars;
};

template <int N>
ParserExecutor<N>
Parser::compileHost () const
{
    if (m_data && m_data->m_parser) {
        AMREX_ALWAYS_ASSERT(N == m_data->m_nvars);

        // Make sure all user functions have been registered.
        for (auto const& [ufname, nargs] : m_ufs) {
            if (m_uf_ptrs.find(ufname) == m_uf_ptrs.cend()) {
                amrex::Abort("amrex::Parser: Unknown user function "+ufname
                             + " in " + m_data->m_expression);
            }
            amrex::ignore_unused(nargs);
        }

        if (!(m_data->m_host_executor)) {
            int stack_size;
            m_data->m_exe_size = static_cast<int>
                (parser_exe_size(m_data->m_parser, m_data->m_max_stack_size,
                                 stack_size));

            if (m_data->m_max_stack_size > AMREX_PARSER_STACK_SIZE) {
                amrex::Abort("amrex::Parser: AMREX_PARSER_STACK_SIZE, "
                             + std::to_string(AMREX_PARSER_STACK_SIZE) + ", is too small for "
                             + m_data->m_expression);
            }
            if (stack_size != 0) {
                amrex::Abort("amrex::Parser: something went wrong with parser stack! "
                             + std::to_string(stack_size));
            }

            m_data->m_host_executor = (char*)The_Pinned_Arena()->alloc(m_data->m_exe_size);
            if (m_data->m_host_executor == nullptr) { // Arena is not ready yet
                m_data->m_host_executor = (char*) std::malloc(m_data->m_exe_size);
                m_data->m_use_arena = false;
            }

            try {
                m_data->m_locals = parser_compile(m_data->m_parser, m_uf_ptrs,
                                                  m_data->m_host_executor);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(std::string(e.what()) + " in Parser expression \""
                                         + m_data->m_expression + "\"");
            }
        }

#ifdef AMREX_USE_GPU
        return ParserExecutor<N>{m_data->m_host_executor, m_data->m_device_executor};
#else
        return ParserExecutor<N>{m_data->m_host_executor};
#endif
    } else {
        return ParserExecutor<N>{};
    }
}

template <int N>
ParserExecutor<N>
Parser::compile () const
{
    auto exe = compileHost<N>();

#ifdef AMREX_USE_GPU
    if (m_data && m_data->m_parser && !(m_data->m_device_executor)
        && m_data->m_use_arena)
    {
        m_data->m_device_executor = (char*)The_Arena()->alloc(m_data->m_exe_size);
        Gpu::htod_memcpy_async(m_data->m_device_executor, m_data->m_host_executor,
                               m_data->m_exe_size);
        Gpu::streamSynchronize();
        exe.m_device_executor = m_data->m_device_executor;
    }
#endif

    return exe;
}

template <int argc, typename F>
void
Parser::register_user_fn (std::string const& name, F fh, F fd)
{
#if !defined(AMREX_USE_GPU)
    // For pure CPU builds, we allow the user to pass the function pointer with either fh or fd
    if (fh == nullptr) { std::swap(fh, fd); }
#endif

    auto it = m_ufs.find(name);
    if (it != m_ufs.end()) {
        if (it->second == argc) {
            m_uf_ptrs[name] = std::make_pair((void*)fh, (void*)fd);
        } else {
            amrex::Abort("Parser::registerUserFn: incorrect no. of arguments in "+name);
        }
    } else {
        amrex::Abort("Parser::registerUserFn: unknown user define function "+name);
    }
}

}

#endif
