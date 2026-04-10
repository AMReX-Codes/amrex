#ifndef AMREX_IPARSER_H_
#define AMREX_IPARSER_H_

/**
 * \file AMReX_IParser.H
 * \brief Integer-only runtime expression parser and executor helpers.
 */

#include <AMReX_Arena.H>
#include <AMReX_Array.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_IParser_Exe.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>
#include <set>

namespace amrex {

/**
 * \struct amrex::IParserExecutor
 * \brief Callable wrapper around an integer parser expression with \p N variables.
 *
 * Instances are returned by IParser::compile() / IParser::compileHost() and expose
 * host/device operator() overloads that evaluate to 64-bit integers.
 */
template <int N>
struct IParserExecutor
{
    /**
     * \brief Evaluate constant expressions (\p N == 0) and return the integer result.
     */
    template <int M=N, std::enable_if_t<M==0,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    long long operator() () const noexcept
    {
        AMREX_IF_ON_DEVICE((return iparser_exe_eval(m_device_executor, nullptr);))
        AMREX_IF_ON_HOST((return iparser_exe_eval(m_host_executor, nullptr);))
    }

    /**
     * \brief Evaluate the expression with \p N integral arguments.
     * \param var Scalar arguments supplied in registration order.
     * \return Result promoted to \c long long.
     */
    template <typename... Ts>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<sizeof...(Ts) == N && std::conjunction_v<std::is_integral<Ts>...>,
                     long long>
    operator() (Ts... var) const noexcept
    {
        amrex::GpuArray<long long,N> l_var{var...};
        AMREX_IF_ON_DEVICE((return iparser_exe_eval(m_device_executor, l_var.data());))
        AMREX_IF_ON_HOST((return iparser_exe_eval(m_host_executor, l_var.data());))
    }

    /**
     * \brief Evaluate using an explicit array of integer arguments.
     * \param var Argument array in registration order.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    long long operator() (GpuArray<long long, N> const& var) const noexcept
    {
        AMREX_IF_ON_DEVICE((return iparser_exe_eval(m_device_executor, var.data());))
        AMREX_IF_ON_HOST((return iparser_exe_eval(m_host_executor, var.data());))
    }

    //! True when this is a valid executor.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    explicit operator bool () const {
        AMREX_IF_ON_DEVICE((return m_device_executor != nullptr;))
        AMREX_IF_ON_HOST((return m_host_executor != nullptr;))
    }

    char* m_host_executor = nullptr; ///< Pointer to host bytecode.
#ifdef AMREX_USE_GPU
    char* m_device_executor = nullptr; ///< Pointer to device bytecode (if copied).
#endif
};

/// \ingroup amrex_utilities
/**
 * \brief Integer-only variant of amrex::Parser.
 *
 * IParser parses expressions composed of integer literals, integer variables,
 * and integer-returning functions. It emits executors that evaluate to long long
 * on host or GPU.
 */
class IParser
{
public:
    /**
     * \brief Construct a parser by immediately parsing \p func_body.
     * \param func_body Expression to parse (empty to defer define()).
     */
    IParser (std::string const& func_body);
    //! Default-construct; call define() before compile().
    IParser () = default;
    /**
     * \brief Parse and own a new integer expression, replacing any previous state.
     * \param func_body Expression text.
     */
    void define (std::string const& func_body);

    //! True when an expression has been parsed successfully.
    explicit operator bool () const;

    /**
     * \brief Bind a named integer constant.
     * \param name Identifier referenced inside the expression.
     * \param c    Constant value substituted during compilation.
     */
    void setConstant (std::string const& name, long long c);

    /**
     * \brief Register the ordered list of integer variables referenced by the expression.
     *
     * \param vars Variable names supplied in parser argument order.
     */
    void registerVariables (Vector<std::string> const& vars);

    //! Print the parse tree to stdout for debugging.
    void print () const;

    /**
     * \brief Return the maximum parse-tree depth.
     * \return Depth measured in tree levels (0 if undefined).
     */
    [[nodiscard]] int depth () const;
    /**
     * \brief Return the maximum parser stack usage.
     * \return Stack size required during execution (0 if undefined).
     */
    [[nodiscard]] int maxStackSize () const;

    /**
     * \brief Return the sanitized expression string.
     * \return Copy of the parsed expression or an empty string if undefined.
     */
    [[nodiscard]] std::string expr () const;

    /**
     * \brief Return the set of symbols referenced by the expression.
     * \return Sorted symbol names.
     */
    [[nodiscard]] std::set<std::string> symbols () const;

    /**
     * \brief Compile the expression into a host/device IParserExecutor.
     * \tparam N Number of registered variables.
     */
    template <int N> [[nodiscard]] IParserExecutor<N> compile () const;

    /**
     * \brief Compile the expression into a host-only executor.
     * \tparam N Number of registered variables.
     */
    template <int N> [[nodiscard]] IParserExecutor<N> compileHost () const;

private:

    /// \cond DOXYGEN_IGNORE
    struct Data {
        std::string m_expression;
        struct amrex_iparser* m_iparser = nullptr;
        int m_nvars = 0;
        bool m_use_arena = true;
        char* m_host_executor = nullptr;
#ifdef AMREX_USE_GPU
        char* m_device_executor = nullptr;
#endif
        int m_max_stack_size = 0;
        int m_exe_size = 0;
        Data () = default;
        ~Data ();
        Data (Data const&) = delete;
        Data (Data &&) = delete;
        Data& operator= (Data const&) = delete;
        Data& operator= (Data &&) = delete;
    };
    /// \endcond

    std::shared_ptr<Data> m_data;
};

template <int N>
IParserExecutor<N>
IParser::compileHost () const
{
    if (m_data && m_data->m_iparser) {
        AMREX_ALWAYS_ASSERT(N == m_data->m_nvars);

        if (!(m_data->m_host_executor)) {
            int stack_size;
            m_data->m_exe_size = static_cast<int>
                (iparser_exe_size(m_data->m_iparser, m_data->m_max_stack_size,
                                  stack_size));

            if (m_data->m_max_stack_size > AMREX_IPARSER_STACK_SIZE) {
                amrex::Abort("amrex::IParser: AMREX_IPARSER_STACK_SIZE, "
                             + std::to_string(AMREX_IPARSER_STACK_SIZE) + ", is too small for "
                             + m_data->m_expression);
            }
            if (stack_size != 0) {
                amrex::Abort("amrex::IParser: something went wrong with iparser stack! "
                             + std::to_string(stack_size));
            }

            m_data->m_host_executor = (char*)The_Pinned_Arena()->alloc(m_data->m_exe_size);
            if (m_data->m_host_executor == nullptr) { // Arena is not ready yet
                m_data->m_host_executor = (char*) std::malloc(m_data->m_exe_size);
                m_data->m_use_arena = false;
            }

            try {
                iparser_compile(m_data->m_iparser, m_data->m_host_executor);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(std::string(e.what()) + " in IParser expression \""
                                         + m_data->m_expression + "\"");
            }
        }

#ifdef AMREX_USE_GPU
        return IParserExecutor<N>{m_data->m_host_executor, m_data->m_device_executor};
#else
        return IParserExecutor<N>{m_data->m_host_executor};
#endif
    } else {
        return IParserExecutor<N>{};
    }
}

template <int N>
IParserExecutor<N>
IParser::compile () const
{
    auto exe = compileHost<N>();

#ifdef AMREX_USE_GPU
    if (m_data && m_data->m_iparser && !(m_data->m_device_executor)
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

}

#endif
