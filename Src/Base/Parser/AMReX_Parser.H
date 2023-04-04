#ifndef AMREX_PARSER_H_
#define AMREX_PARSER_H_

#include <AMReX_Arena.H>
#include <AMReX_Array.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Parser_Exe.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>
#include <set>

namespace amrex {

template <int N>
struct ParserExecutor
{
    template <int M=N, typename std::enable_if_t<M==0,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() () const noexcept
    {
#if AMREX_DEVICE_COMPILE
        return parser_exe_eval(m_device_executor, nullptr);
#else
        return parser_exe_eval(m_host_executor, nullptr);
#endif
    }

    template <typename... Ts>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<sizeof...(Ts) == N && !amrex::Same<float,Ts...>::value, double>
    operator() (Ts... var) const noexcept
    {
        amrex::GpuArray<double,N> l_var{var...};
#if AMREX_DEVICE_COMPILE
        return parser_exe_eval(m_device_executor, l_var.data());
#else
        return parser_exe_eval(m_host_executor, l_var.data());
#endif
    }

    template <typename... Ts>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<sizeof...(Ts) == N &&  amrex::Same<float,Ts...>::value, float>
    operator() (Ts... var) const noexcept
    {
        amrex::GpuArray<double,N> l_var{var...};
#if AMREX_DEVICE_COMPILE
        return static_cast<float>(parser_exe_eval(m_device_executor, l_var.data()));
#else
        return static_cast<float>(parser_exe_eval(m_host_executor, l_var.data()));
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (GpuArray<double,N> const& var) const noexcept
    {
#if AMREX_DEVICE_COMPILE
        return parser_exe_eval(m_device_executor, var.data());
#else
        return parser_exe_eval(m_host_executor, var.data());
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    explicit operator bool () const {
#if AMREX_DEVICE_COMPILE
        return m_device_executor != nullptr;
#else
        return m_host_executor != nullptr;
#endif
    }

    char* m_host_executor = nullptr;
#ifdef AMREX_USE_GPU
    char* m_device_executor = nullptr;
#endif
};

class Parser
{
public:
    Parser (std::string const& func_body);
    Parser () = default;
    void define (std::string const& func_body);

    explicit operator bool () const;

    void setConstant (std::string const& name, double c);

    void registerVariables (Vector<std::string> const& vars);

    void print () const;

    [[nodiscard]] int depth () const;
    [[nodiscard]] int maxStackSize () const;

    [[nodiscard]] std::string expr () const;

    [[nodiscard]] std::set<std::string> symbols () const;

    //! This compiles for both GPU and CPU
    template <int N> [[nodiscard]] ParserExecutor<N> compile () const;

    //! This compiles for CPU only
    template <int N> [[nodiscard]] ParserExecutor<N> compileHost () const;

private:

    struct Data {
        std::string m_expression;
        struct amrex_parser* m_parser = nullptr;
        int m_nvars = 0;
        mutable char* m_host_executor = nullptr;
#ifdef AMREX_USE_GPU
        mutable char* m_device_executor = nullptr;
#endif
        mutable int m_max_stack_size = 0;
        mutable int m_exe_size = 0;
        Data () = default;
        ~Data ();
        Data (Data const&) = delete;
        Data (Data &&) = delete;
        Data& operator= (Data const&) = delete;
        Data& operator= (Data &&) = delete;
    };

    std::shared_ptr<Data> m_data;
};

template <int N>
ParserExecutor<N>
Parser::compileHost () const
{
    if (m_data && m_data->m_parser) {
        AMREX_ASSERT(N == m_data->m_nvars);

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

            try {
                parser_compile(m_data->m_parser, m_data->m_host_executor);
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
    if (m_data && m_data->m_parser && !(m_data->m_device_executor)) {
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
