#ifndef AMREX_EB2_IF_SCALE_H_
#define AMREX_EB2_IF_SCALE_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class ScaleIF
{
public:

    ScaleIF (F&& a_f, const RealArray& a_scalefactor)
        : m_f(std::move(a_f)),
          m_sfinv{AMREX_D_DECL(1.0/a_scalefactor[0],
                               1.0/a_scalefactor[1],
                               1.0/a_scalefactor[2])}
        {}
    ScaleIF (F const& a_f, const RealArray& a_scalefactor)
        : m_f(a_f),
          m_sfinv{AMREX_D_DECL(1.0/a_scalefactor[0],
                               1.0/a_scalefactor[1],
                               1.0/a_scalefactor[2])}
        {}

    ~ScaleIF () {}

    ScaleIF (const ScaleIF& rhs) = default;
    ScaleIF (ScaleIF&& rhs) = default;
    ScaleIF& operator= (const ScaleIF& rhs) = delete;
    ScaleIF& operator= (ScaleIF&& rhs) = delete;

    template <class U=F, class = typename std::enable_if<IsGPUable<U>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return m_f(AMREX_D_DECL(x*m_sfinv.x,
                                y*m_sfinv.y,
                                z*m_sfinv.z));
    }

    inline Real operator() (const RealArray& p) const noexcept
    {
        return m_f({AMREX_D_DECL(p[0]*m_sfinv.x,
                                 p[1]*m_sfinv.y,
                                 p[2]*m_sfinv.z)});
    }

protected:

    F m_f;
    XDim3 m_sfinv;
};

template <class F>
struct IsGPUable<ScaleIF<F>, typename std::enable_if<IsGPUable<F>::value>::type>
    : std::true_type {};

template <class F>
constexpr ScaleIF<typename std::decay<F>::type>
scale (F&&f, const RealArray& scalefactor)
{
    return ScaleIF<typename std::decay<F>::type>(std::forward<F>(f),scalefactor);
}

}}

#endif
