#ifndef AMREX_EB2_IF_SCALE_H_
#define AMREX_EB2_IF_SCALE_H_

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

    Real operator() (const RealArray& p) const
    {
        return m_f({AMREX_D_DECL(p[0]*m_sfinv[0],
                                 p[1]*m_sfinv[1],
                                 p[2]*m_sfinv[2])});
    }

protected:

    F m_f;
    RealArray m_sfinv;
};

template <class F>
constexpr ScaleIF<typename std::decay<F>::type>
scale (F&&f, const RealArray& scalefactor)
{
    return ScaleIF<typename std::decay<F>::type>(std::forward<F>(f),scalefactor);
}

}}

#endif
