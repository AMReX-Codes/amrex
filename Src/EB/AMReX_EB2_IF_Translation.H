#ifndef AMREX_EB2_IF_TRANSLATION_H_
#define AMREX_EB2_IF_TRANSLATION_H_

#include <AMReX_Array.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class TranslationIF
{
public:

    TranslationIF (F&& a_f, const RealArray& a_offset)
        : m_f(std::move(a_f)),
          m_offset(a_offset)
        {}
    TranslationIF (F const& a_f, const RealArray& a_offset)
        : m_f(a_f),
          m_offset(a_offset)
        {}

    ~TranslationIF () {}

    TranslationIF (const TranslationIF& rhs) = default;
    TranslationIF (TranslationIF&& rhs) = default;

    TranslationIF& operator= (const TranslationIF& rhs) = delete;
    TranslationIF& operator= (TranslationIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        return m_f({AMREX_D_DECL(p[0]-m_offset[0],
                                 p[1]-m_offset[1],
                                 p[2]-m_offset[2])});
    }

protected:

    F m_f;
    RealArray m_offset;
};

template <class F>
constexpr TranslationIF<typename std::decay<F>::type>
translate (F&&f, const RealArray& offset)
{
    return TranslationIF<typename std::decay<F>::type>(std::forward<F>(f),offset);
}

}}

#endif
