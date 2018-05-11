#ifndef AMREX_EB2_GEOMETRYSHOP_H_
#define AMREX_EB2_GEOMETRYSHOP_H_

#include <memory>
#include <type_traits>

namespace amrex { namespace EB2 {

template <class F>
class GeometryShop
{
public:

    GeometryShop (F const& f)
        : m_f(f)
        {}

    GeometryShop (F && f)
        : m_f(std::move(f))
        {}

    ~GeometryShop() {}

    GeometryShop (GeometryShop<F> const& rhs)
        : m_f(rhs.m_f)
        {}

    GeometryShop (GeometryShop<F> && rhs)
        : m_f(std::move(rhs.m_f))
        {}

    GeometryShop<F>& operator= (GeometryShop<F> const& rhs) = delete;
    GeometryShop<F>& operator= (GeometryShop<F> && rhs) = delete;

    using FunctionType = F;

    F const& GetImpFunc () const& { return m_f; }
    F&& GetImpFunc () && { return std::move(m_f); }

private:

    F m_f;

};

}}

#endif
