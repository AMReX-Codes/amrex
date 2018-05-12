#ifndef AMREX_EB2_GEOMETRYSHOP_H_
#define AMREX_EB2_GEOMETRYSHOP_H_

#include <AMReX_EB2_Graph.H>
#include <AMReX_Geometry.H>
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

    int getBoxType (const Box& bx, const Geometry & geom) const;

private:

    F m_f;

};

template <class F>
int
GeometryShop<F>::getBoxType (const Box& bx, const Geometry& geom) const
{
    const Real* problo = geom.ProbLo();
    const Real* dx = geom.CellSize();
    const auto& len3 = bx.length3d();
    const int* blo = bx.loVect();
    int nbody = 0, nzero = 0, nfluid = 0;
    for         (int k = 0; k < len3[2]; ++k) {
        for     (int j = 0; j < len3[1]; ++j) {
            for (int i = 0; i < len3[0]; ++i) {
                RealArray xyz {AMREX_D_DECL(problo[0]+(i+blo[0])*dx[0],
                                            problo[1]+(j+blo[1])*dx[1],
                                            problo[2]+(k+blo[2])*dx[2])};
                Real v = m_f(xyz);
                if (v == 0.0) {
                    ++nzero;
                } else if (v > 0.0) {
                    ++nbody;
                } else {
                    ++nfluid;
                }
            }
        }
    }

    if (nbody == 0) {
        return -1;
    } else if (nfluid == 0) {
        return 1;
    } else {
        return 0;
    }
}

}}

#endif
