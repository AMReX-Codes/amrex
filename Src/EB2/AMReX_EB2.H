#ifndef AMREX_EB2_H_
#define AMREX_EB2_H_

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_Level.H>

#include <algorithm>
#include <memory>
#include <type_traits>

namespace amrex { namespace EB2 {

struct Info
{
    int max_coarsening_level = 0;
    int max_grid_size = 64;

    Info& setMaxCoarseningLevel (int n) { max_coarsening_level = n; return *this; }
    Info& setMaxGridSize (int n) { max_grid_size = n; return *this; }
};

class IndexSpace
{
public:
    virtual ~IndexSpace() {}

    static void push (IndexSpace* ispace) { m_instance.emplace_back(ispace); }
    static void pop () { m_instance.pop_back(); }
    static void clear () { m_instance.clear(); }
    static const IndexSpace& top () { return *(m_instance.back()); }
    static bool empty () { return m_instance.empty(); }
    static int size () { return m_instance.size(); }

    virtual const Level& getLevel (const Geometry&) const = 0;

protected:
    static Vector<std::unique_ptr<IndexSpace> > m_instance;
};

template <typename G>
class IndexSpaceImp
    : public IndexSpace
{
public:

    IndexSpaceImp (G&& gshop, const Geometry& geom, const Info& info);

    IndexSpaceImp (IndexSpaceImp<G> const&) = delete;
    IndexSpaceImp (IndexSpaceImp<G> &&) = delete;
    void operator= (IndexSpaceImp<G> const&) = delete;
    void operator= (IndexSpaceImp<G> &&) = delete;

    virtual ~IndexSpaceImp () {}

    virtual const Level& getLevel (const Geometry& geom) const final;

    using F = typename std::remove_reference<G>::type::FunctionType;

private:

    Vector<GShopLevel<G> > m_gslevel;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    std::unique_ptr<F> m_impfunc;
};

#include <AMReX_EB2_IndexSpaceI.H>

template <typename G>
void
Initialize (G&& gshop, const Geometry& geom, const Info& info)
{
    BL_PROFILE("EB2::Initialize()");
    IndexSpace::push(new IndexSpaceImp<G>(std::forward<G>(gshop), geom, info));
}

inline void Finalize ()
{
    IndexSpace::clear();
}

void Initialize (const Geometry& geom, const Info& info);

const Level& getLevel (const Geometry& geom);


}}

#endif
