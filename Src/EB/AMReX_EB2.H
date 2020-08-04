#ifndef AMREX_EB2_H_
#define AMREX_EB2_H_

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_Level.H>

#include <cmath>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <string>

namespace amrex { namespace EB2 {

extern int max_grid_size;

void useEB2 (bool);

void Initialize ();
void Finalize ();

class IndexSpace
{
public:
    virtual ~IndexSpace() {}

    // This function will take the ownership of the IndexSpace
    // pointer, and put it on the top of the stack (i.e., back of the
    // vector).  If the pointer is already in the stack, it will be
    // moved to the top.
    static void push (IndexSpace* ispace);

    // This erases `ispace` from the stack.
    static void erase (IndexSpace* ispace);

    static void pop () noexcept { m_instance.pop_back(); }
    static void clear () noexcept { m_instance.clear(); }
    static const IndexSpace& top () {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_instance.empty(),
                                         "Have you forgot to call EB2::build? It's required even if the geometry is all regular.");
        return *(m_instance.back());
    }
    static bool empty () noexcept { return m_instance.empty(); }
    static int size () noexcept { return m_instance.size(); }

    virtual const Level& getLevel (const Geometry & geom) const = 0;
    virtual const Geometry& getGeometry (const Box& domain) const = 0;
    virtual const Box& coarsestDomain () const = 0;

protected:
    static Vector<std::unique_ptr<IndexSpace> > m_instance;
};

const IndexSpace* TopIndexSpaceIfPresent () noexcept;
inline const IndexSpace* TopIndexSpace () noexcept { return TopIndexSpaceIfPresent(); }

template <typename G>
class IndexSpaceImp
    : public IndexSpace
{
public:

    IndexSpaceImp (const G& gshop, const Geometry& geom,
                   int required_coarsening_level, int max_coarsening_level,
                   int ngrow, bool build_coarse_level_by_coarsening,
                   bool extend_domain_face);

    IndexSpaceImp (IndexSpaceImp<G> const&) = delete;
    IndexSpaceImp (IndexSpaceImp<G> &&) = delete;
    void operator= (IndexSpaceImp<G> const&) = delete;
    void operator= (IndexSpaceImp<G> &&) = delete;

    virtual ~IndexSpaceImp () {}

    virtual const Level& getLevel (const Geometry& geom) const final;
    virtual const Geometry& getGeometry (const Box& dom) const final;
    virtual const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }

    using F = typename G::FunctionType;

private:

    Vector<GShopLevel<G> > m_gslevel;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
    std::unique_ptr<F> m_impfunc;
};

#include <AMReX_EB2_IndexSpaceI.H>

bool ExtendDomainFace ();

template <typename G>
void
Build (const G& gshop, const Geometry& geom,
       int required_coarsening_level, int max_coarsening_level,
       int ngrow = 4, bool build_coarse_level_by_coarsening = true,
       bool extend_domain_face = ExtendDomainFace())
{
    BL_PROFILE("EB2::Initialize()");
    IndexSpace::push(new IndexSpaceImp<G>(gshop, geom,
                                          required_coarsening_level,
                                          max_coarsening_level,
                                          ngrow, build_coarse_level_by_coarsening,
                                          extend_domain_face));
}

void Build (const Geometry& geom,
            int required_coarsening_level,
            int max_coarsening_level,
            int ngrow = 4,
            bool build_coarse_level_by_coarsening = true);

int maxCoarseningLevel (const Geometry& geom);
int maxCoarseningLevel (IndexSpace const* ebis, const Geometry& geom);

}}

#endif
