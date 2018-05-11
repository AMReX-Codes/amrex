#ifndef AMREX_EB2_H_
#define AMREX_EB2_H_

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_ISLevel.H>

#include <memory>
#include <functional>
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

    static IndexSpace& Get() { AMREX_ASSERT(m_pinstance);  return *m_pinstance; }
    static void Reset () { m_finalizer(); }

protected:
    static IndexSpace* m_pinstance;
    static std::function<void()> m_finalizer;
};

template <typename G>
class IndexSpaceImp
    : public IndexSpace
{
public:

    void define (G&& gshop, Geometry const& geom, Info const& info);

    bool isDefined () const { return m_impfunc != nullptr; }

    virtual ~IndexSpaceImp () {}

    IndexSpaceImp (IndexSpaceImp<G> const&) = delete;
    IndexSpaceImp<G>& operator= (IndexSpaceImp<G> const&) = delete;

    static IndexSpaceImp<G>& GetInstance () {
        if (!m_instance) {
            m_instance.reset(new IndexSpaceImp<G>());
            AMREX_ALWAYS_ASSERT(m_pinstance==nullptr);
            m_pinstance = m_instance.get();
            m_finalizer = IndexSpaceImp<G>::ResetInstance;
        }
        return *m_instance;
    }

    static void ResetInstance () {
        m_instance.reset();
        m_pinstance = nullptr;
    }

    using F = typename std::remove_reference<G>::type::FunctionType;

private:

    static std::unique_ptr<IndexSpaceImp<G> > m_instance;
    IndexSpaceImp () {}

    Vector<ISLevel<G> > m_islevel;
    std::unique_ptr<F> m_impfunc;
};

template <typename G> std::unique_ptr<IndexSpaceImp<G> > IndexSpaceImp<G>::m_instance;

#include <AMReX_EB2_IndexSpaceI.H>

template <typename G>
void
Initialize (G&& gshop, Geometry const& geom, Info const& info)
{
    auto& ispace = IndexSpaceImp<G>::GetInstance();
    if (!ispace.isDefined()) {
        ispace.define(std::forward<G>(gshop), geom, info);
    }
}

inline void Finalize ()
{
    IndexSpace::Reset();
}

}}

#endif
