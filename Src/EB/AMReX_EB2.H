#ifndef AMREX_EB2_H_
#define AMREX_EB2_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_Level.H>

#include <cmath>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <string>

namespace amrex::EB2 {

/// \cond DOXYGEN_IGNORE
extern AMREX_EXPORT int max_grid_size;
/// \endcond

//! Initialize EB2’s global state (call before constructing geometries by amrex::Initialize).
void Initialize ();
//! Tear down EB2 global resources (called automatically at shutdown).
void Finalize ();

class IndexSpace
{
public:
    virtual ~IndexSpace() = default;

    IndexSpace () noexcept = default;
    IndexSpace (IndexSpace const&) = delete;
    IndexSpace (IndexSpace &&) = delete;
    IndexSpace& operator= (IndexSpace const&) = delete;
    IndexSpace& operator= (IndexSpace &&) = delete;

    // This function will take the ownership of the IndexSpace
    // pointer, and put it on the top of the stack (i.e., back of the
    // vector).  If the pointer is already in the stack, it will be
    // moved to the top.
    static void push (IndexSpace* ispace);

    static void push (std::unique_ptr<IndexSpace> ispace);

    // This erases `ispace` from the stack.
    static void erase (IndexSpace* ispace);

    static void pop () noexcept { m_instance.pop_back(); }
    static void clear () noexcept { m_instance.clear(); }
    static IndexSpace& top () {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_instance.empty(),
                                         "Have you forgot to call EB2::build? It's required even if the geometry is all regular.");
        return *(m_instance.back());
    }
    static bool empty () noexcept { return m_instance.empty(); }
    static int size () noexcept { return static_cast<int>(m_instance.size()); }

    [[nodiscard]] virtual const Level& getLevel (const Geometry & geom) const = 0;
    [[nodiscard]] virtual const Geometry& getGeometry (const Box& domain) const = 0;
    [[nodiscard]] virtual const Box& coarsestDomain () const = 0;
    virtual void addFineLevels (int num_new_fine_levels) = 0;
    virtual void addRegularCoarseLevels (int num_new_coarse_levels) = 0;

    virtual void setShift (int, int) {
        amrex::Abort("IndexSpace::setShift: not supported");
    }

protected:
    static AMREX_EXPORT Vector<std::unique_ptr<IndexSpace> > m_instance;
};

//! Return the top IndexSpace if one has been built (nullptr otherwise).
const IndexSpace* TopIndexSpaceIfPresent () noexcept;
//! Return the top IndexSpace if one has been built (nullptr otherwise).
inline const IndexSpace* TopIndexSpace () noexcept { return TopIndexSpaceIfPresent(); }

template <typename G>
class IndexSpaceImp
    : public IndexSpace
{
public:

    IndexSpaceImp (const G& gshop, const Geometry& geom,
                   int required_coarsening_level, int max_coarsening_level,
                   int ngrow, bool build_coarse_level_by_coarsening,
                   bool extend_domain_face, int num_coarsen_opt);

    IndexSpaceImp (const G& gshop, const Vector<Geometry>& geom,
                   int ngrow,
                   bool extend_domain_face, int num_coarsen_opt);

    IndexSpaceImp (IndexSpaceImp<G> const&) = delete;
    IndexSpaceImp (IndexSpaceImp<G> &&) = delete;
    void operator= (IndexSpaceImp<G> const&) = delete;
    void operator= (IndexSpaceImp<G> &&) = delete;

    ~IndexSpaceImp () override = default;

    [[nodiscard]] const Level& getLevel (const Geometry& geom) const final;
    [[nodiscard]] const Geometry& getGeometry (const Box& dom) const final;
    [[nodiscard]] const Box& coarsestDomain () const final {
        return m_geom.back().Domain();
    }
    void addFineLevels (int num_new_fine_levels) final;
    void addRegularCoarseLevels (int num_new_coarse_levels) final;

    void setShift (int direction, int ncells) override;

    using F = typename G::FunctionType;

private:

    G m_gshop;
    bool m_build_coarse_level_by_coarsening;
    bool m_extend_domain_face;
    int m_num_coarsen_opt;

    Vector<GShopLevel<G> > m_gslevel;
    Vector<Geometry> m_geom;
    Vector<Box> m_domain;
    Vector<int> m_ngrow;
};

#include <AMReX_EB2_IndexSpaceI.H>

//! Return whether EB2 builders should extend the domain face when outlining geometry.
//! The flag seeds the default `extend_domain_face` argument in EB2::Build.  It returns true
//! unless overridden at runtime via the `eb2.extend_domain_face` parameter, letting users
//! control whether embedded boundaries outside the domain are extrapolated from the domain
//! face (i.e., extended beyond the physical boundary).
bool ExtendDomainFace ();

//! Default number of coarsening options to probe when building the EB hierarchy.
//! The value seeds the default `num_coarsen_opt` argument in EB2::Build and defaults to
//! 0 (no extra coarse probes) unless overridden via the runtime parameter
//! `eb2.num_coarsen_opt`.  It specifies how many successive coarse representations we
//! may inspect to avoid missing fine EB structures inside a sub-box when searching for
//! embedded boundaries.  As soon as EB is detected the hierarchy is rebuilt at full
//! resolution, so varying this value only affects the probability of missing features
//! during the coarse scan—not the final EB geometry if nothing is missed.
int NumCoarsenOpt ();

/**
 * \brief Build EB geometry from a GeometryShop \p gshop and a single AMR Geometry.
 *
 * \tparam G GeometryShop type (see EB2::GeometryShop helpers).
 * \param gshop  Construct describing the implicit surface.
 * \param geom   Finest-level Geometry.
 * \param required_coarsening_level Number of coarse levels required by the application.
 * \param max_coarsening_level      Maximum level the implicit function supports.
 * \param ngrow  Number of EB grow cells to build (default 4).
 * \param build_coarse_level_by_coarsening When true, coarsen from finer levels; otherwise rebuild.
 * \param extend_domain_face Whether to extrapolate the EB description from the domain face to the outside of the domain (defaults to EB2::ExtendDomainFace / runtime param `eb2.extend_domain_face`).
 * \param num_coarsen_opt Number of coarse probes to try (defaults to EB2::NumCoarsenOpt /
 *        runtime `eb2.num_coarsen_opt`) before rebuilding at full resolution.
 */
template <typename G>
void
Build (const G& gshop, const Geometry& geom,
       int required_coarsening_level, int max_coarsening_level,
       int ngrow = 4, bool build_coarse_level_by_coarsening = true,
       bool extend_domain_face = ExtendDomainFace(),
       int num_coarsen_opt = NumCoarsenOpt())
{
    BL_PROFILE("EB2::Initialize()");
    IndexSpace::push(std::make_unique<IndexSpaceImp<G>>
                     (gshop, geom, required_coarsening_level, max_coarsening_level,
                      ngrow, build_coarse_level_by_coarsening,extend_domain_face,
                      num_coarsen_opt));
}

/**
 * \brief Build EB geometry from a GeometryShop and a vector of Geometries (one per level).
 *
 * \tparam G GeometryShop type.
 * \param gshop  Construct describing the implicit surface.
 * \param geom   Vector of Geometries from finest to coarsest (will be sorted by refinement).
 * \param ngrow  Number of EB grow cells to build (default 4).
 * \param extend_domain_face Whether to extrapolate the EB description beyond the domain by extending from the domain face.
 * \param num_coarsen_opt Number of coarse probes to try (defaults to EB2::NumCoarsenOpt /
 *        runtime `eb2.num_coarsen_opt`) before rebuilding at full resolution.
 */
template <typename G>
void
Build (const G& gshop, Vector<Geometry> geom,
       int ngrow = 4,
       bool extend_domain_face = ExtendDomainFace(),
       int num_coarsen_opt = NumCoarsenOpt())
{
    BL_PROFILE("EB2::Initialize()");
    std::sort(geom.begin(), geom.end(), [] (Geometry const& a, Geometry const& b) { return a.Domain().numPts() > b.Domain().numPts(); });
    IndexSpace::push(std::make_unique<IndexSpaceImp<G>>
                     (gshop, geom, ngrow, extend_domain_face, num_coarsen_opt));
}

/**
 * \brief Build EB geometry from the implicit function stored in the Geometry object itself.
 *
 * \param geom   Finest-level Geometry.
 * \param required_coarsening_level Number of coarse levels required by the application.
 * \param max_coarsening_level      Maximum level the implicit function supports.
 * \param ngrow  Number of EB grow cells to build (default 4).
 * \param build_coarse_level_by_coarsening When true, coarsen from finer levels; otherwise rebuild.
 * \param extend_domain_face Whether to extrapolate the EB description beyond the domain by extending from the domain face.
 * \param num_coarsen_opt Number of coarse probes to try (defaults to EB2::NumCoarsenOpt /
 *        runtime `eb2.num_coarsen_opt`) before rebuilding at full resolution.
 * \param support_mvmc Whether to build multi-valued multi-cut support.
 */
void Build (const Geometry& geom,
            int required_coarsening_level,
            int max_coarsening_level,
            int ngrow = 4,
            bool build_coarse_level_by_coarsening = true,
            bool extend_domain_face = ExtendDomainFace(),
            int num_coarsen_opt = NumCoarsenOpt(),
            bool support_mvmc = false);

/**
 * \brief Build EB geometry dedicated to multi-valued multi-cut support.
 *
 * Parameters mirror Build() but always enables multi-valued structures.
 */
void BuildMultiValuedMultiCut (const Geometry& geom,
                               int required_coarsening_level,
                               int max_coarsening_level,
                               int ngrow = 4,
                               bool build_coarse_level_by_coarsening = true,
                               bool extend_domain_face = ExtendDomainFace(),
                               int num_coarsen_opt = NumCoarsenOpt());

/**
 * \brief Rebuild an IndexSpace from a checkpoint file.
 *
 * \param fname   Path to the checkpoint file.
 * \param geom    Geometry describing the finest level.
 * \param required_coarsening_level Number of coarse levels required.
 * \param max_coarsening_level      Maximum supported by the stored implicit function.
 * \param ngrow  Number of grow cells.
 * \param build_coarse_level_by_coarsening Whether to coarsen or rebuild coarse levels.
 * \param extend_domain_face Whether to extrapolate the EB description beyond the domain by extending from the domain face.
 */
void BuildFromChkptFile (std::string const& fname,
                         const Geometry& geom,
                         int required_coarsening_level,
                         int max_coarsening_level,
                         int ngrow = 4,
                         bool build_coarse_level_by_coarsening = true,
                         bool extend_domain_face = ExtendDomainFace());

//! Return the maximum coarsening level supported by the implicit geometry stored in \p geom.
int maxCoarseningLevel (const Geometry& geom);
//! Same as above but for a specific IndexSpace pointer \p ebis.
int maxCoarseningLevel (IndexSpace const* ebis, const Geometry& geom);

//! Add \p num_new_fine_levels to the top of the current IndexSpace hierarchy.
void addFineLevels (int num_new_fine_levels);

//! Add \p num_new_coarse_levels regular levels below the current coarsest level.
void addRegularCoarseLevels (int num_new_coarse_levels);

}

#endif
