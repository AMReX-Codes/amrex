
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_FabFactory.H
 * \brief Abstract factory interface for creating and destroying FAB objects.
 */

#include <AMReX_Box.H>
#include <AMReX_Print.H>
#include <AMReX_MakeType.H>
#include <AMReX_Vector.H>
#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>

namespace amrex
{

class FArrayBox;

/**
 * \ingroup amrex_fab_containers
 * \brief Aggregate classification of a FAB or subregion used by embedded-boundary data.
 */
enum class FabType : int {
    covered = -1,       //!< Every cell in the region is covered.
    regular = 0,        //!< Every cell in the region is regular.
    singlevalued = 1,   //!< No multivalued cells; not entirely regular or covered.
    multivalued = 2,    //!< The region contains at least one multivalued cell.
    undefined = 100     //!< The classification has not been determined.
};

/**
 * \struct amrex::FabInfo
 * \ingroup amrex_fab_containers
 * \brief Construction hints passed to FabFactory::create.
 */
struct FabInfo
{
    bool alloc = true;      //!< Whether the factory should allocate FAB storage.
    bool shared = false;    //!< Whether storage is managed as team-shared memory.
    Arena* arena = nullptr; //!< Arena for storage allocation; null selects the default.

    //! Set the \c alloc flag; returns \c *this for chaining.
    FabInfo& SetAlloc (bool a) noexcept {
        alloc = a;
        return *this;
    }

    //! Set the \c shared flag; returns \c *this for chaining.
    FabInfo& SetShared (bool s) noexcept {
        shared = s;
        return *this;
    }

    //! Set the \c arena pointer; returns \c *this for chaining.
    FabInfo& SetArena (Arena* ar) noexcept {
        arena = ar;
        return *this;
    }
};

/**
 * \class amrex::FabFactory
 * \ingroup amrex_fab_containers
 * \brief Abstract factory interface for creating, aliasing, and destroying FAB objects.
 * \tparam FAB FAB type produced by this factory.
 */
template <class FAB>
class FabFactory
{
public:
    FabFactory () noexcept = default;
    FabFactory (FabFactory const&) noexcept = default;
    FabFactory (FabFactory &&) noexcept = default;
    FabFactory& operator= (FabFactory const&) noexcept = default;
    FabFactory& operator= (FabFactory &&) noexcept = default;
    virtual ~FabFactory () noexcept = default;
    /**
     * \brief Create a FAB covering \p box with \p ncomps components.
     *
     * Storage allocation is controlled by \p info.
     *
     * \param box       Region represented by the FAB.
     * \param ncomps    Number of components per cell.
     * \param info      Allocation hints (arena, shared flag, etc.).
     * \param box_index Global box index supplied by the owning container.
     * \return A newly created FAB owned by the caller.
     */
    [[nodiscard]]
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
    /**
     * \brief Create an alias FAB that views components [\p scomp, \p scomp+\p ncomp) of \p rhs.
     *
     * Default implementation returns nullptr, indicating aliasing is unsupported.
     */
    [[nodiscard]]
    virtual FAB* create_alias (FAB const& /*rhs*/, int /*scomp*/, int /*ncomp*/) const { return nullptr; }
    //! Destroy a FAB previously created by this factory.
    virtual void destroy (FAB* fab) const = 0;
    //! Return a heap-allocated copy of this factory.
    [[nodiscard]] virtual FabFactory<FAB>* clone () const = 0;
    //! Return the storage size for \p ncomps components in \p box, or -1 if unavailable.
    [[nodiscard]] virtual Long nBytes (const Box& box, int ncomps, int /*box_index*/) const {
        if constexpr (IsBaseFab_v<FAB>) {
            return box.numPts() * ncomps * Long(sizeof(typename FAB::value_type));
        } else {
            return -1;
        }
    }
};

/**
 * \class amrex::DefaultFabFactory
 * \ingroup amrex_fab_containers
 * \brief Default FabFactory that calls \c new / \c delete on the FAB type directly.
 * \tparam FAB FAB type to construct.
 */
template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    [[nodiscard]]
    FAB* create (const Box& box, int ncomps, const FabInfo& info, int /*box_index*/) const override
    {
        return new FAB(box, ncomps, info.alloc, info.shared, info.arena);
    }

    [[nodiscard]]
    FAB* create_alias (FAB const& rhs, int scomp, int ncomp) const override
    {
        return new FAB(rhs, amrex::make_alias, scomp, ncomp);
    }

    void destroy (FAB* fab) const override
    {
        delete fab;
    }

    [[nodiscard]]
    DefaultFabFactory<FAB>* clone () const override {
        return new DefaultFabFactory<FAB>();
    }
};

}

#endif
