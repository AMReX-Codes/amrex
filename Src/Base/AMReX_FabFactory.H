
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
 * \brief Classification of FAB cell content, used primarily by the EB (embedded boundary) machinery.
 */
enum class FabType : int {
    covered = -1,       //!< All cells are covered (inside the EB).
    regular = 0,        //!< All cells are regular (no EB intersection).
    singlevalued = 1,   //!< Cells contain a single cut cell.
    multivalued = 2,    //!< Cells contain multiple cut cells.
    undefined = 100     //!< Type has not been determined.
};

/**
 * \struct amrex::FabInfo
 * \ingroup amrex_fab_containers
 * \brief Construction hints passed to FabFactory::create.
 */
struct FabInfo
{
    bool alloc = true;      //!< Whether to allocate storage immediately.
    bool shared = false;    //!< Whether the FAB data is shared (alias).
    Arena* arena = nullptr; //!< Arena to use for allocation (null = default).

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
    //! Allocate and return a new FAB covering \p box with \p ncomps components.
    //! The \p info hints control allocation behavior and \p box_index identifies
    //! the position inside the owning container.
    AMREX_NODISCARD
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
    //! Return an alias FAB that views components [\p scomp, \p scomp+\p ncomp) of \p rhs.
    AMREX_NODISCARD
    virtual FAB* create_alias (FAB const& /*rhs*/, int /*scomp*/, int /*ncomp*/) const { return nullptr; }
    //! Destroy a FAB previously created by this factory; \p fab may be nullptr.
    virtual void destroy (FAB* fab) const = 0;
    //! Return a heap-allocated copy of this factory.
    AMREX_NODISCARD virtual FabFactory<FAB>* clone () const = 0;
    //! Return the number of bytes needed to store \p ncomps components in \p box (optionally using \p box_index for bookkeeping).
    AMREX_NODISCARD virtual Long nBytes (const Box& box, int ncomps, int /*box_index*/) const {
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
    AMREX_NODISCARD
    FAB* create (const Box& box, int ncomps, const FabInfo& info, int /*box_index*/) const override
    {
        return new FAB(box, ncomps, info.alloc, info.shared, info.arena);
    }

    AMREX_NODISCARD
    FAB* create_alias (FAB const& rhs, int scomp, int ncomp) const override
    {
        return new FAB(rhs, amrex::make_alias, scomp, ncomp);
    }

    void destroy (FAB* fab) const override
    {
        delete fab;
    }

    AMREX_NODISCARD
    DefaultFabFactory<FAB>* clone () const override {
        return new DefaultFabFactory<FAB>();
    }
};

}

#endif
