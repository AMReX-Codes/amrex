#ifndef AMREX_PERIODICITY_H_
#define AMREX_PERIODICITY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Periodicity.H
 * \brief Periodicity descriptor for a computational domain.
 */

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <vector>

namespace amrex {

/**
 * \class amrex::Periodicity
 * \ingroup amrex_real_space
 * \brief Periodicity descriptor for a computational domain.
 *
 * Stores the period length in each coordinate direction.  A value of 0
 * in a direction means the domain is not periodic in that direction.
 * It is assumed that the periodic domain starts at index 0.
 */
class Periodicity
{
public:
    //! Construct a non-periodic Periodicity (all periods set to 0).
    Periodicity () noexcept : period(AMREX_D_DECL(0,0,0)) {}

    //! Construct a Periodicity with period lengths given by \p v.
    explicit Periodicity (const IntVect& v) noexcept : period(v) {}

    //! Return true if any direction is periodic.
    [[nodiscard]] bool isAnyPeriodic () const noexcept
        { return AMREX_D_TERM(period[0]>0, || period[1]>0, || period[2]>0); }

    //! Return true if all directions are periodic.
    [[nodiscard]] bool isAllPeriodic () const noexcept
        { return AMREX_D_TERM(period[0]>0, && period[1]>0, && period[2]>0); }

    //! Return true if direction \p dir is periodic.
    [[nodiscard]] bool isPeriodic (int dir) const noexcept
        { return period[dir]>0; }

    //! Return true if both Periodicity objects have the same periods.
    bool operator==(const Periodicity& rhs) const noexcept
        { return period == rhs.period; }

    //! Return a cell-centered domain Box that is "infinitely" long in non-periodic directions.
    [[nodiscard]] Box Domain () const noexcept;

    //! Return the period lengths as an IntVect.
    [[nodiscard]] IntVect const& intVect () const { return period; }

    /**
     * \brief Return the set of shift vectors needed to tile a ghost region.
     *
     * Returns all IntVect shifts (including zero) that map periodic images
     * of the domain into the region grown by \p nghost cells.
     *
     * \param nghost Number of ghost cells in each direction.
     */
    [[nodiscard]] std::vector<IntVect> shiftIntVect (IntVect const& nghost = IntVect(0)) const;

    //! Return a reference to the canonical non-periodic Periodicity singleton.
    static const Periodicity& NonPeriodic () noexcept;

private:
    IntVect period;
};

}

#endif
