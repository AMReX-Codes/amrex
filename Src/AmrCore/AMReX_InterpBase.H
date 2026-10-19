#ifndef AMREX_INTERP_BASE_H_
#define AMREX_INTERP_BASE_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>

#include <typeindex>

/**
 * \file AMReX_InterpBase.H
 * \brief Abstract base for spatial interpolaters used during FillPatch.
 */

namespace amrex {

class InterpBase;

class InterpolaterBoxCoarsener final
    : public BoxConverter
{
public:
    //! Construct a BoxConverter that delegates to \p mapper_ with refinement ratio \p ratio_.
    InterpolaterBoxCoarsener (InterpBase* mapper_, const IntVect& ratio_);
    //! Apply the coarse-box logic to the supplied fine box.
    [[nodiscard]] Box doit (const Box& fine) const override;
    //! Clone this converter (used by FAB factory code).
    [[nodiscard]] BoxConverter* clone () const override;
    //! True if \p other has the same interpolater class and ratio.
    [[nodiscard]] bool equals (const BoxConverter& other) const override;
private:
    InterpBase* mapper;
    IntVect ratio;
    std::type_index m_type;
};

class InterpBase
{
public:
    InterpBase () = default;
    virtual ~InterpBase () = default;
    InterpBase (InterpBase const&) noexcept = default;
    InterpBase (InterpBase &&) noexcept = default;
    InterpBase& operator= (InterpBase const&) noexcept = default;
    InterpBase& operator= (InterpBase &&) noexcept = default;

    /**
    * \brief Return the coarsened box bounding \p fine for scalar refinement ratio \p ratio.
    */
    virtual Box CoarseBox (const Box& fine, int ratio) = 0;

    /**
    * \brief Same as CoarseBox above but accepts per-direction ratios \p ratio.
    *
    * The result must depend only on the interpolater class and \p ratio,
    * because FillPatch caches coarse-box metadata per class and ratio.
    */
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio) = 0;

    //! Convenience helper that returns a BoxConverter bound to this interpolater.
    InterpolaterBoxCoarsener BoxCoarsener (const IntVect& ratio);

    //! Convert BCRec entries to an integer array that interpolaters can consume.
    static Vector<int> GetBCArray (const Vector<BCRec>& bcr);
};

}

#endif
