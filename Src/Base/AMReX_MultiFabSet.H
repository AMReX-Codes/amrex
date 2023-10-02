#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_FabArraySet.H>

namespace amrex {

/**
 * \brief MultiFabSet is a grouping of MultiFabs usually made when each component requires a separate IndexType.
 * Often made for input/output into multi-component linear operators. It can also be used to represent
 * edge-centered (i.e., 1-form or vector) and face-centered (i.e., 2-form or psuedovector) fields. While the
 * MultiFab class allows for multiple components, each component is required to live at exactly the same location
 * (IndexType) on the grid due to the current nature of Array4.
 * 
 * Specifically, each component in a MultiFab shares the same BoxArray and thus IndexType. Unfortunately, 
 * differential forms such as 1-forms (vector fields) and 2-forms (psuedovector fields) are ideally used by finite 
 * difference methods when each component has a distinct IndexType. In particular, 1-form components live on the 
 * lines between nodes (i.e., the edges of the grid cells), and 2-form components live on the faces of grid cells. 
 * This requires each component to have slightly different BoxArray's and associated Array4's. 
 * 
 * In summary, MultiFabSet is intended to be an extension of how edge- and face-centered fields have been 
 * represented thus far, which is a collection of single component MultiFab's that live on slightly different FAB's 
 * with distinct IndexType's (nodal flags).
 */
// template <unsigned long encodedIxTypes>
// class MultiFabSet<encodedIxTypes> : public FabArraySet<MultiFab,encodedIxTypes>
class MultiFabSet : public FabArraySet
{
    //! debug parameter
    static constexpr int N = 10;
public:
    using MF = typename amrex::MultiFab;
    using RT = typename FabArraySet::value_type;
    using FABFactory = typename amrex::FabFactory<fab_type>;
    using DM = typename amrex::DistributionMapping;

    /**
    * \brief Constructs an empty MultiFabSet.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArraySet.
    */
    MultiFabSet () noexcept;

    /**
    * \brief Constructs an empty MultiFabSet.  Data can be defined at a later
    * time using the define member functions inherited from FabArraySet.  If
    * `define` is called later with a nulltpr as MFInfo's arena, the default
    * Arena `a` will be used.  If the arena in MFInfo is not a nullptr, the
    * MFInfo's arena will be used.
    */
    explicit MultiFabSet (Arena* a) noexcept;

    /**
     * \brief Construct a new MultiFabSet object
     * 
     * \param bxs BoxArray is treated as all nodal
     * \param dm 
     * \param ncomp Sets the number of components
     * \param ngrow 
     * \param info 
     * \param factory 
     */
    MultiFabSet (const BoxArray& bxs, const DM& dm, int ncomp, const IntVect& ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo& info, const FABFactory& factory);
#else
                 const MFInfo& info = MFInfo(), const FABFactory& factory = FArrayBoxFactory());
#endif

    //! Generate MultiFabSet from Array<MF,N>
    // template <std::size_t N>
    MultiFabSet (Array<MF*,N>& mfarray);
    //! Generate MultiFabSet from const Array<MF,N>
    // template <std::size_t N>
    MultiFabSet (const Array<MF*,N>& mfarray);

     /**
     * \brief Make an alias MultiFabSet.
     * \param mfs Starting MultiFabSet
     * \param maketype must be amrex::make_alias
     * \param scomp The starting component of the alias
     * \param ncomp The number of components in the new aliasing MultiFabs.
     */
    MultiFabSet (const MultiFabSet& mfs, MakeType maketype, int scomp, int ncomp);

    ~MultiFabSet ();

    MultiFabSet (MultiFabSet&& rhs) noexcept;
    MultiFabSet& operator= (MultiFabSet&& rhs) noexcept = default;

    MultiFabSet (const MultiFabSet& rhs) = delete;
    MultiFabSet& operator= (const MultiFabSet& rhs) = delete;

    /**
    * \brief Negates the value of each cell in the valid region of
    * the MultiFab.  The value of nghost specifies the number of
    * cells in the boundary region that should be modified.
    */
    void negate (int nghost = 0);

};

MultiFabSet::MultiFabSet (const MultiFabSet& src, MakeType maketype, int scomp, int ncomp)
{
    for (int i = 0; i < nSet(); ++i) {
        setElem(i, new MF(src.getElem(i), maketype, scomp, ncomp));
    }
    m_NeedToDelete = true;
}

MultiFabSet::MultiFabSet (const BoxArray& bxs, const DM& dm, int ncomp, const IntVect& ngrow, 
                             const MFInfo& info, const FABFactory& factory)
{
    for (int i = 0; i < nSet(); ++i) {
        setElem(i, new MF(amrex::convert(bxs, m_ixtype_set[i]), dm, ncomp, ngrow, info, factory));
    }
    m_NeedToDelete = true;
}

// template <std::size_t N>
MultiFabSet::MultiFabSet (Array<MF*,N>& mfarray) 
{
    for (int i = 0; i < nSet(); ++i) {
        setElem(i, mfarray[i]);
    }
}

// template <std::size_t N>
MultiFabSet::MultiFabSet (const Array<MF*,N>& mfarray)
{
    for (int i = 0; i < nSet(); ++i) {
        setElem(i, mfarray[i]);
    }
}

}