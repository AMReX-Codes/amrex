#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class MultiFabSet : public MultiFab
{
public:

    /**
    * \brief Constructs an empty MultiFabSet.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    MultiFabSet () noexcept;

    /**
    * \brief Constructs an empty MultiFabSet.  Data can be defined at a later
    * time using the define member functions inherited from FabArray.  If
    * `define` is called later with a nulltpr as MFInfo's arena, the default
    * Arena `a` will be used.  If the arena in MFInfo is not a nullptr, the
    * MFInfo's arena will be used.
    */
    explicit MultiFabSet (Arena* a) noexcept;

    /**
    * \brief
    * Constructs a MultiFab
    * \param bxs a valid region
    * \param dm a DistribuionMapping
    * \param ncomp number of components
    * \param ngrow number of cells the region grows
    * \param info MFInfo

    * The size of the FArrayBox is given by the Box grown by ngrow, and
    * the number of components is given by ncomp. If info is set to
    * not allocating memory, then no FArrayBoxes are allocated at
    * this time but can be defined later.
    */
    MultiFabSet (const BoxArray&              bxs,
                 const DistributionMapping&   dm,
                 const Vector<IndexType>      ixTypes,
                 const IntVect&               ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&                info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&                info = MFInfo(),
                 const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    /**
     * \brief Make an alias MultiFab. maketype must be
     * amrex::make_alias.  scomp is the starting component of the
     * alias and ncomp is the number of components in the new aliasing
     * MultiFab.
     */
    MultiFabSet (const MultiFabSet& rhs, MakeType maketype, int scomp, int ncomp);

    ~MultiFabSet ();

    MultiFabSet (MultiFabSet&& rhs) noexcept;
    MultiFabSet& operator= (MultiFabSet&& rhs) noexcept = default;

    MultiFabSet (const MultiFabSet& rhs) = delete;
    MultiFabSet& operator= (const MultiFabSet& rhs) = delete;

    void define (const BoxArray&              bxs,
                 const DistributionMapping&   dm,
                 const Vector<IndexType>      ixTypes,
                 const IntVect&               ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&                info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&                info = MFInfo(),
                 const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    [[nodiscard]] const BoxArray& boxArray () const noexcept
    {
        BL_PROFILE("MultiFabSet::boxArray()");
        return boxarray.convert(IndexType::TheNodeType());
    }

    IndexType ixType () = delete;

    void LocalCopy (const MultiFabSet& src, int scomp, int dcomp, int ncomp,
                    const IntVect& nghost);

protected:

    Vector<IndexType> m_ixTypes;

};

}

