#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_BaseFabArraySet.H>

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
 * 
 * \tparam N The number of MultiFabs in the set.
 */
template <std::size_t N>
class MultiFabSet
    :
    public BaseMultiFabSet<MultiFab>
{
public:

    using MF = typename amrex::MultiFab;
    using fab_type = typename amrex::FArrayBox; // MF::fab_type;
    using value_type = typename amrex::Real; // MF::value_type;

    // using RT = typename amrex::Real; // typename value_type;
    using RT = typename MF::value_type;

    using FABFactory = typename amrex::FabFactory<fab_type>;
    using DM = typename amrex::DistributionMapping;

    /**
     * \brief Make an alias MultiFabSet. maketype must be
     * 
     * \param mfs Starting MultiFabSet
     * \param maketype maketype must be amrex::make_alias
     * \param scompMF The starting MultiFab component index of the alias
     * \param ncompMF The number of components in each MultiFab the alias
     */
    MultiFabSet (const MultiFabSet<N>& mfs, MakeType maketype, int scompMF = 0, int ncompMF = 1);

    /**
     * \brief Construct a new MultiFabSet<N>
     * 
     * \param ncompMF Number of components that each MultiFab element has
     * \param ngrow 
     * \param info 
     * \param factory
     * \param ncompSet Number of components in the 
     * \param ix_type_array If not null, overrides ngrow for each MF element of the MultiFabSet.
     */
    MultiFabSet (const BoxArray& bxs, const DM& dm, int ncompMF, const IntVect& ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo& info, const FABFactory& factory,
#else
                 const MFInfo& info = MFInfo(), const FABFactory& factory = FArrayBoxFactory(),
#endif
                 const Array<IntVect,N>* ix_type_array = nullptr);

    MultiFabSet (const BoxArray& bxs, const DM& dm, int ncompMF, const Array<IntVect,N>& ngrowSet,
#ifdef AMREX_STRICT_MODE
                 const MFInfo& info, const FABFactory& factory,
#else
                 const MFInfo& info = MFInfo(), const FABFactory& factory = FArrayBoxFactory(),
#endif
                 const Array<IntVect,N>* ix_type_array = nullptr);

    /**
     * \brief Generate MultiFabSet from Array<MultiFab,N>
     */
    MultiFabSet (const Array<const MF&,N>& mfarray);

    
    MF& operator[] (int i) noexcept { return *m_mf_array[i]; }
    const MF& operator[] (int i) const noexcept { return *m_mf_array[i]; }

    MF& getElem (int i) noexcept { return *m_mf_array[i]; }
    const MF& getElem (int i) const noexcept { return *m_mf_array[i]; }

    MF* getElemPtr (int i) noexcept { return m_mf_array[i]; }
    MF const* getElemPtr (int i) const noexcept { return m_mf_array[i]; }
    
    // void setElem (int i, MF* mfptr) { m_mf_array[i] = mfptr; }
    // void setElem (int i, const MF* mfptr) { m_mf_array[i] = mfptr; }
    // void setElem (int i, MF& mf) { m_mf_array[i] = &mf; }
    // void setElem (int i, const MF& mf) { m_mf_array[i] = &mf; }

    void setElem (int i, MF *mfptr) { m_mf_array[i] = *mfptr; }
    void setElem (int i, const MF *mfptr) { m_mf_array[i] = *mfptr; }
    void setElem (int i, MF &mf) { m_mf_array[i] = mf; }
    void setElem (int i, const MF &mf) { m_mf_array[i] = mf; }


    /**
     * \brief Set constant value for selected Multifabs
     *
     * \param val The value to be set
     * \param start_index The starting index in the MultiFabSet
     * \param num_of_MF The number of MultiFabs affected starting from `start_index`
     */
    void setVal (RT val, int start_index = 0, int num_of_MF = N);
    
    /**
     * \brief Return the largest infinity norm amongst the selected MultiFabs
     *
     * \param comp           The starting index in the MultiFabSet
     * \param ncomp          The number of MultiFabs affected starting from `start_index`
     * \param nghost         number of ghost cells
     * \param local          If true, MPI communciation is skipped.
     * \param ignore_covered Ignore covered cells. Only relevant for cell-centered EB data.
     */
    RT norminf (int comp = 0, int ncomp = N, IntVect const& nghost = IntVect(0), bool local = false,
             [[maybe_unused]] bool ignore_covered = false) const;
    /**
     * \brief Perform local copy of MultiFabSet data.
     *
     * \param src source MultiFabSet
     * \param scomp starting component of source
     * \param dcomp starting component of this MultiFabSet
     * \param ncomp number of components
     * \param nghost number of ghost cells
     * \param use_nGrowVect Use nGrowVect from each MultiFab instead of nghost
     */
    void LocalCopy (const MultiFabSet<N>& src, int scomp = 0, int dcomp = 0, int ncomp = -1,
                    IntVect const& nghost = IntVect(0), bool use_nGrowVect = false);

    /**
     * \brief Used in MLLinOpT<MF>::make
     */
    void clear ();

    
    int nComp() { return N; } //Could also return nComp of the first component...

    bool isAllRegular() const noexcept;

private:
    // Array<MF*,N> m_mf_array;

    MF* m_mf_array;

};





}


    // MultiFabSet (const MultiFabSet&, MFInfo&, int, IntVect);
    // MultiFabSet (const MultiFabSet&, MFInfo, int, IxTypes);
    // MultiFabSet (const MultiFabSet&, MFInfo&, int, FabFactory<FArrayBox>&);
    // MultiFabSet (const MultiFabSet&, MFInfo, IntVect, IxTypes);
    // MultiFabSet (const MultiFabSet&, MFInfo&, IntVect, FabFactory<FArrayBox>&);
    // MultiFabSet (const MultiFabSet&, MFInfo, IntVect, IxTypes);
    // MultiFabSet (const MultiFabSet&, MFInfo);
    // MultiFabSet (const MultiFabSet&, MFInfo, FabFactory);
    // MultiFabSet (const MultiFabSet&, MFInfo, IntVect, FabFactory);