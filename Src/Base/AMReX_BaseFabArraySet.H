#include <AMReX_FabArray.H>

namespace amrex {

/**
 * \brief BaseFabArraySet
 * 
 * \tparam MF FabArray or derived class
 * \tparam encodedIxTypes 
 */
// template <class MF, unsigned long encodedIxTypes>
class BaseFabArraySet {
    
public:

    BaseFabArraySet () = default;
    ~BaseFabArraySet () = default;

    BaseFabArraySet (BaseFabArraySet&& rhs) noexcept = default;
    BaseFabArraySet (const BaseFabArraySet& rhs) = default;
    BaseFabArraySet& operator= (const BaseFabArraySet& rhs) = default;
    BaseFabArraySet& operator= (BaseFabArraySet&& rhs) = default;
  
    //! Return the grow factor that defines the region of definition.
    [[nodiscard]] int nGrow (int direction = 0) const noexcept { return m_nGrowVect[direction]; }

    [[nodiscard]] IntVect nGrowVect () const noexcept { return m_nGrowVect; }

    //! Return number of variables (aka components) associated with each point.
    //! Need to decide if nComp reutrns the number of elements in the set or the number of components in each MF
    [[nodiscard]] int nComp () const noexcept { return m_nComp; }

    //Return whether this BaseFabArraySet is empty
    [[nodiscard]] bool empty () const noexcept { return m_boxarray.empty(); }

    /**
    * \brief Return a constant reference to the nodal BoxArray that defines the
    * valid region associated with this BaseFabArraySet.
    */
    [[nodiscard]] const BoxArray& boxArray () const noexcept { return m_boxarray; }

    //! Return constant reference to associated DistributionMapping.
    [[nodiscard]] const DistributionMapping& DistributionMap () const noexcept { return m_distributionMap; }

protected:

    
    BoxArray            m_boxarray;
    DistributionMapping m_distributionMap;
    IntVect             m_nGrowVect;
    int                 m_nComp;
    

};

}