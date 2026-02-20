#ifndef AMREX_LayoutData_H_
#define AMREX_LayoutData_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_LayoutData.H
 * \brief Distributed container holding one object of type \c T per grid box.
 */

#include <AMReX_FabArrayBase.H>
#include <AMReX_MFIter.H>

namespace amrex
{
  /**
   * \class amrex::LayoutData
   * \ingroup amrex_fab_containers
   * \brief Distributed container holding one object of type \c T per grid box.
   *
   * Each MPI rank owns only the objects corresponding to the boxes it owns
   * according to the DistributionMapping.  Elements are accessed by iterating
   * with an \c MFIter or by global box index.
   */
  template <class T>
  class LayoutData: public FabArrayBase
  {
  public:

    //! Default-construct an empty (undefined) LayoutData.
    LayoutData() = default;

    //! Construct and define over \p a_grids with distribution \p a_dm.
    LayoutData(const BoxArray            & a_grids,
               const DistributionMapping & a_dm)
      {
        define(a_grids, a_dm);
      }

    //! Define (or redefine) over \p a_grids with distribution \p a_dm.
    void define(const BoxArray           & a_grids,
                const DistributionMapping & a_dm)
      {
        FabArrayBase::define(a_grids, a_dm, 1, 0);
        m_need_to_clear_bd = true;
        addThisBD();
        m_data.resize(local_size());
      }

    ~LayoutData () { if (m_need_to_clear_bd) { clearThisBD(); } }

    LayoutData (const LayoutData& rhs) noexcept
        : FabArrayBase(rhs),
          m_data(rhs.m_data),
          m_need_to_clear_bd(rhs.m_need_to_clear_bd)
    {
        if (m_need_to_clear_bd) { addThisBD(); }
    }

    LayoutData (LayoutData&& rhs) noexcept
        : FabArrayBase(static_cast<FabArrayBase&&>(rhs)),
          m_data(std::move(rhs.m_data)),
          m_need_to_clear_bd(rhs.m_need_to_clear_bd)
    {
        rhs.m_need_to_clear_bd = false;
    }

    LayoutData& operator= (const LayoutData& rhs) noexcept
    {
        if (&rhs != this)
        {
            m_data.clear();
            if (m_need_to_clear_bd) { clearThisBD(); }

            FabArrayBase::operator=(rhs);
            m_data = rhs.m_data;
            m_need_to_clear_bd = rhs.m_need_to_clear_bd;

            if (m_need_to_clear_bd) { addThisBD(); }
        }
        return *this;
    }

    LayoutData& operator= (LayoutData&& rhs) noexcept
    {
        if (&rhs != this)
        {
            m_data.clear();
            if (m_need_to_clear_bd) { clearThisBD(); }

            FabArrayBase::operator=(static_cast<FabArrayBase&&>(rhs));
            m_data = std::move(rhs.m_data);
            m_need_to_clear_bd = rhs.m_need_to_clear_bd;

            rhs.m_need_to_clear_bd = false;
        }
        return *this;
    }

    //! Return a reference to the element owned by the current \p a_mfi tile.
    T& operator[](const MFIter& a_mfi) noexcept
      {
        int local_index = a_mfi.LocalIndex();
        BL_ASSERT(local_index >= 0 && local_index < m_data.size() &&
                  DistributionMap() == a_mfi.DistributionMap());
        return m_data[local_index];
      }


    //! Return a const reference to the element owned by the current \p a_mfi tile.
    const T& operator[](const MFIter& a_mfi) const noexcept
      {
        int local_index = a_mfi.LocalIndex();
        BL_ASSERT(local_index >= 0 && local_index < m_data.size() &&
                  DistributionMap() == a_mfi.DistributionMap());
        return m_data[local_index];
      }

    //! Return a reference to the element for global box index \p a_box_index.
    T& operator[](int a_box_index) noexcept
      {
        int local_index = this->localindex(a_box_index);
        BL_ASSERT(local_index >= 0 && local_index < m_data.size());
        return m_data[local_index];
      }


    //! Return a const reference to the element for global box index \p a_box_index.
    const T& operator[](int a_box_index) const noexcept
      {
        int local_index = this->localindex(a_box_index);
        BL_ASSERT(local_index >= 0 && local_index < m_data.size());
        return m_data[local_index];
      }

    //! Return a const pointer to the local element array.
    const T* data () const noexcept { return m_data.data(); }
    //! Return a pointer to the local element array.
    T*       data ()       noexcept { return m_data.data(); }

  private:
      Vector<T> m_data;
      bool m_need_to_clear_bd = false;
  };
}
#endif
