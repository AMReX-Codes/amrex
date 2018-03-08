#ifndef _EBNORMALIZEBYVOLUMEFRACTION_H___
#define _EBNORMALIZEBYVOLUMEFRACTION_H___

#include "AMReX_EBCellFAB.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_LayoutData.H"
#include "AMReX_FabArray.H"
#include "AMReX_AggStencil.H"


namespace amrex
{

///
/** 
 *  This class represents an operation that takes a volume-fraction-weighted
 *  quantity kappa*Q (where kappa is the volume fraction) and computes Q
 *  from it by normalizing it w.r.t. the volume fractions of cells within a
 *  monotone path radius. The radius may either be specified explicitly or
 *  computed from a given EBLevelGrid.
 */
  class EBNormalizeByVolumeFraction
  {
  public:

    ///
    /** Construct an EBNormalizeByVolumeFraction object associated with the
     *  associated with the given EBLevelGrid.
     *  \param EBLevelGrid The level grid for which data can be normalized
     *           by the resulting object.
     *  nghost is the number of ghost cells in data
     */
    EBNormalizeByVolumeFraction(const EBLevelGrid          & a_eblg,
                                const int                  & a_nghost,
                                const int                    a_radius = 1)
      {
        define(a_eblg, a_nghost, a_radius);
      }


    ///
    void define(const EBLevelGrid          & a_eblg,
                const int                  & a_ghost,
                const int                    a_radius = 1);
  
    /** Destructor. */
    ~EBNormalizeByVolumeFraction()
      {
      }
  
    EBNormalizeByVolumeFraction()
      {
        m_isDefined = false;
      }


    ///set Q = sum(kappa Q)/sum(kappa)
    void normalize(FabArray<EBCellFAB>&       a_Q,
                   const FabArray<EBCellFAB>& a_kappaQ,
                   int idst, int inco);

  protected:

    void getLocalStencil(VoFStencil      & a_stencil, 
                         const VolIndex  & a_vof, 
                         const MFIter    & a_dit);

    EBNormalizeByVolumeFraction(const EBNormalizeByVolumeFraction&);
    EBNormalizeByVolumeFraction& operator=(const EBNormalizeByVolumeFraction&);

  private:
    bool m_isDefined;
    int m_radius;
    int m_ghost;
    LayoutData<std::shared_ptr<AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
    EBLevelGrid  m_eblg;
  };

}
#endif
