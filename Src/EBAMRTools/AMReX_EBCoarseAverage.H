#ifndef _EBCOARSEAVERAGE_H___
#define _EBCOARSEAVERAGE_H___

#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_FabArray.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_LayoutData.H"
#include "AMReX_AggStencil.H"


namespace amrex
{
/// replaces coarse level data with an average of fine level data.
/**
   This class replaces data at a coarse level of refinement with an
   average of data at a finer level of refinement, in areas where fine
   data is present.  Coarse level data is not modified where fine
   level data is not present.
   For the discussion below, kappa is the volume fraction.
   This can be used for both data where the data is orignally kappa-weighted or not.
   If you are going to use this for kappa weighted data 
   (say you are averaging down kappa*residual in multigrid), then  the kappa weight is not included in the stencil.
   If you are averaging down a non-kappa weighted quantity (density, for example), then you should turn kappa weighting on.
   This makes phi_c = (sum phi_f \kappa_f )/(sum kappa_f).   If you turn kappa weighing off (the input data is already psi_f = kappa phi_f), 
   psi_c = sum(psi_c)/(R^D) where R is the refinement ratio and D is the dimension (R^D = the number of fine cells per coarse cell).
*/
  class EBCoarseAverage
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    EBCoarseAverage()
      {
        m_isDefined = false;
      }

    ///
    ~EBCoarseAverage()
      {
      }


    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().

       {\bf Arguments:}\\
       eblgFine, eblgCoar: The fine and coarse EBLevelGrids \\
       nref: The refinement ratio between the two levels. \\
       nvar: The number of variables contained in the data
       at each VoF.
       
       IF THE DATA IS ALREADY KAPPA-WEIGHTED (think kappa*residual), SET a_useKappaWeightingInStencil to FALSE
       IF THE DATA IS *NOT* ALREADY KAPPA-WEIGHTED (think density),  SET a_useKappaWeightingInStencil to TRUE
    */
    EBCoarseAverage(const EBLevelGrid & a_eblgFine,
                    const EBLevelGrid & a_eblgCoar,
                    const int         & a_refinementRatio,
                    const int         & a_ghostCellsInData,
                    const bool        & a_useKappaWeightingInStencil,
                    bool enableFaceAveraging = false);


    ///
    /**
       Defines this object.  Existing information is overriden.

       {\bf Arguments:}\\
       dblFine, dblCoar: The fine and coarse layouts
       of the data.\\
       ebislFine, ebislCoar: The fine and coarse layouts
       of the geometric description.\\
       nref: The refinement ratio between the two levels. \\
       nvar: The number of variables contained in the data
       at each VoF.
       IF THE DATA IS ALREADY KAPPA-WEIGHTED (think kappa*residual), SET a_useKappaWeightingInStencil to FALSE
       IF THE DATA IS *NOT* ALREADY KAPPA-WEIGHTED (think density),  SET a_useKappaWeightingInStencil to TRUE
    */
    void define(const EBLevelGrid & a_eblgFine,
                const EBLevelGrid & a_eblgCoar,
                const int         & a_refinementRatio,
                const int         & a_ghost,
                const bool        & a_useKappaWeightingInStencil,
                bool enableFaceAveraging = false);


    ///
    /**
       Returns true if this object was created with the defining
       constructor or if define() has been called.
    */
    bool isDefined() const
      {
        return m_isDefined;
      }

    ///
    /**
       Replaces a_coarse_data with the average of a_fine_data within the
       coarsening of the domain of the fine level.  Elsewhere, the
       coarse data is unchanged.  It is an error to call if not
       this->isDefined(). The domain of a_fine_data should be
       the same as the fine domain specified in the most recent call to
       define().  It is expected that the coarse and fine level's
       domains are properly nested.  Both a_coarse_data and a_fine_data
       should have the same number of components specified in the most
       recent call to define().

       {\bf Arguments:}\\
       coarData: The data over the coarse layout.\\
       fineData: The data over the fine layout.  \\
       Fine and coarse data must
       have the same number of variables.\\
       variables:  The variables to average.  Those not
       in this range will be left alone.   This range of variables
       must be in both the coarse and fine data.
    */
    void
    average(FabArray<EBCellFAB>&       a_coarData,
            const FabArray<EBCellFAB>& a_fineData,
            int isrc, int idst , int inco);


    ///same deal but for faces
    void
    average(FabArray<EBFluxFAB>&       a_coarData,
            const FabArray<EBFluxFAB>& a_fineData,
            int isrc, int idst , int inco);

  protected:
    void defineStencils();
    void definePointStencil(VoFStencil& a_sten, const VolIndex& vofCoFi, const MFIter& a_mfi);
    
    //only called if face averaging enabled
    void defineFaceStencils();
    void definePointStencilFace(FaceStencil& a_sten, const FaceIndex& facCoFi, const MFIter& a_mfi);
    void definePointStencilIrreg(VoFStencil& a_sten, const VolIndex&  vofCoFi, const MFIter& a_mfi, const IntVectSet& a_ivsIrreg);

    bool m_enableFaceAveraging;
    EBLevelGrid  m_eblgFine;
    EBLevelGrid  m_eblgCoar;
    EBLevelGrid  m_eblgCoFi;
    bool         m_useKappaWeightingInStencil;

    bool m_isDefined;

    LayoutData<std::shared_ptr<AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;

    //only defined if face averaging enabled
    LayoutData<std::shared_ptr<AggStencil<EBFaceFAB, EBFaceFAB> > > m_faceStencil[BL_SPACEDIM];

    //only defined if face averaging enabled
    LayoutData<std::shared_ptr<AggStencil<IrregFAB,  IrregFAB > > > m_irregStencil;

    int m_refRat;
    int m_ghost;

  private:
    //disallowed for all the usual reasons
    EBCoarseAverage(const EBCoarseAverage& ebcin);
    void operator=(const EBCoarseAverage& fabin);

  };
}

#endif
