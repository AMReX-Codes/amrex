/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_EBArith.H"
#include "AMReX_EBFastFR.H"

/*******************/
EBFastFR::
EBFastFR(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const int&         a_refRat,
         const int&         a_nvar,
         bool a_forceNoEBCF)
{
  setDefaultValues();
  define(a_eblgFine, a_eblgCoar, a_refRat, a_nvar, a_forceNoEBCF);
}
/*******************/
bool
EBFastFR::computeHasEBCF()
{
  CH_TIME("EBFastFR::computeHasEBCF");
  const ProblemDomain&            domai =   m_eblgFine.getDomain();
  const EBISLayout&               ebisl =   m_eblgFine.getEBISL();
  const DisjointBoxLayout&        grids =   m_eblgFine.getDBL();
  const LayoutData<IntVectSet>&   cfivs = *(m_eblgFine.getCFIVS());
  bool doEBCrossing = false;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          Box grid  = grids.get(dit());
          IntVect ivgrow = IntVect::Zero;
          IntVectSet ebcfivsLo, ebcfivsHi;
          EBCFData::getEBCFIVSGrid(ebcfivsLo, grid, idir, Side::Lo, ivgrow, domai, cfivs[dit()], ebisl[dit()]);
          EBCFData::getEBCFIVSGrid(ebcfivsHi, grid, idir, Side::Hi, ivgrow, domai, cfivs[dit()], ebisl[dit()]);

          //keep track if this  crossing ever happens
          if ((!ebcfivsLo.isEmpty()) ||
             (!ebcfivsHi.isEmpty()))
            {
              doEBCrossing = true;
            }
        }
    }

  //in the case of parallel, need to check if ANY of the procs
  //have ebcrossing
#ifdef CH_MPI
  int gatherint = 0;
  if (doEBCrossing) gatherint = 1;
  //  pout() << "before gather = " << gatherint << endl;

  int idoebcf;
  MPI_Allreduce(&gatherint, &idoebcf, 1, MPI_INT,
                MPI_MAX, MPI_COMM_WORLD);
  doEBCrossing = (idoebcf==1);


#endif

  return doEBCrossing;
}
/*******************/
void
EBFastFR::
define(const EBLevelGrid& a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_refRat,
       const int&         a_nvar,
       bool a_forceNoEBCF)
{
  CH_TIME("EBFastFR::define");
  clear();
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (a_forceNoEBCF)
    {
      m_hasEBCF = false;
    }
  else
    {
      m_hasEBCF = computeHasEBCF();
    }
  if(!m_hasEBCF)
  {
    amrex::Error("EB/CF bit not ported yet for flux register");
  }
  if (!m_eblgFine.coarsenable(a_refRat))
    {
      MayDay::Error("EBFastFR::define: dbl not coarsenable by refrat");
    }

  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());

#if (CH_SPACEDIM == 2)
  m_nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  m_nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif

  //set all registers to zero to start.
  setToZero();
  m_isDefined = true;
}

/*********/
EBFastFR::
EBFastFR()
{
  setDefaultValues();
}
/*********/
void
EBFastFR::
clear()
{
  m_isDefined = false;
}
/*******************/
EBFastFR::~EBFastFR()
{
  clear();
}
/*******************/
bool
EBFastFR::
isDefined() const
{
  return m_isDefined;
}
/*******************/
void
EBFastFR::
setToZero()
{
  CH_TIME("EBFastFR::setToZero");
  m_levelFluxReg.setToZero();

}
/*******************/
void
EBFastFR::
incrementCoarse(const EBFaceFAB&      a_coarFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_coarDatInd,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd)
{
  incrementCoarRegul(a_coarFlux,
                     a_scale,
                     a_coarDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

  incrementCoarIrreg(a_coarFlux,
                     a_scale,
                     a_coarDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

}
/*************/
void
EBFastFR::
incrementCoarRegul(const EBFaceFAB& a_coarFlux,
                   const Real&      a_scale,
                   const DataIndex& a_coarDatInd,
                   const Interval&  a_variables,
                   const int&       a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_assert(m_isDefined);
  CH_TIME("EBFastFR::incrementCoarseBoth");
  FArrayBox& coarFluxFAB = (FArrayBox&)(a_coarFlux.getFArrayBox());

  //increment  as if there were no EB
  m_levelFluxReg->incrementCoarse(coarFluxFAB, a_scale, a_coarDatInd,
                                  a_variables, a_variables, a_dir, a_sd);
}
/*******************/
void
EBFastFR::
incrementCoarIrreg(const EBFaceFAB&       a_coarFlux,
                   const Real&            a_scale,
                   const DataIndex&       a_coarDatInd,
                   const Interval&        a_variables,
                   const int&             a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementCoarseIrreg");
  CH_assert(m_isDefined);
  if (m_hasEBCF)
    {

      const EBISBox& ebisBox     = m_eblgCoar.getEBISL()[a_coarDatInd];
      int iindex = index(a_dir, a_sd);

      Vector<VoFIterator>& vofits = (m_vofiCoar[iindex])[a_coarDatInd];

      for (int ivofits = 0; ivofits < vofits.size(); ivofits++)
        {
          VoFIterator& vofit = vofits[ivofits];

          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //we are increment the coarse side so we have to get the faces on the
              //flip of sd
              Vector<FaceIndex> faces = ebisBox.getFaces(vof, a_dir, flip(a_sd));
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBox.areaFrac(face);
                      Real flux = a_coarFlux(face, ivar);
                      Real areaSum = area*flux;
                      //consistent with non-eb version
                      //              Real oldFlux = m_delUCoar[a_coarDatInd][a_dir](face, ivar);
                      areaSum *= -a_scale;
                      m_delUCoar[a_coarDatInd](vof, ivar) +=  areaSum;
                    } //end loop over variables
                }//end is this face in the holders
            }//loop over vofs
        }//Be bloody, bold and resolute
    } //for no man of woman born shall harm macbeth
}
/*******************/
void
EBFastFR::
incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                  const Real&           a_scale,
                  const DataIndex&      a_fineDatInd,
                  const Interval&       a_variables,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineBoth");
  incrementFineRegul(a_fineFlux,
                     a_scale,
                     a_fineDatInd,
                     a_variables,
                     a_dir,
                     a_sd);


  incrementFineIrreg(a_fineFlux,
                     a_scale,
                     a_fineDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

}
/*******************/
void
EBFastFR::
incrementFineRegul(const EBFaceFAB&      a_fineFlux,
                   const Real&           a_scale,
                   const DataIndex&      a_fineDatInd,
                   const Interval&       a_variables,
                   const int&            a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineRegul");
  //increment  as if there were no EB
  FArrayBox& fineFluxFAB = (FArrayBox&)(a_fineFlux.getFArrayBox());
  m_levelFluxReg->incrementFine(fineFluxFAB, a_scale, a_fineDatInd,
                                a_variables, a_variables, a_dir, a_sd);

  incrementFineSparse(a_fineFlux, a_scale, a_fineDatInd, a_variables,
                      a_dir, a_sd, true);
}
/*******************/
void
EBFastFR::
incrementFineIrreg(const EBFaceFAB&      a_fineFlux,
                   const Real&           a_scale,
                   const DataIndex&      a_fineDatInd,
                   const Interval&       a_variables,
                   const int&            a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineIrreg");
  incrementFineSparse(a_fineFlux, a_scale, a_fineDatInd, a_variables,
                      a_dir, a_sd, false);

}// mmm deep loops

void
EBFastFR::
compareFineSparse(const EBFaceFAB&      a_fluxOld,
                  const EBFaceFAB&      a_fluxNew,
                  const DataIndex&      a_fineDatInd,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineSparse");
  if (m_hasEBCF)
    {
      int iindex = index(a_dir, a_sd);

      VoFIterator& vofit       = (m_vofiCoFi[iindex])[a_fineDatInd];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          //remember the registers live at the coar level
          const VolIndex& coarVoF = vofit();
          const EBISBox& ebisBoxCoFi = m_eblgCoFi.getEBISL()[a_fineDatInd];

          Vector<FaceIndex> facesCoar = ebisBoxCoFi.getFaces(coarVoF, a_dir, flip(a_sd));
          for (int ifacec = 0; ifacec < facesCoar.size(); ifacec++)
            {
              Vector<FaceIndex> facesFine =
                m_eblgCoFi.getEBISL().refine(facesCoar[ifacec], m_refRat, a_fineDatInd);
              for (int ifacef = 0; ifacef < facesFine.size(); ifacef++)
                {
                  VolIndex vofFine = facesFine[ifacef].getVoF(flip(a_sd));
                  for (int ivar = 0; ivar < a_fluxOld.nComp(); ivar++)
                    {
                      Real fluxOld  = a_fluxOld(facesFine[ifacef], ivar);
                      Real fluxNew  = a_fluxNew(facesFine[ifacef], ivar);
                      Real eps = 1.0e-9;
                      if (Abs(fluxOld-fluxNew) > eps)
                        {
                          pout() << "hurm at fine face" << facesFine[ifacef] << endl;
                        }
                    } //this is the way the world ends
                }//this is the way the world ends
            }//this is the way the world ends
        }//not with a bang
    } //but with a
}//whimper
/*******************/
void
EBFastFR::
incrementFineSparse(const EBFaceFAB&      a_fineFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_fineDatInd,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd,
                    bool a_doingFineRegular)
{
  CH_TIME("EBFastFR::incrementFineSparse");
  if (m_hasEBCF)
    {
      //sign because of which side of the
      //vof we are on is opposite of which side
      // of the fine grid we are on
      //int isign = -sign(a_sd);
      //Real rsign = isign;
      Real newScale = a_scale/m_nrefdmo;
      int iindex = index(a_dir, a_sd);

      const EBISBox& ebisBoxFine= m_eblgFine.getEBISL() [a_fineDatInd];
      VoFIterator& vofit       = (m_vofiCoFi[iindex])[a_fineDatInd];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          //remember the registers live at the coar level
          const VolIndex& coarVoF = vofit();
          const EBISBox& ebisBoxCoFi = m_eblgCoFi.getEBISL()[a_fineDatInd];

          EBCellFAB& delUCoFi = m_delUCoFi[a_fineDatInd];
          //we are on the coarse side of the face so we have to look back on it---ergo flip
          Vector<FaceIndex> facesCoar = ebisBoxCoFi.getFaces(coarVoF, a_dir, flip(a_sd));
          for (int ifacec = 0; ifacec < facesCoar.size(); ifacec++)
            {
              Vector<FaceIndex> facesFine =
                m_eblgCoFi.getEBISL().refine(facesCoar[ifacec], m_refRat, a_fineDatInd);
              for (int ifacef = 0; ifacef < facesFine.size(); ifacef++)
                {
                  Real area = ebisBoxFine.areaFrac(facesFine[ifacef]);
                  VolIndex vofFine = facesFine[ifacef].getVoF(flip(a_sd));
                  bool vofIsRegular = ebisBoxFine.isRegular(vofFine.gridIndex());
                  if (( a_doingFineRegular &&  vofIsRegular) || (!a_doingFineRegular && !vofIsRegular))
                    {
                      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                        {
                          Real flux          = a_fineFlux(facesFine[ifacef], ivar);
                          delUCoFi(coarVoF, ivar) +=  newScale*area*flux;
                        } //end loop over variables
                    }
                } //end loop over fine faces
            }// end loop over coarse faces
        }//end loop over EBCF vofs
    } //end if we have EBCF in this box
}// mmm deep loops
/*******************/
void
EBFastFR::
reflux(LevelData<EBCellFAB>& a_uCoar,
       const Interval&       a_variables,
       const Real&           a_scale,
       bool a_multByKappaOneMinusKappa)
{
  CH_TIME("EBFastFR::reflux");
  LevelData<FArrayBox> uCoarLDF;
  //save initial values of ucoar because the non-eb flux  reg will
  //change it in its ignorance
  if (m_hasEBCF)
    {
      cacheOldSolution(a_uCoar, a_variables);
    }

  //reflux as if there were no EB
  aliasEB(uCoarLDF, a_uCoar);
  m_levelFluxReg->reflux(uCoarLDF, a_variables, a_variables, a_scale);

  //correct at irregular cells
  if (m_hasEBCF)
    {
      restoreOldSolution(a_uCoar, a_variables);

      irregReflux(a_uCoar, a_variables, a_scale, a_multByKappaOneMinusKappa);
    }
}
/*******************/
void
EBFastFR::
irregReflux(LevelData<EBCellFAB>& a_uCoar,
            const Interval&       a_variables,
            const Real&           a_scale,
            bool a_multByKappaOneMinusKappa)
{
  CH_assert(m_isDefined);
  CH_TIME("EBFastFR::irregReflux");
  if (m_hasEBCF)
    {
      //coming into this routine, coar holds -coarFlux*area
      //and fine holds area*fineFlux
      //make diff = area(fineFlux-coarFlux)  (with scaling stuff)

      EBLevelDataOps::setVal(m_delUDiff, 0.0);
      m_delUCoar.copyTo(a_variables, m_delUDiff, a_variables);

      EBAddOpEBFFR op;
      m_delUCoFi.copyTo(a_variables, m_delUDiff, a_variables, m_reverseCopier, op);

      //add refluxDivergence to solution u -= a_scale*(area*(fineFlux-coarFlux))
      incrementByRefluxDivergence(a_uCoar, m_delUDiff, a_variables, a_scale, false,
                                  a_multByKappaOneMinusKappa);

      //dumpEBFFR(a_uCoar, string("uCoar holds after reflux divergence"));

    }
}
/*******************/
void
EBFastFR::
incrementByRefluxDivergence(LevelData<EBCellFAB>& a_uCoar,
                            LevelData<EBCellFAB>& a_delUDiff,
                            const Interval      & a_variables,
                            const Real          & a_newScale,
                            bool a_multByOneMinusKappa,
                            bool a_multByKappaOneMinusKappa)
{
  if (m_hasEBCF)
    {
      //cannot both be true
      CH_assert(!a_multByOneMinusKappa  || !a_multByKappaOneMinusKappa);
      for (DataIterator dit= m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int iindex = index(idir, sit());
                  const EBISBox& ebisBox       = m_eblgCoar.getEBISL()[dit()];
                  Vector<VoFIterator> vofitVec = m_vofiCoar[iindex][dit()];
                  for (int ivec = 0; ivec < vofitVec.size(); ivec++)
                    {
                      VoFIterator& vofit = vofitVec[ivec];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          const VolIndex& vof = vofit();
                          Real scale = sign(sit())*a_newScale;
                          Real volFrac = ebisBox.volFrac(vof);
                          if (a_multByOneMinusKappa)
                            {
                              scale *= (1.0-volFrac);
                            }
                          if (a_multByKappaOneMinusKappa)
                            {
                              scale *= volFrac*(1.0-volFrac);
                            }
                          for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                            {
                              Real udelDif = a_delUDiff[dit()](vof, ivar);
                              Real soluVal = a_uCoar[dit()](vof, ivar);
                              Real soluNew = soluVal - scale*udelDif;

                              // Just used so now zero it just in case a
                              // vof is in more than one CF interface
                              a_delUDiff[dit()](vof, ivar) = 0.0;

                              a_uCoar[dit()](vof, ivar) = soluNew;
                            }
                        }
                    }
                }
            }
        }
    }

}
/*******************/
void
EBFastFR::
restoreOldSolution(LevelData<EBCellFAB>&       a_uCoar,
                   const Interval&             a_variables)
{

  CH_TIME("EBFastFR::saveOldSolution");
  for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int iindex = index(idir, sit());
              for (int ivec = 0; ivec < (m_vofiCoar[iindex])[dit()].size(); ivec++)
                {
                  VoFIterator& vofit = (m_vofiCoar[iindex])[dit()][ivec];
                  for (vofit.reset(); vofit.ok(); ++vofit)
                    {
                      for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                        {
                          a_uCoar[dit()](vofit(), icomp) =  m_saveCoar[dit()](vofit(), icomp) ;
                        }
                    }
                }
            }
        }
    }
}
/*******************/
void
EBFastFR::
cacheOldSolution(const LevelData<EBCellFAB>& a_uCoar,
                 const Interval&             a_variables)
{
  CH_TIME("EBFastFR::saveOldSolution");
  a_uCoar.copyTo(a_variables, m_saveCoar, a_variables);
}
/*******************/
void
EBFastFR::
incrementDensityArray(LevelData<EBCellFAB>& a_coarDense,
                      const Interval&       a_variables,
                      const Real&           a_scale)
{

  //only do stuff on irregular cells because on non-irregular cells 1-kappa=0
  //does not need to be fast because this is a debugging tool
  if (m_hasEBCF)
    {
      //coming into this routine, coarFlux holds -coarFlux*area
      //and fineFlux holds area*fineFlux
      //make fluxDiff = area(fineFlux-coarFlux)  (with scaling stuff)
      EBLevelDataOps::setToZero(m_delUDiff);

      m_delUCoar.copyTo(a_variables, m_delUDiff, a_variables);
      EBAddOpEBFFR op;

      m_delUCoFi.copyTo(a_variables, m_delUDiff, a_variables, m_reverseCopier, op);

      //add refluxDivergence to solution u = (1-kappa)*a_scale*(area*(fineFlux-coarFlux))
      incrementByRefluxDivergence(a_coarDense, m_delUDiff, a_variables, a_scale, true, false);
    }
}
/*************************/
#include "NamespaceFooter.H"
