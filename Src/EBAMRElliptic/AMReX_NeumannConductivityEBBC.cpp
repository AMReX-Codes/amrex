#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "NeumannConductivityEBBC.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"
/*****************/
void NeumannConductivityEBBC::getEBFlux(Real&                         a_flux,
                                        const VolIndex&               a_vof,
                                        const LevelData<EBCellFAB>&   a_phi,
                                        const LayoutData<IntVectSet>& a_cfivs,
                                        const DataIndex&              a_dit,
                                        const RealVect&               a_probLo,
                                        const RealVect&               a_dx,
                                        const bool&                   a_useHomogeneous,
                                        const Real&                   a_time,
                                        const pair<int,Real>*         a_cacheHint )
{
  m_bc.getEBFlux(a_flux,
                 a_vof,
                 a_phi,
                 a_cfivs,
                 a_dit,
                 a_probLo,
                 a_dx,
                 a_useHomogeneous,
                 a_time,
                 a_cacheHint );

  Real bcoef = (*m_bcoe)[a_dit](a_vof,0);
  a_flux *= bcoef;
}
/*****************/
NeumannConductivityEBBC::NeumannConductivityEBBC(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  :m_bc(a_domain, a_layout, a_dx)
{
  m_dataBased = false;
}
/*****************/
NeumannConductivityEBBC::~NeumannConductivityEBBC()
{
}
/*****************/
void NeumannConductivityEBBC::setValue(Real a_value)
{
  m_bc.setValue(a_value);
}
/*****************/
void NeumannConductivityEBBC::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_bc.setFunction(a_flux);
}
/*****************/
void NeumannConductivityEBBC::applyEBFlux(EBCellFAB&                    a_lphi,
                                          const EBCellFAB&              a_phi,
                                          VoFIterator&                  a_vofit,
                                          const LayoutData<IntVectSet>& a_cfivs,
                                          const DataIndex&              a_dit,
                                          const RealVect&               a_probLo,
                                          const RealVect&               a_dx,
                                          const Real&                   a_factor,
                                          const bool&                   a_useHomogeneous,
                                          const Real&                   a_time)
{
  CH_TIME("NeumannConductivityEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp()  == 1);

  Real flux = 0.0;

  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      if (m_dataBased)
        {
//          if ((*m_data)[a_dit].getIVS().contains(vof.gridIndex()))
//             {
//               flux = (*m_data)[a_dit](vof, 0);
//             }
//          else
//            {
//              flux = 0.;
//            }
               flux = (*m_data)[a_dit](vof, 0);
        }
      else if (m_bc.m_isFunction)
        {
          const RealVect& centroid = ebisBox.bndryCentroid(vof);
          const RealVect&   normal = ebisBox.normal(vof);

          Real value = m_bc.m_flux->value(vof,centroid,normal,a_dx,a_probLo,a_dit,a_time,0);
          flux = -value;
        }
      else
        {
          if (m_bc.m_onlyHomogeneous)
            {
              MayDay::Error("NeumannConductivityEBBC::getFaceFlux called with undefined inhomogeneous BC");
            }

          flux = m_bc.m_value;
        }

      const Real& areaFrac = ebisBox.bndryArea(vof);
      flux *= areaFrac;
      Real bcoef = (*m_bcoe)[a_dit](vof,0);
      flux *= m_beta*bcoef;
      a_lphi(vof,0) += flux * a_factor;
    }
}
/*****************/
NeumannConductivityEBBCFactory::NeumannConductivityEBBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();
  m_dataBased = false;
  m_onlyHomogeneous = true;
  m_isFunction = false;
}
/*****************/
NeumannConductivityEBBCFactory::~NeumannConductivityEBBCFactory()
{
}
/*****************/
void NeumannConductivityEBBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/*****************/
void NeumannConductivityEBBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/*****************/
NeumannConductivityEBBC* NeumannConductivityEBBCFactory::create(const ProblemDomain& a_domain,
                                                                const EBISLayout&    a_layout,
                                                                const RealVect&      a_dx,
                                                                const IntVect*       a_ghostCellsPhi /*=0*/,
                                                                const IntVect*       a_ghostCellsRhs /*=0*/)
{
  NeumannConductivityEBBC* fresh = new NeumannConductivityEBBC(a_domain,a_layout,a_dx);

  if (!m_onlyHomogeneous)
    {
      if (m_dataBased)
        {
          fresh->setData(m_data);
        }
      else if (!m_isFunction)
        {
          fresh->setValue(m_value);
        }
      else
        {
          fresh->setFunction(m_flux);
        }
    }

  return fresh;
}
#include "NamespaceFooter.H"
