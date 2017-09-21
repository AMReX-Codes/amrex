
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

#ifndef _EBFLUXFAB_H_
#define _EBFLUXFAB_H_

#include "AMReX_Box.H"
#include "AMReX_EBFaceFAB.H"
#include "AMReX_IrregFAB.H"


namespace amrex
{
/// A EBFaceFAB-like container for edge-centered fluxes
/**
   This is a class to contain face-centered fluxes on a box.
*/
  class EBFluxFAB
  {

  public:

    ///
    /**
       This stuff required by LevelData in parallel:
    */
    static bool preAllocatable()
      {
        return true; // symmetric messaging
      }
    static bool isCopyOMPSafe () { return true; }

    /// default constructor
    EBFluxFAB();

    /// destructor
    ~EBFluxFAB();


    ///invalid but necessary for fabarray to compile in certain cases
    EBFluxFAB(const Box& a_region,
              int        n = 1,
              bool       alloc = true,
              bool       shared = false)
      {
        //this strictly for checkpointing 
        m_region = a_region;
      }

    //same as above. 
    void resize (const Box& b,
                 int        N = 1)
      {
        m_region = b;
      }
    ///
    EBFluxFAB& operator+=(const EBFluxFAB& a_ebfab)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_fluxes[idir] += a_ebfab.m_fluxes[idir];
        }
        m_irrFlux += a_ebfab.m_irrFlux;
        return *this;
      }


    ///
    EBFluxFAB& operator-=(const EBFluxFAB& a_ebfab)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_fluxes[idir] -= a_ebfab.m_fluxes[idir];
        }
        return *this;
      }


    ///
    EBFluxFAB& operator*=(const EBFluxFAB& a_ebfab)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_fluxes[idir] *= a_ebfab.m_fluxes[idir];
        }
        m_irrFlux *= a_ebfab.m_irrFlux;
        return *this;
      }


    ///
    EBFluxFAB& operator/=(const EBFluxFAB& a_ebfab)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_fluxes[idir] /= a_ebfab.m_fluxes[idir];
        }
        m_irrFlux /= a_ebfab.m_irrFlux;
        return *this;
      }
  
    const Box& box() const
      {
        return getRegion();
      }

    /// {\bf access functions}

    /// number of components
    int nComp() const;

    /// returns cell-centered box which defines fluxBox
    const Box& getRegion() const;

    /// returns EBFaceFAB in direction dir
    EBFaceFAB& operator[] (const int dir);

    /// constant version
    const EBFaceFAB& operator[] (const int dir) const;

    IrregFAB& getEBFlux() 
      {
        return m_irrFlux;
      }

    const IrregFAB& getEBFlux() const
      {
        return m_irrFlux;
      }

    void setVal(const Real& val);

    ///defines this from input and copies the values
    void clone(const EBFluxFAB& a_input);

    ///
    EBFluxFAB& copy(const EBFluxFAB&  a_source,
                    const Box&        a_srcbox,
                    int               a_srccomp,
                    const Box&        a_destbox,
                    int               a_destcomp,
                    int               a_numcomp);


    EBFluxFAB& plus(const EBFluxFAB&  src,
                    const Box&        srcbox,
                    const Box&        destbox,
                    int               srccomp,
                    int               destcomp,
                    int               numcomp=1)
      {
        amrex::Error("this is not really a thing");
        return *this;
      }

    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;


    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);
    ///
    bool isDefined() const;

    ///
    EBFluxFAB(const EBISBox& a_ebisBox, const Box& bx, int n);

    ///
    void define(const EBISBox& a_ebisBox, const Box& bx, int n);


    /// Invalid but necessary for LevelData<EBFluxFAB> to compile
    EBFluxFAB(const Box& bx, int n)
      {
        amrex::Error("invalid constructor called for EBFluxFAB");
      }
    EBISBox getEBISBox() const
      {
        return m_ebisBox;
      }

  protected:

    ///
    Box m_region;

    ///
    EBISBox m_ebisBox;

    ///
    int m_nComp;

    bool m_isDefined;

    /// CH_SPACEDIM EBFaceFABes which hold fluxes
    EBFaceFAB m_fluxes[BL_SPACEDIM];
    IrregFAB  m_irrFlux;


  private:
    /// these are disallowed

    EBFluxFAB (const EBFluxFAB&);
    void operator = (const EBFluxFAB&);

  };
}
#endif
