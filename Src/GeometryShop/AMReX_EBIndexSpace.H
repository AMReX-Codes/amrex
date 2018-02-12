
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


#ifndef _EBINDEXSPACE_H_
#define _EBINDEXSPACE_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_EBISBox.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_GeometryService.H"
#include "AMReX_EBISLevel.H"

namespace amrex
{

  class EBISLevel;
///
/**
   EBIndexSpace represents the geometric information
   of the domain.    It should be generated through
   the AMReX_EBIS class.  This follows the singleton
   pattern.
*/
  class EBIndexSpace
  {
  public:

    ///
    /**
       If a_ncellMax is set, that is the max width of
       an internal grid.  Otherwise use defaults
       of (16 in 3D, 64 in 2d)
    */
    void
    define(const Box             & a_domain,
           const RealVect        & a_origin,
           const Real            & a_dx,
           const GeometryService & a_geoserver,
           int                     a_nCellMax = -1,
           int                     a_maxCoarsenings = -1);


    void buildFirstLevel(const Box& a_domain,
                         const RealVect& a_origin,
                         const Real& a_dx,
                         const GeometryService& a_geoserver,
                         int a_nCellMax,
                         int a_maxCoarsenings);

    void buildNextLevel(const GeometryService & a_geoserver, 
                        const int             & a_whichlev);



    ///
    ~EBIndexSpace()
      {
        clear();
      }


    ///
    /**
       return level index of domain.
       return -1 if a_domain does not
       correspond to any refinement of EBIS.  public for sanity checks
    */
    int getLevel(const Box& a_domain) const;


    ///
    /**
       Return the finest level number and domain that contains a non-zero number of multivalued cells.
       Returns an invalid box and -1 if no level has multivalued cells
    */
    void getFinestLevelWithMultivaluedCells(Box& a_domain, int& a_levelNumber) const;

    int getNCellMax() const
      {
        return m_nCellMax;
      }

    ///
    void fillEBISLayout(EBISLayout     & a_ebis,
                        const BoxArray & a_grids,
                        const DistributionMapping & a_dm,
                        const Box      & a_domain,
                        const int      & a_nghost) const;

    ///
    /**
       Return true if the define function has been called.
    */
    bool isDefined() const
      {
        return m_isDefined;
      }

    ///
    void clear();

    //only AMReX_EBIS can make one
    friend class AMReX_EBIS;

    ///write to file
    void write(const string& a_filename) const;

    ///define from file
    void  read(const string& a_filename);
    ///write to file
    void writeHeader(const string& a_filename) const;

    ///define from file
    void  readHeader(const string& a_filename);

    //for testing
    int getNumLevels() const
      {
        return m_nlevels;
      }
    //for testing
    const Vector<Box>& getDomains() const
      {
        return m_domainLevel;
      }

    const EBISLevel& getEBISLevel (int lev) const { return *m_ebisLevel[lev]; }

    ///fills 0th component with value from implit function
    void fillNodeFarrayBoxFromImplicitFunction(FArrayBox& a_fab, const RealVect& a_dx, RealVect a_origin= RealVect::Zero) const;

  private:

    //only AMReX_EBIS can make one
    EBIndexSpace() = default;

    int  m_nCellMax;
    bool m_isDefined;

    Vector<std::unique_ptr<EBISLevel> > m_ebisLevel;
    //only defined if the incoming geometryservice has one.
    std::shared_ptr<const BaseIF> m_implicitFunction;
    Vector<Box> m_domainLevel;

    int m_nlevels;

    void operator=(const EBIndexSpace& ebiin) = delete;
    EBIndexSpace(const EBIndexSpace& ebiin) = delete;
  };

}

#endif
