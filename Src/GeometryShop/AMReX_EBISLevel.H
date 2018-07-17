#ifndef AMREX_EBISLEVEL_H_
#define AMREX_EBISLEVEL_H_
#include <map>
#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BoxArray.H"
#include "AMReX_FabArray.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBData.H"
#include "AMReX_EBISBox.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_GeometryService.H"
#include "AMReX_LayoutData.H"
#include "AMReX_GeomIntersectUtils.H"

namespace amrex
{

  extern void null_deleter_fab_ebg(FabArray<EBGraph>* a_ptr);

  class EBIndexSpace;

///singleton version of EBIndexSpace.
  class AMReX_EBIS
  {
  private:
    ///
    /**
       This is NULL until the first time
       instance() is called.
       On the first time it is called it
       is set to  new EBIndexSpace();
       After that it points to the same
       thing.
    */
    static EBIndexSpace* s_instance;

    AMReX_EBIS()
      {;}

    ~AMReX_EBIS()
      {;}

  public:
    ///
    /**
       Return a pointer to the singleton
       EBIndexSpace.  Until this is called
       once, the singleton is null.  See
       Design patterns for why this is done
       this way.
    */
    static EBIndexSpace* instance();

    // delete the singleton
    static void reset ();

  };

  ///
  /**
     EBISLevel represents the geometric information
     of the domain at a particular level of refinement.
  */
  class EBISLevel
  {
  public:

    EBISLevel() = default;

    ///
    /**
       Defines the level.  Calls the geoserver
       to fill the layout with geometric information.
    */
    EBISLevel(const Box             & a_domain,
              const RealVect        & a_origin,
              const Real            & a_dx,
              const int             & a_ncellMax,
              const GeometryService & a_geoserver);
    ///
    /**
       Defines the level from a finer level
       (factor of two refinement ratio).
       Fills in all the coarse-fine information
       in both levels.
    */
    EBISLevel(EBISLevel             & a_finerLevel,
              const GeometryService & a_geoserver);



    ///
    void fillEBISLayout(EBISLayout     &         a_ebis,
                        const BoxArray &         a_grids,
                        const DistributionMapping & a_dm,
                        const int      &         a_nghost) const;

    ~EBISLevel()
      {
      }

    const Box& getDomain() const
      {
        return m_domain;
      }

    ///return true if this refinement contains multivalued cells.
    /**
       if this has not been called before, this generates a blocking communication
     */
    bool hasMultiValuedCells() const
      {
        if(!m_alreadyCheckedForMVCells)
        {
          checkForMultiValuedCells();
        }
        return m_hasMultiValuedCells;
      }
    
    ///
    void checkForMultiValuedCells() const;

    ///
    void fixRegularNextToMultiValued();
         
    ///write to file
    void write(const string& a_dirname) const;

    ///define from file
    void  read(const string& a_dirname);


    ///write to file
    void writeHeader(const string& a_dirname) const;

    ///define from file
    void  readHeader(const string& a_dirname);

    const FabArray<EBGraph>& getGraph () const { return m_graph; }

    ///get dx
    Real getDx () const { return m_dx; }

  protected:

    void buildEBSurface(const GeometryService & a_geoserver);

  private:
    int m_nCellMax;
     
    ///
    BoxArray m_grids;
     
    DistributionMapping m_dm;
    ///
    Box m_domain;
     
    ///
    RealVect m_origin;
     
    ///
    Real m_dx;
     
    ///
    FabArray<EBGraph> m_graph;
     
    ///
    FabArray<EBData>  m_data;

    ///
    LayoutData<NodeMap> m_intersections;

    ///
    bool m_build_eb_surface;

    bool m_hasMoments;

    //only generated if needed (ergo mutable).
    mutable bool m_hasMultiValuedCells;

    mutable bool m_alreadyCheckedForMVCells = false;
    void defineFromGeometryService(const GeometryService &  a_geoserver);
    void coarsenVoFsAndFaces(EBISLevel& a_fineEBIS);

    //these are disallowed
    void operator=(const EBISLevel& ebiin) = delete;
    EBISLevel(const EBISLevel& ebiin) = delete;
  };
}
#endif
