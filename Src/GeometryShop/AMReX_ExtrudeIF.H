#ifndef AMREX_EXTRUDEIF_H_
#define AMREX_EXTRUDEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{


  ///
  /**
     This implicit function takes 2D implicit functions (in x,y) and uses them
     to produce an extruded implicit function along z direction.
     In 2D, the implicit function is maintained as it is.
    
  */
  class ExtrudeIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying one implicit function to be extruded and whether
       the domain is on the inside (a_inside), i.e. where the function is
       negative.
    */
    ExtrudeIF(const BaseIF& a_impFunc1,
            const bool&   a_inside);
    
    /// Copy constructor

    
    /// Destructor
    virtual ~ExtrudeIF();
    
    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;
    
    
    virtual BaseIF* newImplicitFunction() const;
  protected:
    BaseIF* m_impFunc1; // implicit function to extrude
    
    bool  m_inside;   // inside flag
    
  private:
    ExtrudeIF(const ExtrudeIF& a_inputIF);    
    ExtrudeIF();
    void operator=(const ExtrudeIF& a_inputIF);

  };
}
    
#endif
