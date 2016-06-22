
#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

#include <InterpBndryData.H>

class MacBndry
    :
    public InterpBndryData
{
public:

    MacBndry ();

    MacBndry (const BoxArray& _grids,
              int             _ncomp,
              const Geometry& _geom,
	      ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());

    virtual ~MacBndry ()  BL_OVERRIDE;

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio) BL_OVERRIDE;

    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
                                int            comp=0) BL_OVERRIDE;
    void setHomogValues (const BCRec&   bc,
                         const IntVect& ratio);

    int phys_bc_lo (int dir) const {return m_phys_bc.lo(dir);}
    int phys_bc_hi (int dir) const {return m_phys_bc.hi(dir);}

private:

    // save phys_bc for FMultiGrid
    BCRec m_phys_bc;
};

#endif
