
#ifndef AMREX_MACBNDRYDATA_H_
#define AMREX_MACBNDRYDATA_H_

#include <AMReX_InterpBndryData.H>

namespace amrex {

class MacBndry
    :
    public InterpBndryData
{
public:

    MacBndry ();

    MacBndry (const BoxArray& _grids,
	      const DistributionMapping& _dmap,
              int             _ncomp,
              const Geometry& _geom);

    virtual ~MacBndry ()  override;

    MacBndry (MacBndry&& rhs) = delete;
    MacBndry (const MacBndry& rhs) = delete;
    MacBndry& operator= (const MacBndry& rhs) = delete;
    MacBndry& operator= (MacBndry&& rhs) = delete;

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio) override;

    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
                                int            comp=0) override;
    void setHomogValues (const BCRec&   bc,
                         const IntVect& ratio);

    int phys_bc_lo (int dir) const noexcept {return m_phys_bc.lo(dir);}
    int phys_bc_hi (int dir) const noexcept {return m_phys_bc.hi(dir);}

    const BCRec& PhysBC () const noexcept { return m_phys_bc; }

private:

    //! save phys_bc for FMultiGrid
    BCRec m_phys_bc;
};

}

#endif
