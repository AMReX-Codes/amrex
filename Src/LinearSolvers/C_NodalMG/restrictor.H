#ifndef _RESTRICTOR_H_
#define _RESTRICTOR_H_

#include "interface.H"
#include "boundary.H"

// This is the base class for all amr_restrictors.  The basic functions
// box and fill work in the interiors of patches:  Give box the fine
// grid dimensions and it will tell you what coarse patch you can fill
// from it; fill then does the work.  The interface routine is for
// restrictors like node-based full-weighting that do special stuff
// at the edges of grids.

class amr_restrictor_class 
{
public:
    amr_restrictor_class() 
    { 
    }
    virtual ~amr_restrictor_class() 
    { 
    }
    virtual Box box(const Box&, const IntVect&) const;
    virtual Box rebox(const Box&, const IntVect&) const;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const IntVect&) const = 0;
    virtual void fill_interface(MultiFab&, MultiFab&, const level_interface&, const amr_boundary_class*, const IntVect&) const 
    {
    }
};

class cell_average_restrictor_class : public amr_restrictor_class 
{
public:
    explicit cell_average_restrictor_class(int i) 
	: integrate(i) 
    { 
    }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
protected:
    int integrate;
};

class terrain_velocity_restrictor_class : public cell_average_restrictor_class 
{
public:
    explicit terrain_velocity_restrictor_class(int i) 
	: cell_average_restrictor_class(i) 
    { 
    }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
};

class injection_restrictor_class : public amr_restrictor_class 
{
public:
    injection_restrictor_class() 
    { 
    }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
};

class default_restrictor : public amr_restrictor_class 
{
public:
    default_restrictor() 
    { 
    }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
};

// Vanilla bilinear (full-weighting) restrictor class.
// Restricts fine-fine interfaces, even those at the boundary.
// Does nothing at coarse-fine interfaces.
class bilinear_restrictor_class : public amr_restrictor_class 
{
public:
    bilinear_restrictor_class(int i, bool hg_terrain);
    virtual Box box(const Box& fb, const IntVect& rat) const;
    virtual Box rebox(const Box& cb, const IntVect& rat) const;
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
    virtual void fill_interface(MultiFab&, MultiFab&, const level_interface&, const amr_boundary_class*, const IntVect&) const;
protected:
    int integrate;
    bool m_hg_terrain;
};

#endif
