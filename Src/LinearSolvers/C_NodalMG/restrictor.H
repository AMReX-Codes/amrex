//BL_COPYRIGHT_NOTICE

#ifndef _RESTRICTOR_H_
#define _RESTRICTOR_H_

#include "interface.H"
#include "boundary.H"

//
// This is the base class for all amr_restrictors.  The basic functions
// box and fill work in the interiors of patches:  Give box the fine
// grid dimensions and it will tell you what coarse patch you can fill
// from it; fill then does the work.  The interface routine is for
// restrictors like node-based full-weighting that do special stuff
// at the edges of grids.
//

class amr_restrictor
{
public:
    virtual ~amr_restrictor () = 0;
    virtual Box box (const Box&,
		     const IntVect&) const;
    virtual Box rebox (const Box&,
		       const IntVect&) const;
    virtual void fill (FArrayBox&,
		       const Box&,
		       const FArrayBox&,
		       const IntVect&) const = 0;
    virtual void fill_interface (MultiFab&, MultiFab&,
				 const level_interface&,
				 const amr_boundary*,
				 const IntVect&) const;
};

class cell_average_restrictor
    : public amr_restrictor
{
public:
    explicit cell_average_restrictor (int i) ;
    virtual void fill (FArrayBox& patch, const Box& region,
		       const FArrayBox& fgr, const IntVect& rat) const;
protected:
    int integrate;
};

inline
cell_average_restrictor::cell_average_restrictor (int i)
    :
    integrate(i)
{
}

class terrain_velocity_restrictor
    : public cell_average_restrictor
{
public:
    explicit terrain_velocity_restrictor (int i) ;
    virtual void fill (FArrayBox& patch, const Box& region,
		       const FArrayBox& fgr, const IntVect& rat) const;
};

inline
terrain_velocity_restrictor::terrain_velocity_restrictor (int i)
    : cell_average_restrictor(i)
{
    BL_ASSERT(i >= 0 && i < BL_SPACEDIM);
}

class injection_restrictor
    : public amr_restrictor
{
public:
    virtual void fill (FArrayBox& patch, const Box& region,
		       const FArrayBox& fgr, const IntVect& rat) const;
};

class default_restrictor
    : public amr_restrictor
{
public:
    virtual void fill (FArrayBox& patch, const Box& region,
		       const FArrayBox& fgr, const IntVect& rat) const;
};

//
// Vanilla bilinear (full-weighting) restrictor class.
// Restricts fine-fine interfaces, even those at the boundary.
// Does nothing at coarse-fine interfaces.
//

class bilinear_restrictor
    : public amr_restrictor
{
public:
    bilinear_restrictor (int i,
			 bool hg_terrain);
    virtual Box box (const Box& fb,
		     const IntVect& rat) const;
    virtual Box rebox (const Box& cb,
		       const IntVect& rat) const;
    virtual void fill (FArrayBox& patch,
		       const Box& region,
		       const FArrayBox& fgr,
		       const IntVect& rat) const;
    virtual void fill_interface (MultiFab&,
				 MultiFab&,
				 const level_interface&,
				 const amr_boundary*,
				 const IntVect&) const;
protected:
    int  integrate;
    bool m_hg_dense;
};

#endif /*_RESTRICTOR_H_*/
