#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "interface.H"
#include "RegType.H"

// This is the base class for all amr_boundaries.  The default function dir
// should always be acceptable.  The default box and check_against_boundary
// functions are appropriate for wall-type boundaries.  The default
// boundary_mesh function is appropriate for all boundary types I've
// created so far, since it is little more than a front end for the
// internal work routine check_against_boundary.

// The domain (2nd Box) arguments to dir and box should be cell-based
// regardless of the type of the region.

class amr_boundary_class 
{
public:
    amr_boundary_class() { }
    virtual ~amr_boundary_class() { }

    int dir(const Box&, const Box&) const;
    void boundary_mesh(BoxArray& exterior_mesh, int*& grid_ref, const BoxArray& interior_mesh, const Box& domain) const;

    virtual Box box(const Box&, const Box&, int) const = 0;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, int, const Box&) const = 0;
    virtual void sync_borders(MultiFab&, const level_interface&) const = 0;
    virtual void fill_borders(MultiFab& r, const level_interface& lev_interface, int w) const = 0;
    virtual void duplicate(List<Box>&, const Box&) const = 0;
    virtual bool singular() const  = 0;
protected:
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const = 0;
    virtual void check_against_boundary(BoxList& bl, List<int>& il, const Box& b, int ib, const Box& d, int dim1) const = 0;
};

inline int amr_boundary_class::dir(const Box& region, const Box& domain) const
{
    for ( int i = 0; i < BL_SPACEDIM; ++i )
    {
	if ( region.bigEnd(i) < domain.smallEnd(i) )
	    return -(i+1);
	if ( region.smallEnd(i) > domain.bigEnd(i) )
	    return +(i+1);
    }
    for ( int i = 0; i < BL_SPACEDIM; ++i )
    {
	if ( region.bigEnd(i) == domain.smallEnd(i) )
	    return -(i+1);
	if ( region.smallEnd(i) == domain.bigEnd(i) )
	    return +(i+1);
    }
    BoxLib::Error("amr_boundary_class::dir---boundary box not outside domain.");
    return 0;
}

// This is the base class for physical boundary conditions of the
// types used in the AMR code.  Members return appropriate
// amr_boundary objects for various flow variables.  The defaults
// are all error_boundary, so this class should only be used as
// a base for derivation.

class amr_fluid_boundary_class 
{
public:
    amr_fluid_boundary_class();
    virtual ~amr_fluid_boundary_class() { }
    
    virtual const amr_boundary_class* velocity(int i) const 
    {
	return v[i];
    }
    virtual const amr_boundary_class* scalar() const 
    {
	return s;
    }
    virtual const amr_boundary_class* pressure() const 
    {
	return p;
    }
    virtual const amr_boundary_class* terrain_sigma() const 
    {
	return ts;
    }
protected:
    amr_boundary_class *v[BL_SPACEDIM];
    amr_boundary_class *s;
    amr_boundary_class *p;
    amr_boundary_class *ts;
};

// This is the object returned by the inviscid_fluid_boundary_class class.
// It chooses among periodic, reflection and negation boundary behavior
// based on the boundary types accessed through the pointer and the
// flow direction flowdim.
// flowdim = -1 indicates an advected scalar
// flowdim = -2 indicates a pressure
// flowdim = -3 indicates a streamfunction
// flowdim = -4 indicates the terrain sigma (all components)

class inviscid_fluid_boundary_class;

class mixed_boundary_class : public amr_boundary_class 
{
public:
    mixed_boundary_class(inviscid_fluid_boundary_class *Ptr, int idim)  : ptr(Ptr), flowdim(idim) {}
    
    virtual Box box(const Box&, const Box&, int) const;
    // virtual void fill(FArrayBox&, const Box&, const MultiFab&, int, const Box&) const;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, int, const Box&) const;
    virtual void sync_borders(MultiFab&, const level_interface&) const;
    virtual void fill_borders(MultiFab& r, const level_interface& lev_interface, int w) const;
    virtual void duplicate(List<Box>& bl, const Box& domain) const;
    virtual bool singular() const;
protected:
    const inviscid_fluid_boundary_class* ptr;
    const int flowdim;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const;
    virtual void check_against_boundary(BoxList& bl, List<int>& il, const Box& b, int ib, const Box& d, int dim1) const;
};

// This is a simple boundary condition with slip walls.
// Takes an array of RegType specifiers, one for each side of
// the domain.  Currently only periodic and refWall type
// boundaries are supported.

class inviscid_fluid_boundary_class : public amr_fluid_boundary_class 
{
    friend mixed_boundary_class;
public:
    inviscid_fluid_boundary_class(RegType Bc[BL_SPACEDIM][2]);
    virtual ~inviscid_fluid_boundary_class();
protected:
    RegType bc[BL_SPACEDIM][2];
};

#endif
