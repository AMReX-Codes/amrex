#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "interface.H"
#include "RegType.H"

#ifdef HG_USE_CACHE
class copy_cache;
#endif

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
protected:
    virtual void check_against_boundary(BoxList& bl, List<int>& il,
	const Box& b, int ib,
	const Box& d, int dim1) const;
public:
    amr_boundary_class() { }
    virtual ~amr_boundary_class() { }
    int dir(const Box&, const Box&) const;
    virtual bool defined() const { return false; }
    virtual Box box(const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const
    { BoxLib::Error("amr_boundary_class::fill---boundary fill function not defined."); }
    virtual void fill(FArrayBox&, const Box&, const MultiFab&, int, const Box&) const
    { BoxLib::Error("amr_boundary_class::fill---boundary fill function not defined."); }
    virtual void sync_borders(MultiFab&, const level_interface&) const { }
    virtual void fill_borders(MultiFab& r, const level_interface& lev_interface, int w = -1) const
    { BoxLib::Error("amr_boundary_class::fill---fill_borders function not defined."); }
#ifdef HG_USE_CACHE
    virtual void set_sync_cache(copy_cache*, int, int&, MultiFab&, const level_interface&) const { }
    virtual void set_border_cache(copy_cache*, int, int&, MultiFab&, const level_interface&, int) const 
    { BoxLib::Error("amr_boundary_class::set_border_cache---boundary condition not defined."); }
#endif
    virtual void boundary_mesh(BoxArray& exterior_mesh, int *&grid_ref, const BoxArray& interior_mesh, const Box& domain) const;
    virtual void duplicate(List<Box>&, const Box&) const { }
    // virtual int fill_patch_special(FArrayBox&, const Box&, MultiFab&, const Box&, int) const { return 0; }
    virtual bool singular() const 
    {
	BoxLib::Error("amr_boundary_class::singular---boundary condition not defined.");
	return false;
    }
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
	    return (i+1);
    }
    BoxLib::Error("amr_boundary_class::dir---boundary box not outside domain.");
    return 0;
}

class reflection_boundary_class : public amr_boundary_class 
{
public:
    reflection_boundary_class() { }
    virtual bool defined() const { return true; }
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&, const MultiFab&, int, const Box&) const
    { BoxLib::Error("reflection_boundary_class::fill---lev_interface version not defined"); }
    virtual bool singular() const { return true; }
};

class negation_boundary_class : public amr_boundary_class 
{
public:
    negation_boundary_class() { }
    virtual bool defined() const { return true; }
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&, const MultiFab&, int, const Box&) const
    { BoxLib::Error("negation_boundary_class::fill---lev_interface version not defined"); }
    virtual bool singular() const { return false; }
};

class periodic_boundary_class : public amr_boundary_class 
{
protected:
    void check_against_boundary(BoxList& bl, List<int>& il,
	const Box& b, int ib,
	const Box& d, int dim1) const;
public:
    periodic_boundary_class() { }
    virtual bool defined() const { return true; }
    virtual Box box(const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&, const FArrayBox&, const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&, const MultiFab&, int, const Box&) const
    { BoxLib::Error("periodic_boundary_class::fill---lev_interface version not defined"); }
    virtual void duplicate(List<Box>& bl, const Box& domain) const;
    virtual bool singular() const { return true; }
};

const amr_boundary_class error_boundary;

// This is the base class for physical boundary conditions of the
// types used in the AMR code.  Members return appropriate
// amr_boundary objects for various flow variables.  The defaults
// are all error_boundary, so this class should only be used as
// a base for derivation.

class amr_fluid_boundary_class 
{
protected:
    const amr_boundary_class *v[BL_SPACEDIM];
    const amr_boundary_class *s;
    const amr_boundary_class *p;
    const amr_boundary_class *ts;
public:
    amr_fluid_boundary_class();
    virtual ~amr_fluid_boundary_class() { }
    
    virtual const amr_boundary_class& velocity(int i) const 
    {
	return *v[i];
    }
    virtual const amr_boundary_class& scalar() const 
    {
	return *s;
    }
    virtual const amr_boundary_class& pressure() const 
    {
	return *p;
    }
    virtual const amr_boundary_class& terrain_sigma() const 
    {
	return *ts;
    }
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
protected:
    inviscid_fluid_boundary_class *ptr;
    int flowdim;
    void check_against_boundary(BoxList& bl, List<int>& il,
	const Box& b, int ib,
	const Box& d, int dim1) const;
public:
    mixed_boundary_class(inviscid_fluid_boundary_class *Ptr, int idim) 
    {
	ptr = Ptr; flowdim = idim;
    }
    
    virtual bool defined() const { return true; }
    virtual Box box(const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&,
	const FArrayBox&, const Box&, const Box&, int) const;
    virtual void fill(FArrayBox&, const Box&,
	const MultiFab&, int, const Box&) const;
    virtual void sync_borders(MultiFab&, const level_interface&) const;
    virtual void fill_borders(MultiFab& r, const level_interface& lev_interface, int w = -1) const;
#ifdef HG_USE_CACHE
    virtual void set_sync_cache(copy_cache* cache, int nsets, int& iset,
		      MultiFab& r,
		      const level_interface& lev_interface) const;
    virtual void set_border_cache(copy_cache* cache, int nsets, int& iset,
	MultiFab& r,
	const level_interface& lev_interface,
	int w) const;
#endif
    virtual void duplicate(List<Box>& bl, const Box& domain) const;
    // int fill_patch_special(FArrayBox&, const Box&, MultiFab&, const Box&, int flags) const;
    virtual bool singular() const;
};

// This is a simple boundary condition with slip walls.
// Takes an array of RegType specifiers, one for each side of
// the domain.  Currently only periodic and refWall type
// boundaries are supported.

class inviscid_fluid_boundary_class : public amr_fluid_boundary_class 
{
    friend mixed_boundary_class;
protected:
    RegType bc[BL_SPACEDIM][2];
public:
    inviscid_fluid_boundary_class(RegType Bc[BL_SPACEDIM][2]);
    ~inviscid_fluid_boundary_class();
};

#endif
