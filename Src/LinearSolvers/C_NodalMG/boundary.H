//BL_COPYRIGHT_NOTICE

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "interface.H"
#include "RegType.H"

//
// This is the base class for all amr_boundaries.  The default function dir
// should always be acceptable.  The default box and check_against_boundary
// functions are appropriate for wall-type boundaries.  The default
// boundary_mesh function is appropriate for all boundary types I've
// created so far, since it is little more than a front end for the
// internal work routine check_against_boundary.
//
//
// The domain (2nd Box) arguments to dir and box should be cell-based
// regardless of the type of the region.
//

class amr_boundary 
{
public:
    amr_boundary () {}

    virtual ~amr_boundary ();

    int dir (const Box&,
             const Box&) const;

    void boundary_mesh (BoxArray&       exterior_mesh,
                        int*&           grid_ref,
                        const BoxArray& interior_mesh,
                        const Box&      domain) const;

    virtual Box box (const Box&,
                     const Box&,
                     int) const = 0;

    virtual void fill (FArrayBox&,
                       const Box&,
                       const FArrayBox&,
                       const Box&) const = 0;

    virtual void sync_borders (MultiFab&,
                               const level_interface&) const = 0;

    virtual void fill_borders (MultiFab& r,
                               const level_interface& lev_interface,
                               int w) const = 0;

    virtual void duplicate (List<Box>&,
                            const Box&) const = 0;

    virtual bool singular () const  = 0;

    virtual Box image (const Box& region,
                       const Box& srcbox,
                       const Box& domain) const = 0;

protected:

    virtual void fill (FArrayBox&,
                       const Box&,
                       const FArrayBox&,
                       const Box&,
                       const Box&,
                       int) const = 0;

    virtual void check_against_boundary (BoxList&   bl,
                                         List<int>& il,
                                         const Box& b,
                                         int        ib,
                                         const Box& d,
                                         int        dim1) const = 0;
};

//
// This is the base class for physical boundary conditions of the
// types used in the AMR code.  Members return appropriate
// amr_boundary objects for various flow variables.  The defaults
// are all error_boundary, so this class should only be used as
// a base for derivation.
//

class amr_fluid_boundary 
{
public:

    amr_fluid_boundary ();

    virtual ~amr_fluid_boundary ();
    
    virtual const amr_boundary* velocity (int i) const;

    virtual const amr_boundary* scalar() const;

    virtual const amr_boundary* pressure() const;

    virtual const amr_boundary* terrain_sigma() const;

protected:
    amr_boundary *v[BL_SPACEDIM];
    amr_boundary *s;
    amr_boundary *p;
    amr_boundary *ts;
};

//
// This is the object returned by the inviscid_fluid_boundary class.
// It chooses among periodic, reflection and negation boundary behavior
// based on the boundary types accessed through the pointer and the
// flow direction flowdim.
// flowdim = -1 indicates an advected scalar
// flowdim = -2 indicates a pressure
// flowdim = -3 indicates a streamfunction
// flowdim = -4 indicates the terrain sigma (all components)
//

class inviscid_fluid_boundary;

class mixed_boundary
    : public amr_boundary 
{
public:

    mixed_boundary (inviscid_fluid_boundary* Ptr,
		    int                      idim)  
	:
        ptr(Ptr),
        flowdim(idim) 
    {}
    
    virtual Box box (const Box&,
                     const Box&,
                     int) const;

    virtual void fill (FArrayBox&,
                       const Box&,
                       const FArrayBox&,
                       const Box&) const;

    virtual void sync_borders (MultiFab&,
                               const level_interface&) const;

    virtual void fill_borders (MultiFab&              r,
                               const level_interface& lev_interface,
                               int                    w) const;

    virtual void duplicate (List<Box>& bl,
                            const Box& domain) const;

    virtual bool singular () const;

    virtual Box image (const Box& region,
                       const Box& srcbox,
                       const Box& domain) const;

protected:

    virtual void fill (FArrayBox&,
                       const Box&,
                       const FArrayBox&,
                       const Box&,
                       const Box&,
                       int) const;

    virtual void check_against_boundary (BoxList&   bl,
                                         List<int>& il,
                                         const Box& b,
                                         int        ib,
                                         const Box& d,
                                         int        dim1) const;
    //
    // The data.
    //
    const inviscid_fluid_boundary* ptr;
    const int flowdim;
};

// This is a simple boundary condition with slip walls.
// Takes an array of RegType specifiers, one for each side of
// the domain.  Currently only periodic and refWall type
// boundaries are supported.

class inviscid_fluid_boundary
    : public amr_fluid_boundary 
{
    friend class mixed_boundary;
public:
    explicit inviscid_fluid_boundary(RegType Bc[BL_SPACEDIM][2]);
    virtual ~inviscid_fluid_boundary();
    RegType getLoBC(int idim) const;
    RegType getHiBC(int idim) const;
protected:
    RegType bc[BL_SPACEDIM][2];
};

#endif /*_BOUNDARY_H_*/
