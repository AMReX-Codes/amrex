//BL_COPYRIGHT_NOTICE

#ifndef _AMR_MULTI_H_
#define _AMR_MULTI_H_

#include "fill_patch.H"

class amr_multigrid
{
public:

    amr_multigrid (const Array<BoxArray>&    Mesh,
                   const Array<IntVect>&     Gen_ratio,
                   int                       Lev_min_min,
                   int                       Lev_min_max,
                   int                       Lev_max_max,
                   const amr_boundary* Boundary,
                   int                       Pcode) ;
    virtual ~amr_multigrid ();

#if BL_SPACEDIM==2
    enum CoordSys { cartesian = 0, rz = 1 };
#else
    enum CoordSys { cartesian = 0 };
#endif
    CoordSys setCoordSys (CoordSys c_sys_);
    CoordSys getCoordSys () const;

    static void mesh_read (Array<BoxArray>& m,
                           Array<IntVect>&  r,
                           Array<Box>&      d,
                           istream&         is);

    static void mesh_write (const Array<BoxArray>& m,
                            const Array<Box>&      d,
                            ostream&               os);

    static void mesh_write (const Array<BoxArray>& m,
                            const Array<IntVect>&  r,
                            Box                    fd,
                            ostream&               os);

    void solve (Real reltol,
                Real abstol,
                int  i1,
                int  i2);

    virtual void mg_interpolate_level (int lto,
                                       int lfrom) = 0;

    virtual void mg_restrict_level (int lto,
                                    int lfrom) = 0;

    virtual void level_residual (MultiFab& r,
                                 MultiFab& s,
                                 MultiFab& d,
                                 int       mglev,
                                 bool      iclear,
                                 int       for_sync_resid) = 0;

    virtual void interface_residual (int mglev,
                                     int mgc) = 0;

    virtual void relax (int  mglev,
                        int  i1,
                        bool is_zero) = 0;

    virtual void cgsolve (int mglev) = 0;

protected:

    Real ml_cycle (int lev,
                   int mglev,
                   int i1,
                   int i2,
                   Real tol,
                   Real res_norm_fine);

    void mg_cycle (int mglev,
                   int i1,
                   int i2);

    Real ml_residual (int mglev,
                      int lev);

    void build_index ();

    int get_amr_level (int mglev) const ;
    //
    // All of these have built-in defaults except can_coarsen().
    //
    void build_mesh (const Box& fdomain);

    int build_down (const BoxArray& l_mesh,
                    const Box&      l_domain,
                    int             flev,
                    IntVect         rat,
                    int             nlev);

    void make_coarser_level (BoxArray& mesh,
                             Box&      domain,
                             int&      flev,
                             IntVect&  rat);
    //
    // Before using the multigrid, a derived class must call build_mesh()
    // followed by alloc(), or their equivalents.  The intended use is
    // that build_mesh() be called by the derived class constructor, while
    // alloc(), which actually creates large data objects, and clear(),
    // which deletes them, be called when the solver is actually being used.
    // These are made non-virtual so that a derived class can define a more
    // complicated lev_interface.
    //
    void alloc_amr_multi (PArray<MultiFab>& Dest,
			  PArray<MultiFab>& Source,
			  PArray<MultiFab>& Coarse_source,
			  int               Lev_min,
			  int               Lev_max);

    void clear_amr_multi ();

    virtual void sync_interfaces () = 0;

    virtual bool can_coarsen (const BoxArray& mesh,
                              const Box&      domain) const = 0;
    //
    // The data.
    //
    const int                 lev_min_min;
    const int                 lev_min_max;
    const int                 lev_max_max;
    int                       lev_min;
    int                       lev_max;
    int                       mglev_max;
    Array<IntVect>            gen_ratio;          // assoc with ml_mesh
    Array<int>                ml_index;           // assoc with ml_mesh
    Array<Box>                mg_domain;          // assoc with mg_mesh
    Array<BoxArray>           ml_mesh;
    Array<BoxArray>           mg_mesh;
    Array< Array<Box> >       mg_domain_array;
    Array< Array<BoxArray> >  mg_mesh_array;
    const amr_boundary* mg_boundary;
    level_interface*          lev_interface;
    PArray<MultiFab>          source;             // on active levels of ml_mesh
    PArray<MultiFab>          dest;
    PArray<MultiFab>          resid;              // on mg_mesh
    PArray<MultiFab>          corr;
    PArray<MultiFab>          work;
    PArray<MultiFab>          save;               // on middle levels of ml_mesh
    PArray<MultiFab>          coarse_source;      // on ml_mesh, if it exists


    CoordSys                  c_sys;

    int                       pcode;
    int                       integrate;
private:
    level_interface**         interface_array;

};

#endif /*_AMR_MULTI_H_*/
