//BL_COPYRIGHT_NOTICE

#include "hg_projector.H"

#include <RunStats.H>

#if defined( BL_FORT_USE_UNDERSCORE )
#define   FORT_HGDIV		hgdiv_
#define   FORT_HGDIV_TERRAIN    hgdiv_terrain_
#define   FORT_HGFDIV		hgfdiv_
#define   FORT_HGFDIV_TERRAIN   hgfdiv_terrain_
#define   FORT_HGEDIV		hgediv_
#define   FORT_HGEDIV_TERRAIN   hgediv_terrain_
#define   FORT_HGCDIV		hgcdiv_
#define   FORT_HGCDIV_TERRAIN   hgcdiv_terrain_
#define   FORT_HGGRAD		hggrad_
#define   FORT_HGGRAD_TERRAIN   hggrad_terrain_
#define   FORT_HGAVG		hgavg_
#define   FORT_HGFAVG		hgfavg_
#define   FORT_HGEAVG		hgeavg_
#define   FORT_HGCAVG		hgcavg_
#elif defined( BL_FORT_USE_UPPERCASE )
#define   FORT_HGDIV		HGDIV
#define   FORT_HGDIV_TERRAIN    HGDIV_TERRAIN
#define   FORT_HGFDIV		HGFDIV
#define   FORT_HGFDIV_TERRAIN   HGFDIV_TERRAIN
#define   FORT_HGEDIV		HGEDIV
#define   FORT_HGEDIV_TERRAIN   HGEDIV_TERRAIN
#define   FORT_HGCDIV		HGCDIV
#define   FORT_HGCDIV_TERRAIN   HGCDIV_TERRAIN
#define   FORT_HGGRAD		HGGRAD
#define   FORT_HGGRAD_TERRAIN   HGGRAD_TERRAIN
#define   FORT_HGAVG		HGAVG
#define   FORT_HGFAVG		HGFAVG
#define   FORT_HGEAVG		HGEAVG
#define   FORT_HGCAVG		HGCAVG
#elif defined( BL_FORT_USE_LOWERCASE )
#define   FORT_HGDIV		hgdiv
#define   FORT_HGDIV_TERRAIN    hgdiv_terrain
#define   FORT_HGFDIV		hgfdiv
#define   FORT_HGFDIV_TERRAIN   hgfdiv_terrain
#define   FORT_HGEDIV		hgediv
#define   FORT_HGEDIV_TERRAIN   hgediv_terrain
#define   FORT_HGCDIV		hgcdiv
#define   FORT_HGCDIV_TERRAIN   hgcdiv_terrain
#define   FORT_HGGRAD		hggrad
#define   FORT_HGGRAD_TERRAIN   hggrad_terrain
#define   FORT_HGAVG		hgavg
#define   FORT_HGFAVG		hgfavg
#define   FORT_HGEAVG		hgeavg
#define   FORT_HGCAVG		hgcavg
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C" 
{
    
#if (BL_SPACEDIM == 1)
#error not relevant
#elif (BL_SPACEDIM == 2 || BL_SPACEDIM == 3)
    void FORT_HGGRAD_TERRAIN(RealPS, intS, const Real*, intS, intS, CRealPS);
    void FORT_HGGRAD        (RealPS, intS, const Real*, intS, intS, CRealPS, const int*);

    void FORT_HGDIV_TERRAIN (Real*,  intS, CRealPS, intS, intS, CRealPS);
    void FORT_HGFDIV_TERRAIN(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGEDIV_TERRAIN(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGCDIV_TERRAIN(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
    

    void FORT_HGDIV         (Real*, intS, CRealPS, intS, intS, CRealPS, const int*, const int*);

    void FORT_HGEDIV        (Real*, intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
#if (BL_SPACEDIM == 2)
    void FORT_HGAVG         (Real*, intS, const Real*, intS, intS, const Real*, const int*, const int*);
    void FORT_HGFAVG        (Real*, intS, const Real*, intS, const Real*, intS, intS, intRS, const int*, const int*, const Real*, const int*, const int*);
    void FORT_HGCAVG        (Real*, intS, const Real*, intS, const Real*, intS, intS, intRS, const int*, const int*, const Real*, const int*, const int*);
    void FORT_HGFDIV        (Real*, intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
    void FORT_HGCDIV        (Real*, intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
#elif (BL_SPACEDIM == 3)
    void FORT_HGFDIV        (Real*, intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGCDIV        (Real*, intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGAVG         (Real*, intS, const Real*, intS, intS);
    void FORT_HGFAVG        (Real*, intS, const Real*, intS, const Real*, intS, intS, intRS, const int*, const int*);
    void FORT_HGEAVG        (Real*, intS, const Real*, intS, const Real*, intS, intS, intRS, const int*, const int*);
    void FORT_HGCAVG        (Real*, intS, const Real*, intS, const Real*, intS, intS, intRS, const int*, const int*);
#endif
#endif
}

extern "C"
{
#if BL_SPACEDIM==2
  typedef void (*FECFUNC)(Real*, intS, const Real*, intS, const Real*, intS,
			  intS, intRS, const int*, const int*, const Real*, const int*, const int*);
#else
  typedef void (*FECFUNC)(Real*, intS, const Real*, intS, const Real*, intS,
			  intS, intRS, const int*, const int*);
#endif
}

class task_fecavg : public task_fec_base
{
public:

    task_fecavg (FECFUNC         f_,
                 task_list&      tl_,
                 MultiFab&       s_,
                 const MultiFab& S_,
                 int             igrid_,
                 task_fab*       tf_,
                 const Box&      creg_,
                 const IntVect&  rat_,
                 int             idim_,
                 int             idir_
#if BL_SPACEDIM == 2
                 , Real hx_, int isRZ_, int imax_
#endif
	);

    virtual bool ready ();

private:

    void doit ();

    FECFUNC         f;
    const MultiFab& S;
    const Box       creg;
    const IntVect   rat;
    const int       idim;
    const int       idir;
#if BL_SPACEDIM == 2
    Real            hx;
    int             isRZ;
    int             imax;
#endif
};

task_fecavg::task_fecavg (FECFUNC         f_,
                          task_list&      tl_,
                          MultiFab&       s_,
                          const MultiFab& S_,
                          int             igrid_,
                          task_fab*       tf_,
                          const Box&      creg_,
                          const IntVect&  rat_,
                          int             idim_,
                          int             idir_
#if BL_SPACEDIM == 2
                          , Real hx_, int isRZ_, int imax_
#endif
    )
    :
    task_fec_base(tl_,s_,igrid_),
    f(f_),
    S(S_),
    creg(creg_),
    rat(rat_),
    idim(idim_),
    idir(idir_)
#if BL_SPACEDIM == 2
    , hx(hx_), isRZ(isRZ_), imax(imax_)
#endif
{
    push_back(tf_);

    if (is_local_target() && dependencies.empty()) doit();
}

bool
task_fecavg::ready ()
{
    BL_ASSERT(!done);

    if (is_local_target()) doit();

    return true;
}

void
task_fecavg::doit ()
{
    BL_ASSERT(!done);
    BL_ASSERT(is_local_target());
    BL_ASSERT(dependencies.empty());

    done = true;

    const int        igrid    = grid_number();
    FArrayBox&       sfab     = target_fab();
    const Box&       sfab_box = sfab.box();
    const FArrayBox& cfab     = task_fab_result(0);
    const Box&       cfab_box = cfab.box();
    const FArrayBox& Sfab     = S[igrid];
    const Box&       Sfab_box = Sfab.box();

    (*f)(sfab.dataPtr(), DIMLIST(sfab_box), cfab.dataPtr(), DIMLIST(cfab_box), Sfab.dataPtr(), DIMLIST(Sfab_box), DIMLIST(creg), D_DECL(rat[0], rat[1], rat[2]), &idim, &idir
#if BL_SPACEDIM == 2
         , &hx, &isRZ, &imax
#endif
        );
}

class task_fecavg_2 : public task_fec_base
{
public:

    task_fecavg_2 (FECFUNC           f_,
                   task_list&        tl_,
                   MultiFab&         s_,
                   int               igrid_,
                   task_fab*         Sfp_,
                   task_fab*         Scp_,
                   const Box&        creg_,
                   const IntVect&    rat_,
                   const Array<int>& ga_,
                   const IntVect&    t_
#if BL_SPACEDIM==2
                   , Real hx_, int isRZ_, int imax_
#endif
    );

    virtual bool ready ();

private:

    void doit ();

    FECFUNC          f;
    const Box        creg;
    const IntVect    rat;
    const IntVect    t;
    const Array<int> ga;
#if BL_SPACEDIM == 2
    Real             hx;
    int              isRZ;
    int              imax;
#endif
};

task_fecavg_2::task_fecavg_2 (FECFUNC           f_,
                              task_list&        tl_,
                              MultiFab&         s_,
                              int               igrid_,
                              task_fab*         Sfp_,
                              task_fab*         Scp_,
                              const Box&        creg_,
                              const IntVect&    rat_,
                              const Array<int>& ga_,
                              const IntVect&    t_
#if BL_SPACEDIM==2
                              , Real hx_, int isRZ_, int imax_
#endif
    )
    :
    task_fec_base(tl_,s_,igrid_),
    f(f_),
    creg(creg_),
    rat(rat_),
    t(t_),
    ga(ga_) 
#if BL_SPACEDIM == 2
    , hx(hx_), isRZ(isRZ_), imax(imax_)	
#endif
{
    push_back(Sfp_);
    push_back(Scp_);

    if (is_local_target() && dependencies.empty()) doit();
}

bool
task_fecavg_2::ready ()
{
    BL_ASSERT(!done);

    if (is_local_target()) doit();

    return true;
}

void
task_fecavg_2::doit ()
{
    BL_ASSERT(!done);
    BL_ASSERT(is_local_target());
    BL_ASSERT(dependencies.empty());

    done = true;

    const int        igrid      = grid_number();
    FArrayBox&       sfab       = target_fab();
    const Box&       sfab_box   = sfab.box();
    const FArrayBox& Sf_fab     = task_fab_result(0);
    const Box&       Sf_fab_box =  Sf_fab.box();
    const FArrayBox& Sc_fab     = task_fab_result(1);
    const Box&       Sc_fab_box = Sc_fab.box();

    (*f)(sfab.dataPtr(), DIMLIST(sfab_box), Sc_fab.dataPtr(), DIMLIST(Sc_fab_box), Sf_fab.dataPtr(), DIMLIST(Sf_fab_box), DIMLIST(creg), D_DECL( rat[0], rat[1], rat[2]), ga.dataPtr(), t.getVect()
#if BL_SPACEDIM == 2
         , &hx, &isRZ, &imax
#endif
        );
}

extern "C"
{
    typedef void (*FECDIV)(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int* 
#if BL_SPACEDIM == 2
			 , const int*, const int *
#endif 
			 );
}

class task_fecdiv : public task_fec_base
{
public:
    task_fecdiv (FECDIV         f_,
                 task_list&     tl_,
                 MultiFab&      s_,
                 MultiFab*      upt_[],
                 int            igrid_,
                 task_fab*      ucp_[],
                 const Box&     creg_,
                 const Real*    h_,
                 const IntVect& rat_,
                 int            idim_,
                 int            idir_
#if BL_SPACEDIM == 2
                 ,int isRZ_, int imax_
#endif
	);

    virtual ~task_fecdiv ();
    virtual bool ready ();

private:

    void doit ();

    FECDIV        f;
    MultiFab*     upt[BL_SPACEDIM];
    const Box     creg;
    Real          h[BL_SPACEDIM];
    const IntVect rat;
    const int     idim;
    const int     idir;
#if BL_SPACEDIM  == 2
    const int     isRZ;
    const int     imax;
#endif
};

task_fecdiv::task_fecdiv (FECDIV         f_,
                          task_list&     tl_,
                          MultiFab&      s_,
                          MultiFab*      upt_[],
                          int            igrid_,
                          task_fab*      ucp_[],
                          const Box&     creg_,
                          const Real*    h_,
                          const IntVect& rat_,
                          int            idim_,
                          int            idir_
#if BL_SPACEDIM == 2
                          ,int isRZ_, int imax_
#endif
    )
    :
    task_fec_base(tl_,s_,igrid_),
    f(f_),
    creg(creg_),
    rat(rat_),
    idim(idim_),
    idir(idir_)
#if BL_SPACEDIM == 2
    , isRZ(isRZ_), imax(imax_)
#endif
{
    for (int i = 0; i  < BL_SPACEDIM; ++i)
    {
        upt[i] = upt_[i];
        push_back(ucp_[i]);
        h[i] = h_[i];
    }

    if (is_local_target() && dependencies.empty()) doit();
}

bool
task_fecdiv::ready ()
{
    BL_ASSERT(!done);

    if (is_local_target()) doit();

    return true;
}

void
task_fecdiv::doit ()
{
    BL_ASSERT(!done);
    BL_ASSERT(is_local_target());
    BL_ASSERT(dependencies.empty());

    done = true;

    const int   igrid           = grid_number();
    FArrayBox&  s               = target_fab();
    const Box&  s_box           = target_fab().box();
    const Real* up[BL_SPACEDIM] = { D_DECL( upt[0]->operator[](igrid).dataPtr(), upt[1]->operator[](igrid).dataPtr(), upt[2]->operator[](igrid).dataPtr() ) };
    const Box&  up_box          = upt[0]->operator[](igrid).box();
    const Real* uc[BL_SPACEDIM] = { D_DECL( task_fab_result(0).dataPtr(), task_fab_result(1).dataPtr(), task_fab_result(2).dataPtr() ) };
    const Box&  uc_box          = task_fab_result(0).box();

    (*f)(s.dataPtr(), DIMLIST(s_box), D_DECL( uc[0], uc[1], uc[2]), DIMLIST(uc_box), D_DECL(up[0], up[1], up[2]), DIMLIST(up_box), DIMLIST(creg), D_DECL(&h[0], &h[1], &h[2]), D_DECL(rat[0], rat[1], rat[2]), &idim, &idir
#if BL_SPACEDIM == 2
         , &isRZ, &imax
#endif
        );
}

task_fecdiv::~task_fecdiv ()
{
    HG_DEBUG_OUT("destroying a task_fecdiv " << this << endl);
}

extern "C"
{
  typedef void (*FECDIV_2)(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
}

class task_fecdiv_2 : public task_fec_base
{
public:
    task_fecdiv_2 (FECDIV_2           f_,
                   task_list&        tl_,
                   MultiFab&         s_,
                   int               igrid_,
                   task_fab*         ufp_[],
                   task_fab*         ucp_[],
                   const Box&        creg_,
                   const Real*       h_,
                   const IntVect&    rat_,
                   const Array<int>& ga_,
                   const IntVect&    t_ = IntVect());

    virtual bool ready ();
private:

    void doit ();

    FECDIV_2         f;
    const Box        creg;
    Real             h[BL_SPACEDIM];
    const IntVect    rat;
    const Array<int> ga;
    const IntVect    t;
};

task_fecdiv_2::task_fecdiv_2 (FECDIV_2            f_,
                              task_list&        tl_,
                              MultiFab&         s_,
                              int               igrid_,
                              task_fab*         ufp_[],
                              task_fab*         ucp_[],
                              const Box&        creg_,
                              const Real*       h_,
                              const IntVect&    rat_,
                              const Array<int>& ga_,
                              const IntVect&    t_)
    :
    task_fec_base(tl_,s_,igrid_),
    f(f_),
    creg(creg_),
    rat(rat_),
    ga(ga_),
    t(t_)
{
    for (int i = 0; i  < BL_SPACEDIM; ++i)
    {
        push_back(ucp_[i]);
    }
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        push_back(ufp_[i]);
        h[i]   = h_[i];
    }

    if (is_local_target() && dependencies.empty()) doit();
}

bool
task_fecdiv_2::ready ()
{
    BL_ASSERT(!done);

    if (is_local_target()) doit();

    return true;
}

void
task_fecdiv_2::doit ()
{
    BL_ASSERT(!done);
    BL_ASSERT(is_local_target());
    BL_ASSERT(dependencies.empty());

    done = true;

    const int   igrid           = grid_number();
    FArrayBox&  s               = target_fab();
    const Box&  s_box           = target_fab().box();
    const Real* uf[BL_SPACEDIM] = { D_DECL( task_fab_result(BL_SPACEDIM).dataPtr(), task_fab_result(BL_SPACEDIM+1).dataPtr(), task_fab_result(BL_SPACEDIM+2).dataPtr() ) };
    const Box&  uf_box          = task_fab_result(BL_SPACEDIM).box();
    const Real* uc[BL_SPACEDIM] = { D_DECL( task_fab_result(0).dataPtr(), task_fab_result(1).dataPtr(), task_fab_result(2).dataPtr() ) };
    const Box&  uc_box          = task_fab_result(0).box();

    (*f)(s.dataPtr(), DIMLIST(s_box), D_DECL( uc[0], uc[1], uc[2]), DIMLIST(uc_box), D_DECL(uf[0], uf[1], uf[2]), DIMLIST(uf_box), DIMLIST(creg), D_DECL(&h[0], &h[1], &h[2]), D_DECL(rat[0], rat[1], rat[2]), ga.dataPtr(), t.getVect());
}

extern "C"
{
  typedef void (*FECDIV_3)(Real*,  intS, CRealPS, intS, CRealPS, intS, intS, CRealPS, intRS, const int*, const int*);
}

class task_fecdiv_3 : public task_fec_base
{
public:
    task_fecdiv_3 (FECDIV_3          f_,
                   task_list&        tl_,
                   MultiFab&         s_,
                   int               igrid_,
                   task_fab*         ufp_[],
                   task_fab*         ucp_[],
                   const Box&        creg_,
                   const Real*       h_,
                   const IntVect&    rat_,
                   const Array<int>& ga_,
                   const int         isrz_);

    virtual bool ready ();
private:

    void doit ();

    FECDIV_3         f;
    const Box        creg;
    Real             h[BL_SPACEDIM];
    const IntVect    rat;
    const Array<int> ga;
    int              isrz;
};

task_fecdiv_3::task_fecdiv_3 (FECDIV_3          f_,
                              task_list&        tl_,
                              MultiFab&         s_,
                              int               igrid_,
                              task_fab*         ufp_[],
                              task_fab*         ucp_[],
                              const Box&        creg_,
                              const Real*       h_,
                              const IntVect&    rat_,
                              const Array<int>& ga_,
                              int               isrz_)
    :
    task_fec_base(tl_,s_,igrid_),
    f(f_),
    creg(creg_),
    rat(rat_),
    ga(ga_),
    isrz(isrz_)
{
    for (int i = 0; i  < BL_SPACEDIM; ++i)
    {
        push_back(ucp_[i]);
    }
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        push_back(ufp_[i]);
        h[i]   = h_[i];
    }

    if (is_local_target() && dependencies.empty()) doit();
}

bool
task_fecdiv_3::ready ()
{
    BL_ASSERT(!done);

    if (is_local_target()) doit();

    return true;
}

void
task_fecdiv_3::doit ()
{
    BL_ASSERT(!done);
    BL_ASSERT(is_local_target());
    BL_ASSERT(dependencies.empty());

    done = true;

    const int   igrid           = grid_number();
    FArrayBox&  s               = target_fab();
    const Box&  s_box           = target_fab().box();
    const Real* uf[BL_SPACEDIM] = { D_DECL( task_fab_result(BL_SPACEDIM).dataPtr(), task_fab_result(BL_SPACEDIM+1).dataPtr(), task_fab_result(BL_SPACEDIM+2).dataPtr() ) };
    const Box&  uf_box          = task_fab_result(BL_SPACEDIM).box();
    const Real* uc[BL_SPACEDIM] = { D_DECL( task_fab_result(0).dataPtr(), task_fab_result(1).dataPtr(), task_fab_result(2).dataPtr() ) };
    const Box&  uc_box          = task_fab_result(0).box();

    (*f)(s.dataPtr(), DIMLIST(s_box), D_DECL( uc[0], uc[1], uc[2]), DIMLIST(uc_box), D_DECL(uf[0], uf[1], uf[2]), DIMLIST(uf_box), DIMLIST(creg), D_DECL(&h[0], &h[1], &h[2]), D_DECL(rat[0], rat[1], rat[2]), ga.dataPtr(), &isrz);
}

PArray<MultiFab> null_amr_real;

void
holy_grail_amr_projector::project (PArray<MultiFab>* u,
                                   PArray<MultiFab>& p,
                                   PArray<MultiFab>& Coarse_source,
                                   PArray<MultiFab>& Sigma,
                                   Real              H[],
                                   Real              tol,
                                   int               Lev_min,
                                   int               Lev_max,
                                   Real              scale)
{
    static RunStats hg_proj("hg_project");

    hg_proj.start();
    if (Lev_min < 0)
	Lev_min = lev_min_max;
    if (Lev_max < 0)
	Lev_max = Lev_min;
    
    BL_ASSERT(Sigma.length() > 0);
    
    BL_ASSERT(u[      0      ][Lev_min].nGrow() == 1);
    BL_ASSERT(u[      1      ][Lev_min].nGrow() == 1);
    BL_ASSERT(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);
    
    alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
    right_hand_side(u, null_amr_real);
    if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) 
    {
	sparse_node_source_adjustment(coarse_source);
    }
    
    solve(tol, scale, 2, 2);
    form_solution_vector(u, Sigma);
    clear();
    hg_proj.end();
}

void
holy_grail_amr_projector::sync_project (PArray<MultiFab>* u,
                                        PArray<MultiFab>& p,
                                        PArray<MultiFab>& Coarse_source,
                                        PArray<MultiFab>& Sigma,
                                        Real              H[],
                                        Real              tol,
                                        int               Lev_min,
                                        int               Lev_max,
                                        Real              scale)
{
    if (Lev_min < 0)
	Lev_min = lev_min_max;
    if (Lev_max < 0)
	Lev_max = Lev_min;
    
    BL_ASSERT(Sigma.length() > 0);
    
    BL_ASSERT(u[      0      ][Lev_min].nGrow() == 1);
    BL_ASSERT(u[      1      ][Lev_min].nGrow() == 1);
    BL_ASSERT(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);
    
    alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
    sync_right_hand_side(u);
    if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) 
    {
	sparse_node_source_adjustment(coarse_source);
    }
    
    solve(tol, scale, 2, 2);
    form_solution_vector(u, Sigma);
    clear();
}

void
holy_grail_amr_projector::manual_project (PArray<MultiFab>* u,
                                          PArray<MultiFab>& p,
                                          PArray<MultiFab>& rhs,
                                          PArray<MultiFab>& Coarse_source,
                                          PArray<MultiFab>& Sigma,
                                          bool              use_u,
                                          Real              H[],
                                          Real              tol,
                                          int               Lev_min,
                                          int               Lev_max,
                                          Real              scale)
{
    static RunStats hg_man_proj("hg_manual_project");

    hg_man_proj.start();
    if (Lev_min < 0)
	Lev_min = lev_min_max;
    if (Lev_max < 0)
	Lev_max = Lev_min;
    
    BL_ASSERT(Sigma.length() > 0);
    
    if (use_u) 
    {
	BL_ASSERT(u[      0      ][Lev_min].nGrow() == 1);
	BL_ASSERT(u[      1      ][Lev_min].nGrow() == 1);
	BL_ASSERT(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);
	
	alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
	if (rhs.length() > 0) 
	{
	    if (type(rhs[Lev_min]) == IntVect::TheNodeVector()) 
	    {
		right_hand_side(u, null_amr_real);
		for (int lev = Lev_min; lev <= Lev_max; lev++)
		    source[lev].plus(rhs[lev], 0, 1, 0);
		if (singular && make_sparse_node_source_solvable) 
		{
                    //
		    // Note:  You don't want to do this if rhs is not sparse!
                    //
		    sparse_node_source_adjustment(rhs);
		}
	    }
	    else 
	    {
		BL_ASSERT(rhs[Lev_min].nGrow() == 1);
		right_hand_side(u, rhs);
	    }
	}
	else 
	{
	    right_hand_side(u, null_amr_real);
	}
    }
    else 
    {
	BL_ASSERT(rhs.length() > 0);
	BL_ASSERT(rhs[Lev_min].nGrow() == 1);
	
	if (type(rhs[Lev_min]) == IntVect::TheNodeVector()) 
	{
	    alloc(p, rhs, Coarse_source, Sigma, H, Lev_min, Lev_max);
	    if (singular && make_sparse_node_source_solvable) 
	    {
                //
		// Note:  You don't want to do this if rhs is not sparse!
                //
		sparse_node_source_adjustment(rhs);
	    }
	}
	else 
	{
	    alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
            //
	    // Source is set to 0 at this point.
            //
	    right_hand_side(0, rhs);
	}
    }
    if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) 
    {
	sparse_node_source_adjustment(coarse_source);
    }
    
    solve(tol, scale, 2, 2);
    form_solution_vector(u, Sigma);
    clear();
    hg_man_proj.end();
}

void
holy_grail_amr_projector::sparse_node_source_adjustment (PArray<MultiFab>& sparse_source)
{
    //
    // This routine takes advantage of the sparse structure of
    // the sync source to avoid costly restriction operations.  It
    // is necessary to use the inner_product routine, which weights
    // boundary nodes, since the coarse-fine lev_interface can touch
    // the boundary.  Otherwise a call to sum would suffice.
    //
    // Note that the correction is applied to source, not to
    // sparse_source, since the sparse structure of the latter
    // may need to be preserved.
    //
    BL_ASSERT(singular);
    BL_ASSERT(make_sparse_node_source_solvable);
    
    Real adjust = 0.0;
    for (int lev = lev_max; lev >= lev_min; lev--) 
    {
	if (sparse_source.defined(lev)) 
	{
	    int mglev = ml_index[lev];
	    corr[mglev].setVal(1.0);
	    adjust += inner_product(sparse_source[lev], corr[mglev]);
	}
	if (lev > lev_min && adjust != 0.0) 
	{
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		adjust /= gen_ratio[lev-1][i];
	    }
	}
    }
    if (adjust != 0.0) 
    {
	adjust /= mg_domain[ml_index[lev_min]].numPts();
	
	if (pcode >= 2  && ParallelDescriptor::IOProcessor())
	{
	    cout << "HG: Sparse-source solvability adjustment: " << adjust << endl;
	}
	
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    source[lev].plus(-adjust, 0);
	}
    }
}

//
// This is a combination routine which combines sources from a divergence
// and from a cell-based right hand side S in the proper sequence.  The
// key feature is that both "grid" routines must be called before starting
// the lev_interface calculation, since they trash some lev_interface points.
//

void
holy_grail_amr_projector::right_hand_side (PArray<MultiFab>* u,
                                           PArray<MultiFab>& S)
{
    if (u) 
    {
	grid_divergence(u);
    }
    if (S.length() > 0) 
    {
	grid_average(S);
    }
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	const int mglev = ml_index[lev];
	
	clear_part_interface(source[lev], lev_interface[mglev]);
	
	if (lev > lev_min) 
	{
	    if (u) 
	    {
		interface_divergence(u, lev);
	    }
   
	    if (S.length() > 0) 
	    {
		interface_average(S, lev);
	    }
	}
    }
}

//
// Averages cell-based S onto node-based source conservatively
// across the composite mesh.  S must be passed in with a one
// cell wide border.  At inflow boundaries border values should
// be set to correct inflow condition.  Other border values passed
// in may be meaningless, but should not be NaNs.
//
// This routine will modify the borders of S.  Also, if the problem
// being solved is singular, S will be adjusted so that it integrates
// to 0 to maximum precision.  (It is assumed that any additional
// contribution to the right hand side will also integrate to 0.)
//
// This is an incomplete routine---interface_average must also be called.
//

void
holy_grail_amr_projector::grid_average (PArray<MultiFab>& S)
{
    BL_ASSERT(S[lev_min].nGrow() == 1);
    
    if (singular) 
    {
	Real adjust = 0.0;
	for (int lev = lev_max; lev > lev_min; lev--) 
	{
	    restrict_level(S[lev-1], S[lev], gen_ratio[lev-1], default_restrictor(), default_level_interface, 0);
	}
	for (MultiFabIterator S_mfi(S[lev_min]); S_mfi.isValid(); ++S_mfi)
	{
	    adjust += S_mfi->sum(S_mfi.validbox(), 0);
	}
	ParallelDescriptor::ReduceRealSum(adjust);
	adjust /= mg_domain[ml_index[lev_min]].numPts();
	
	if (pcode >= 2  && ParallelDescriptor::IOProcessor())
	{
	    cout << "HG: Cell-source solvability adjustment: " << adjust << endl;
	}
	
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    S[lev].plus(-adjust, 0);
	}
    }
    
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	const int mglev = ml_index[lev];
	
	fill_borders(S[lev], lev_interface[mglev], boundary.scalar(), -1, m_stencil == terrain);
	
	for (MultiFabIterator s_mfi(source[lev]); s_mfi.isValid(); ++s_mfi)
	{
	    DependentMultiFabIterator S_dmfi(s_mfi, S[lev]);
	    const Box& sbox = s_mfi->box();
	    const Box& fbox = S_dmfi->box();
	    const Box& freg = lev_interface[mglev].part_fine(s_mfi.index());
	    Real* sptr = s_mfi->dataPtr();
	    const Real* csptr = S_dmfi->dataPtr();
#if (BL_SPACEDIM == 2)
	    const Real hx = h[mglev][0];
	    const int isRZ = IsRZ();
	    const int imax = mg_domain[mglev].bigEnd(0) + 1;
	    FORT_HGAVG(sptr, DIMLIST(sbox), csptr, DIMLIST(fbox), DIMLIST(freg), &hx, &isRZ, &imax);
#else
	    FORT_HGAVG(sptr, DIMLIST(sbox), csptr, DIMLIST(fbox), DIMLIST(freg));
#endif
	}
    }
}

//
// This is an incomplete routine---interface_divergence must also be called.
//

void
holy_grail_amr_projector::grid_divergence (PArray<MultiFab>* u)
{ 
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	const int mglev = ml_index[lev];
	const Real hx = h[mglev][0];
	const Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
	const Real hz = h[mglev][2];
#endif
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    fill_borders(u[i][lev], lev_interface[mglev], boundary.velocity(i), -1, m_stencil == terrain);
	}
	
	for (MultiFabIterator s_mfi(source[lev]); s_mfi.isValid(); ++s_mfi)
	{
	    DependentMultiFabIterator u_dmfi_0(s_mfi, u[0][lev]);
	    DependentMultiFabIterator u_dmfi_1(s_mfi, u[1][lev]);
	    const Box& sbox = s_mfi->box();
	    const Box& fbox = u_dmfi_0->box();
	    const Box& freg = lev_interface[mglev].part_fine(s_mfi.index());
	    Real* sptr = s_mfi->dataPtr();
	    Real* u0ptr = u_dmfi_0->dataPtr();
	    Real* u1ptr = u_dmfi_1->dataPtr();
#if (BL_SPACEDIM == 3)
	    DependentMultiFabIterator u_dmfi_2(s_mfi, u[2][lev]);
	    Real* u2ptr = u_dmfi_2->dataPtr();
#endif
	    if (m_stencil == terrain)
	    {
		FORT_HGDIV_TERRAIN(sptr, DIMLIST(sbox), D_DECL(u0ptr, u1ptr, u2ptr), DIMLIST(fbox), DIMLIST(freg), D_DECL(&hx, &hy, &hz));
	    }
	    else
	    {
#if (BL_SPACEDIM == 2)
		const int isRZ = IsRZ();
		const int imax = mg_domain[mglev].bigEnd(0) + 1;
		FORT_HGDIV(sptr, DIMLIST(sbox), D_DECL(u0ptr, u1ptr, u2ptr), DIMLIST(fbox), DIMLIST(freg), D_DECL(&hx, &hy, &hz), &isRZ, &imax);
#else
		FORT_HGDIV(sptr, DIMLIST(sbox), D_DECL(u0ptr, u1ptr, u2ptr), DIMLIST(fbox), DIMLIST(freg), D_DECL(&hx, &hy, &hz), 0, 0);
#endif
	    }
	}
    }
}

//
// Obsolete:
//

void
holy_grail_amr_projector::sync_right_hand_side (PArray<MultiFab>* u)
{ 
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	source[lev].setVal(0.0);
    }
    
    const int mglev0 = ml_index[lev_min];
    interface_divergence(u, lev_min+1);
    
    if (singular) 
    {
	const int mglev1 = ml_index[lev_min+1];
	restrict_level(source[lev_min], source[lev_min+1], gen_ratio[lev_min], 
	    bilinear_restrictor_class(0, m_stencil == terrain), lev_interface[mglev1], mg_boundary);
	work[mglev0].setVal(1.0);
	Real adjustment = inner_product(source[lev_min], work[mglev0]) /
	    mg_domain[ml_index[lev_min]].volume();
	if (pcode >= 2  && ParallelDescriptor::IOProcessor())
	{
	    cout << "HG: Solvability adjustment is " << adjustment << endl;
	}
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    source[lev].plus(-adjustment, 0);
	}
    }
}

void
holy_grail_amr_projector::interface_average (PArray<MultiFab>& S, int lev)
{ 
    const int mglev = ml_index[lev];
    const int mgc = ml_index[lev-1];
    
    const IntVect& rat = gen_ratio[lev-1];
    task_list tl;
    for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
    {
        //
	// Find a fine grid touching this face.
        //
	int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	if (igrid < 0)
	    igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
	    continue;
        //
	// Fine grid on just one side.
        //
	const int idim = lev_interface[mglev].fdim(iface);
	const int idir = (geo & level_interface::LOW) ? -1 : 1;
	Box cbox = lev_interface[mglev].box(level_interface::FACEDIM, iface);
	const IntVect t = cbox.type();
	if (idir > 0)
	    cbox.growLo(idim, rat[idim]);
	else
	    cbox.growHi(idim, rat[idim]);
	cbox.convert(IntVect::TheCellVector()).coarsen(rat);
	task_fab* Scp = new task_fill_patch(tl, source[lev], igrid, cbox, S[lev-1], lev_interface[mgc], boundary.scalar(), -1, -1);
	Box creg = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	creg.coarsen(rat).grow(t - 1);
#if (BL_SPACEDIM == 2)
	const Real hx = h[mglev][0];
	const int isRZ = IsRZ();
	const int imax = mg_domain[mglev].bigEnd(0) + 1;
	tl.add_task(new task_fecavg(&FORT_HGFAVG, tl, source[lev], S[lev], igrid, Scp, creg, rat, idim, idir, hx, isRZ, imax));
#else
	tl.add_task(new task_fecavg(&FORT_HGFAVG, tl, source[lev], S[lev], igrid, Scp, creg, rat, idim, idir));
#endif
    }
    tl.execute("holy_grail_amr_projector::interface_average(1)");
#if (BL_SPACEDIM == 3)
    for (int iedge = 0; iedge < lev_interface[mglev].nboxes(1); iedge++) 
    {
        //
	// Find a fine grid touching this edge.
        //
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(1); i++) 
	{
	    igrid = lev_interface[mglev].grid(1, iedge, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(1, iedge);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(1, iedge) )
	    continue;
        //
	// Fine grid on just one side.
        //
	Box cbox = lev_interface[mglev].box(1, iedge);
	const IntVect t = cbox.type();
	cbox.coarsen(rat).grow(t).convert(IntVect::TheCellVector());
	Box fbox = cbox;
	fbox.refine(rat);
	task_fab* Scp = new task_fill_patch(tl, source[lev], igrid, cbox, S[lev-1], lev_interface[mgc],   boundary.scalar(), -1, -1);
	task_fab* Sfp = new task_fill_patch(tl, source[lev], igrid, fbox, S[lev],   lev_interface[mglev], boundary.scalar(), 1, iedge);
	Box creg = lev_interface[mglev].node_box(1, iedge);
	creg.coarsen(rat).grow(t - 1);
	Array<int> ga = lev_interface[mglev].geo_array(1, iedge);
	task::task_proxy tp = tl.add_task(new task_fecavg_2(&FORT_HGEAVG, tl, source[lev], igrid, Sfp, Scp, creg, rat, ga, t));
        //
	// Fill in the grids on the other sides, if any.
        //
	const Box& freg = lev_interface[mglev].node_box(1, iedge);
	for (int i = 1; i < lev_interface[mglev].ngrids(1); i++) 
	{
	    const int jgrid = lev_interface[mglev].grid(1, iedge, i);
	    if (jgrid >= 0 && jgrid != igrid)
	    {
                tl.add_task(new task_copy(tl,source[lev],jgrid,source[lev],igrid,freg,tp));
	    }
	}
    }
    tl.execute("holy_grail_amr_projector::interface_average(2)");
#endif
    for (int icor = 0; icor < lev_interface[mglev].nboxes(0); icor++) 
    {
        //
	// Find a fine grid touching this corner.
        //
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(0); i++) 
	{
	    igrid = lev_interface[mglev].grid(0, icor, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(0, icor);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(0, icor) )
	    continue;
        //
	// Fine grid on just one side.
        //
	Box cbox = lev_interface[mglev].box(0, icor);
	cbox.coarsen(rat).grow(1).convert(IntVect::TheCellVector());
	Box fbox = cbox;
	fbox.refine(rat);
	task_fab* Scp = new task_fill_patch(tl, source[lev], igrid, cbox, S[lev-1], lev_interface[mgc],   boundary.scalar(), -1, -1);
	task_fab* Sfp = new task_fill_patch(tl, source[lev], igrid, fbox, S[lev],   lev_interface[mglev], boundary.scalar(), 0, icor);
	Box creg = lev_interface[mglev].box(0, icor);
	creg.coarsen(rat);
	Array<int> ga = lev_interface[mglev].geo_array(0, icor);
#if (BL_SPACEDIM == 2)
	const Real hx = h[mglev][0];
	const int isRZ = IsRZ();
	const int imax = mg_domain[mglev].bigEnd(0) + 1;
	task::task_proxy tp =
	  tl.add_task(new task_fecavg_2(&FORT_HGCAVG, tl, source[lev], igrid, Sfp, Scp, creg, rat, ga, IntVect(), hx, isRZ, imax));
#else
	task::task_proxy tp =
	  tl.add_task(new task_fecavg_2(&FORT_HGCAVG, tl, source[lev], igrid, Sfp, Scp, creg, rat, ga, IntVect()));
#endif
        //
	// Fill in the grids on the other sides, if any.
        //
	const Box& freg = lev_interface[mglev].box(0, icor);
	for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
	{
	    const int jgrid = lev_interface[mglev].grid(0, icor, i);
	    if (jgrid >= 0 && jgrid != igrid)
	    {
                tl.add_task(new task_copy(tl,source[lev],jgrid,source[lev],igrid,freg,tp));
	    }
	}
    }
    tl.execute("holy_grail_amr_projector::interface_average(3)");
}

void
holy_grail_amr_projector::interface_divergence (PArray<MultiFab>* u,
                                                int               lev)
{
    const int mglev = ml_index[lev];
    const int mgc = ml_index[lev-1];
    
    const IntVect& rat = gen_ratio[lev-1];
    task_list tl;
    HG_TEST_NORM(source[lev], "interface_divergence,b");
    for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
    {
        //
	// Find a fine grid touching this face.
        //
	int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	if (igrid < 0)
	    igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
	    continue;
        //
	// Fine grid on just one side.
        //
	const int idim = lev_interface[mglev].fdim(iface);
	const int idir = (geo & level_interface::LOW) ? -1 : 1;
	Box cbox = lev_interface[mglev].box(level_interface::FACEDIM, iface);
	const IntVect t = cbox.type();
	if (idir > 0)
	    cbox.growLo(idim, rat[idim]);
	else
	    cbox.growHi(idim, rat[idim]);
	cbox.convert(IntVect::TheCellVector()).coarsen(rat);
	task_fab* ucp[BL_SPACEDIM];
	MultiFab* upt[BL_SPACEDIM];
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    ucp[i] = new task_fill_patch(tl, source[lev], igrid, cbox, u[i][lev-1], lev_interface[mgc], boundary.velocity(i), -1, -1);
	    upt[i] = &u[i][lev];
	}
	Box creg = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	creg.coarsen(rat).grow(t - 1);
	if (m_stencil == terrain)
	{
	    BoxLib::Error("Not Ready");
	    // tl.add_task(new task_fecdiv(&FORT_HGFDIV_TERRAIN, tl, source[lev], upt, igrid, ucp, creg, h[mglev], rat, idim, idir));
	}
	else
	{
#if (BL_SPACEDIM == 2)
	    const int isRZ = IsRZ();
	    const int imax = mg_domain[mglev].bigEnd(0) + 1;
	    tl.add_task(new task_fecdiv(&FORT_HGFDIV, tl, source[lev], upt, igrid, ucp, creg, h[mglev], rat, idim, idir, isRZ, imax));
#else
	    tl.add_task(new task_fecdiv(&FORT_HGFDIV, tl, source[lev], upt, igrid, ucp, creg, h[mglev], rat, idim, idir));
#endif
	}
    }
    tl.execute("holy_grail_amr_projector::interface_divergence(1)");

    HG_TEST_NORM(source[lev], "interface_divergence,a");
    
#if (BL_SPACEDIM == 3)
    
    for (int iedge = 0; iedge < lev_interface[mglev].nboxes(1); iedge++) 
    {
        //
	// Find a fine grid touching this edge.
        //
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(1); i++) 
	{
	    igrid = lev_interface[mglev].grid(1, iedge, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(1, iedge);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(1, iedge) )
	    continue;
        //
	// Fine grid on just one side.
        //
	Box cbox = lev_interface[mglev].box(1, iedge);
	const IntVect t = cbox.type();
	cbox.coarsen(rat).grow(t).convert(IntVect::TheCellVector());
	Box fbox = cbox;
	fbox.refine(rat);
	task_fab* ucp[BL_SPACEDIM];
	task_fab* ufp[BL_SPACEDIM];
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    ucp[i] = new task_fill_patch(tl, source[lev], igrid, cbox, u[i][lev-1], lev_interface[mgc],   boundary.velocity(i), -1, -1);
	    ufp[i] = new task_fill_patch(tl, source[lev], igrid, fbox, u[i][lev],   lev_interface[mglev], boundary.velocity(i), 1, iedge);
	}
	Box creg = lev_interface[mglev].node_box(1, iedge);
	creg.coarsen(rat).grow(t - 1);
	Array<int> ga = lev_interface[mglev].geo_array(1, iedge);
	task::task_proxy tp;
	if (m_stencil == terrain)
	{
	    tp = tl.add_task(new task_fecdiv_2(&FORT_HGEDIV_TERRAIN, tl, source[lev], igrid, ufp, ucp, creg, h[mglev], rat, ga, t));
	}
	else
	{
	    tp = tl.add_task(new task_fecdiv_2(&FORT_HGEDIV,         tl, source[lev], igrid, ufp, ucp, creg, h[mglev], rat, ga, t));
	}
        //
	// Fill in the grids on the other sides, if any.
        //
	const Box& freg = lev_interface[mglev].node_box(1, iedge);
	for (int i = 1; i < lev_interface[mglev].ngrids(1); i++) 
	{
	    const int jgrid = lev_interface[mglev].grid(1, iedge, i);
	    if (jgrid >= 0 && jgrid != igrid)
	    {
                tl.add_task(new task_copy(tl,source[lev],jgrid,source[lev],igrid,freg,tp));
	    }
	}
    }
    tl.execute("holy_grail_amr_projector::interface_divergence(2)");    
    HG_TEST_NORM(source[lev], "interface_divergence,a1");
#endif
    for (int icor = 0; icor < lev_interface[mglev].nboxes(0); icor++) 
    {
        //
	// Find a fine grid touching this corner.
        //
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(0); i++) 
	{
	    igrid = lev_interface[mglev].grid(0, icor, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(0, icor);
        //
	// Reject fine-fine interfaces and those without an interior fine grid
        //
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(0, icor) )
	    continue;
        //
	// Fine grid on just one side.
        //
	Box cbox = lev_interface[mglev].box(0, icor);
	cbox.coarsen(rat).grow(1).convert(IntVect::TheCellVector());
	Box fbox = cbox;
	fbox.refine(rat);
	task_fab* ucp[BL_SPACEDIM];
	task_fab* ufp[BL_SPACEDIM];
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    ucp[i] = new task_fill_patch(tl, source[lev], igrid, cbox, u[i][lev-1], lev_interface[mgc],   boundary.velocity(i), -1, -1);
	    ufp[i] = new task_fill_patch(tl, source[lev], igrid, fbox, u[i][lev],   lev_interface[mglev], boundary.velocity(i), 0, icor);
	}
	Box creg = lev_interface[mglev].box(0, icor);
	creg.coarsen(rat);
	Array<int> ga = lev_interface[mglev].geo_array(0, icor);
	task::task_proxy tp;
	if (m_stencil == terrain)
	{
	    tp = tl.add_task(new task_fecdiv_2(&FORT_HGCDIV_TERRAIN, tl, source[lev], igrid, ufp, ucp, creg, h[mglev], rat, ga));
	}
	else
	{
#if (BL_SPACEDIM == 2)
            const int isRZ = IsRZ();
	    tp = tl.add_task(new task_fecdiv_3(&FORT_HGCDIV,         tl, source[lev], igrid, ufp, ucp, creg, h[mglev], rat, ga, isRZ));
#elif (BL_SPACEDIM == 3)
	    tp = tl.add_task(new task_fecdiv_2(&FORT_HGCDIV,         tl, source[lev], igrid, ufp, ucp, creg, h[mglev], rat, ga));
#endif
	}
        //
	// Fill in the grids on the other sides, if any.
        //
	const Box& freg = lev_interface[mglev].box(0, icor);
	for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
	{
	    const int jgrid = lev_interface[mglev].grid(0, icor, i);
	    if (jgrid >= 0 && jgrid != igrid)
	    {
                tl.add_task(new task_copy(tl,source[lev],jgrid,source[lev],igrid,freg,tp));
	    }
	}
    }
    tl.execute("holy_grail_amr_projector::interface_divergence(3)");
    HG_TEST_NORM(source[lev], "interface_divergence,a2");
}

void
holy_grail_amr_projector::form_solution_vector (PArray<MultiFab>*       u,
                                                const PArray<MultiFab>& sigma_in)
{
    BL_ASSERT(u != 0);

    if (u)
    {
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    int mglev = ml_index[lev];
	    Real hxyz[BL_SPACEDIM] = { D_DECL( h[mglev][0], h[mglev][1], h[mglev][2] ) };
            FArrayBox gp[BL_SPACEDIM];
	    for (MultiFabIterator d_mfi(dest[lev]); d_mfi.isValid(); ++d_mfi)
	    {
		const Box& gbox = ml_mesh[lev][d_mfi.index()];
		const Box& dbox = d_mfi->box();
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    gp[i].resize(gbox);
		}
		if (m_stencil == terrain)
		{
		    FORT_HGGRAD_TERRAIN(D_DECL(gp[0].dataPtr(), gp[1].dataPtr(), gp[2].dataPtr()), DIMLIST(gbox),
					d_mfi->dataPtr(), DIMLIST(dbox),
					DIMLIST(gbox), D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]));
		}
		else
		{
#if (BL_SPACEDIM == 2)
		    const int isRZ = IsRZ();
		    FORT_HGGRAD(D_DECL(gp[0].dataPtr(), gp[1].dataPtr(), gp[2].dataPtr()), DIMLIST(gbox),
				d_mfi->dataPtr(), DIMLIST(dbox),
				DIMLIST(gbox), D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]), &isRZ);
#else
		    FORT_HGGRAD(D_DECL(gp[0].dataPtr(), gp[1].dataPtr(), gp[2].dataPtr()), DIMLIST(gbox),
				d_mfi->dataPtr(), DIMLIST(dbox),
				DIMLIST(gbox), D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]), 0);
#endif
		}		
		if (m_stencil != terrain)
		{
		    DependentMultiFabIterator s_dmfi(d_mfi, sigma_in[lev]);
		    for (int i = 0; i < BL_SPACEDIM; i++) 
		    {
			gp[i].mult(s_dmfi());
			DependentMultiFabIterator u_dmfi(d_mfi, u[i][lev]);
			u_dmfi->minus(gp[i]);
		    }
		}
		else
		{
		    DependentMultiFabIterator s_dmfi(d_mfi, sigma_in[lev]);
		    FArrayBox cross(gbox);
		    for (int i = 0; i < BL_SPACEDIM; i++) 
		    {
			cross.copy(gp[i]);
			cross.mult(s_dmfi(), i, 0, 1);
			DependentMultiFabIterator u_dmfi(d_mfi, u[i][lev]);
			u_dmfi->minus(cross);
		    }
		    DependentMultiFabIterator ul_dmfi(d_mfi, u[BL_SPACEDIM-1][lev]);
		    for (int i = 0; i < BL_SPACEDIM - 1; i++) 
		    {
			cross.copy(gp[BL_SPACEDIM-1]);
			cross.mult(s_dmfi(), BL_SPACEDIM+i, 0, 1);
			DependentMultiFabIterator u_dmfi(d_mfi, u[i][lev]);
			u_dmfi->plus(cross);
			cross.copy(gp[i]);
			cross.mult(s_dmfi(), BL_SPACEDIM+i, 0, 1);
			ul_dmfi->plus(cross);
		    }
		}
	    }
	}
	
	sync_periodic_interfaces();
	for (int lev = lev_max; lev > lev_min; lev--) 
	{
	    const IntVect& rat = gen_ratio[lev-1];
	    restrict_level(dest[lev-1], dest[lev], rat, injection_restrictor_class(), default_level_interface, 0);
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		if (m_stencil != terrain)
		{
		    restrict_level(u[i][lev-1], u[i][lev], rat, default_restrictor(), default_level_interface, 0);
		}
		else
		{
		    restrict_level(u[i][lev-1], u[i][lev], rat, terrain_velocity_restrictor_class(i), default_level_interface, 0);
		}
	    }
	}
    }
    else 
    {
	sync_periodic_interfaces();
	for (int lev = lev_max; lev > lev_min; lev--) 
	{
	    restrict_level(dest[lev-1], dest[lev], gen_ratio[lev-1], injection_restrictor_class(), default_level_interface, 0);
	}
    }
}
