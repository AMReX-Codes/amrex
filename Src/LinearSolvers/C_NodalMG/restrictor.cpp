//BL_COPYRIGHT_NOTICE

#include "restrictor.H"
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FACRST1  acrst1_
#define FORT_FANRST1  anrst1_
#define FORT_FANRST2  anrst2_
#define FORT_FANFR2   anfr2_
#define FORT_FANER2   aner2_
#define FORT_FANCR2   ancr2_
#define FORT_FANOR2   anor2_
#define FORT_FANIR2   anir2_
#define FORT_FANDR2   andr2_
#else
#define FORT_FACRST1  ACRST1
#define FORT_FANRST1  ANRST1
#define FORT_FANRST2  ANRST2
#define FORT_FANFR2   ANFR2
#define FORT_FANER2   ANER2
#define FORT_FANCR2   ANCR2
#define FORT_FANOR2   ANOR2
#define FORT_FANIR2   ANIR2
#define FORT_FANDR2   ANDR2
#endif

extern "C" 
{
    void FORT_FACRST1(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANRST1(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    //
    // Used in the parallel loops, most of these routines have bogus elements
    // in their calling sequences.
    //
    void FORT_FANRST2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANFR2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANER2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANCR2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
#if BL_SPACEDIM == 2
    void FORT_FANOR2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANIR2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANDR2 (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
#endif
}

class task_fab_get : public task_fab
{
public:
    task_fab_get (task_list&      tl_,
                  const MultiFab& d_,
                  int             dgrid_,
                  const Box&      bx,
                  const MultiFab& s_,
                  int             sgrid_);

    virtual ~task_fab_get ();
    virtual bool work_to_do () const;
private:
    //
    // The data.
    //
    task::task_proxy tf;
    const MultiFab&  s;
    const int        sgrid;
    const Box        bx;
};

task_fab_get::task_fab_get (task_list&      tl_,
                            const MultiFab& d_,
                            int             dgrid_,
                            const Box&      bx_,
                            const MultiFab& s_,
                            int             sgrid_) 
    :
    task_fab(tl_,d_,dgrid_,bx_,s_.nComp()),
    s(s_),
    sgrid(sgrid_),
    bx(bx_),
    tf(0)
{
    depend_on(tf = m_task_list.add_task(new task_copy_local(m_task_list,target,target_proc_id(),bx,s,sgrid)));
}

task_fab_get::~task_fab_get () {}

bool
task_fab_get::work_to_do () const
{
    return ParallelDescriptor::MyProc() == m_target_proc_id || !tf.null();
}
 
typedef void (*RESTFUN)(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*); 

struct task_restriction_fill : public task
{
    task_restriction_fill (const RESTFUN  ref_,
                           task_list&     tl_,
                           MultiFab&      m_,
                           int            ind_,
                           const Box&     cbox_,
                           task_fab*      tf_,
                           const IntVect& rat_,
                           int            integrate_,
                           int            i1_ = 0,
                           int            i2_ = 0);

    task_restriction_fill (const RESTFUN     ref_,
                           task_list&        tl_,
                           MultiFab&         m_,
                           int               ind_,
                           const Box&        cbox_,
                           task_fab*         tf_,
                           const IntVect&    rat_,
                           int               integrate_,
                           const Array<int>& i1_);

    task_restriction_fill (const RESTFUN     ref_,
                           task_list&        tl_,
                           MultiFab&         m_,
                           int               ind_,
                           const Box&        cbox_,
                           task_fab*         tf_,
                           const IntVect&    rat_,
                           int               integrate_,
                           const IntVect&    i1_,
                           const Array<int>& i2_);

    virtual bool ready ();
    virtual bool work_to_do () const;

private:
    //
    // The data.
    //
    const RESTFUN ref;
    task_proxy    tf;
    MultiFab&     m;
    int           ind;
    const Box     cbox;
    const IntVect rat;
    const int     integrate;
    Array<int>    arg1;
    Array<int>    arg2;
};

task_restriction_fill::task_restriction_fill (const RESTFUN  ref_,
                                              task_list&     tl_,
                                              MultiFab&      m_,
                                              int            ind_,
                                              const Box&     cbox_,
                                              task_fab*      tf_,
                                              const IntVect& rat_,
                                              int            integrate_,
                                              int            i1_,
                                              int            i2_) 
    :
    task(tl_),
    ref(ref_),
    m(m_),
    ind(ind_),
    cbox(cbox_),
    rat(rat_),
    integrate(integrate_),
    arg1(1),
    arg2(1) 
{
    depend_on(tf = m_task_list.add_task(tf_));
    arg1[0] = i1_;
    arg2[0] = i2_;
}

task_restriction_fill::task_restriction_fill (const RESTFUN     ref_,
                                              task_list&        tl_,
                                              MultiFab&         m_,
                                              int               ind_,
                                              const Box&        cbox_,
                                              task_fab*         tf_,
                                              const IntVect&    rat_,
                                              int               integrate_,
                                              const Array<int>& i1_) 
    :
    task(tl_),
    ref(ref_),
    m(m_),
    ind(ind_),
    cbox(cbox_),
    rat(rat_),
    integrate(integrate_),
    arg1(i1_),
    arg2(1) 
{
    depend_on(tf = m_task_list.add_task(tf_));
    arg2[0] = 0;
}

task_restriction_fill::task_restriction_fill (const RESTFUN     ref_,
                                              task_list&        tl_,
                                              MultiFab&         m_,
                                              int               ind_,
                                              const Box&        cbox_,
                                              task_fab*         tf_,
                                              const IntVect&    rat_,
                                              int               integrate_,
                                              const IntVect&    i1_,
                                              const Array<int>& i2_) 
    :
    task(tl_),
    ref(ref_),
    m(m_),
    ind(ind_),
    cbox(cbox_),
    rat(rat_),
    integrate(integrate_),
    arg1(i1_.getVect(),BL_SPACEDIM),
    arg2(i2_) 
{
    depend_on(tf = m_task_list.add_task(tf_));
}

bool
task_restriction_fill::ready ()
{
    assert(!tf.null());
    assert(tf->ready());
    if (is_local(m, ind))
    {
        task_fab* tff = dynamic_cast<task_fab*>(tf.get());
        assert(tff != 0);
        const Box& fb = tff->fab().box();
        const Box& pb = m[ind].box();
        (*ref)(m[ind].dataPtr(), DIMLIST(pb), DIMLIST(cbox), tff->fab().dataPtr(), DIMLIST(fb), D_DECL(rat[0], rat[1], rat[2]), m.nComp(), &integrate, arg1.dataPtr(), arg2.dataPtr());
    }
    return true;
}

bool
task_restriction_fill::work_to_do () const
{
    return is_local(m,ind) || !tf.null();
}

amr_restrictor_class::~amr_restrictor_class () {}

Box
amr_restrictor_class::box (const Box&     fb,
                           const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat);
}

Box
amr_restrictor_class::rebox (const Box&     cb,
                             const IntVect& rat) const
{
    Box retbox(cb);
    return retbox.refine(rat);
}

void
amr_restrictor_class::fill_interface (MultiFab&,
                                      MultiFab&,
                                      const level_interface&,
                                      const amr_boundary_class*,
                                      const IntVect&) const 
{}

void
cell_average_restrictor_class::fill (FArrayBox&       patch,
                                     const Box&       region,
                                     const FArrayBox& fgr,
                                     const IntVect&   rat) const
{
    assert(patch.box().cellCentered());
    assert(patch.nComp() == 1);
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
    D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), &integrate, 0, 0);
}

void
terrain_velocity_restrictor_class::fill (FArrayBox&       patch,
                                         const Box&       region,
                                         const FArrayBox& fgr,
                                         const IntVect&   rat) const
{
    assert(patch.box().cellCentered());
    assert(patch.nComp() == 1);
    const int integ = 1;
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
    D_DECL(rat[0], rat[1], rat[2]), 1, &integ, 0, 0);
    Real fac = 1.0 / rat[integrate];
    patch.mult(fac, region);
}

void
injection_restrictor_class::fill (FArrayBox&       patch,
                                  const Box&       region,
                                  const FArrayBox& fgr,
                                  const IntVect&   rat) const
{
    assert(patch.box().type() == IntVect::TheNodeVector());
    assert(patch.nComp() == fgr.nComp());
    FORT_FANRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
    D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), 0, 0, 0);
}

void
default_restrictor::fill (FArrayBox&       patch,
                          const Box&       region,
                          const FArrayBox& fgr,
                          const IntVect&   rat) const
{
    assert(patch.box().cellCentered() || patch.box().type() == IntVect::TheNodeVector());
    assert(patch.nComp() == fgr.nComp());
    if (patch.box().cellCentered())
    {
        cell_average_restrictor_class(0).fill(patch, region, fgr, rat);
    }
    else if (patch.box().type() == IntVect::TheNodeVector())
    {
        injection_restrictor_class().fill(patch, region, fgr, rat);
    }
}

bilinear_restrictor_class::bilinear_restrictor_class (int i,
                                                      bool hg_terrain)
    :
    integrate(i),
    m_hg_terrain(hg_terrain)
{
    assert(i == 0 || i == 1);
}

Box
bilinear_restrictor_class::box (const Box&     fb,
                                const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat).grow(-1);
}

Box
bilinear_restrictor_class::rebox (const Box&     cb,
                                  const IntVect& rat) const
{
    Box retbox(cb);
    return retbox.refine(rat).grow(rat-1);
}

void
bilinear_restrictor_class::fill (FArrayBox&       patch,
                                 const Box&       region,
                                 const FArrayBox& fgr,
                                 const IntVect&   rat) const
{
    assert(patch.box().type() == IntVect::TheNodeVector());
    assert(patch.nComp() == fgr.nComp());
    FORT_FANRST2(patch.dataPtr(), DIMLIST(patch.box()), 
                 DIMLIST(region), 
                 fgr.dataPtr(), DIMLIST(fgr.box()),
                 D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), &integrate, 0, 0);
}

void
bilinear_restrictor_class::fill_interface (MultiFab&                 dest,
                                           MultiFab&                 fine,
                                           const level_interface&    lev_interface,
                                           const amr_boundary_class* bdy,
                                           const IntVect&            rat) const
{
    assert(type(dest) == IntVect::TheNodeVector());
    assert(dest.nComp() == fine.nComp());

    int ratmax = rat[0];
    for (int i = 1; i < BL_SPACEDIM; ++i)
        ratmax = (rat[i] > ratmax) ? rat[i] : ratmax;
    
    if (fine.nGrow() >= ratmax - 1) 
        fill_borders(fine, lev_interface, bdy, ratmax - 1, m_hg_terrain);

    task_list tl;
    for (int jgrid = 0; jgrid < dest.length(); jgrid++)
    {
        const Box& region = dest.box(jgrid);
        //
        // Interface restriction is sufficiently rare and specialized that
        // we will let the restrictor handle it---at least for now.
        //
        // This assertion difficult in BoxLib since r.mesh() is not cc:
        //assert(r.mesh() == lev_interface.interior_mesh());
        //
        const Box regplus = ::grow(region,1);
    
        for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
        {
            if (lev_interface.flag(level_interface::FACEDIM, iface))
                continue;
            Box cbox = lev_interface.node_box(level_interface::FACEDIM, iface);
            const IntVect t = lev_interface.box(level_interface::FACEDIM, iface).type();
            const unsigned int geo = lev_interface.geo(level_interface::FACEDIM, iface);
            cbox.coarsen(rat);
            if (region.intersects(cbox)) 
            {
                //
                // This extends fine face by one coarse cell past coarse face:
                //
                cbox &= regplus;
                int idim = lev_interface.fdim(iface);
                cbox.grow(t - 1);
                if (geo == level_interface::ALL) 
                {
                    //
                    // Fine grid on both sides.
                    //
                    if (fine.nGrow() >= ratmax - 1) 
                    {
                        int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
                        if (igrid < 0)
                            igrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
                        //
                        // FIXME--want minimal box.
                        //
                        const Box& fb = ::grow(fine.box(igrid), fine.nGrow());
                        assert(is_remote(fine,igrid) || fb == fine[igrid].box());
                        // const Box& fb = fine[igrid].box();
                        task_fab* tfab = new task_fab_get(tl,dest,jgrid,fb,fine,igrid);
                        tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
                    }
                    else
                    {
                        Box fbox = ::grow(refine(cbox, rat), rat - 1);
                        task_fab* tfab = new task_fill_patch(tl,dest,jgrid,fbox,fine,lev_interface,bdy,level_interface::FACEDIM,iface);
                        tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
                    }
                }
                else 
                {
                    //
                    // Fine grid on just one side.
                    //
                    const int idir = (geo & level_interface::LOW) ? -1 : 1;
                    const int igrid = (idir < 0) ? lev_interface.grid(level_interface::FACEDIM, iface, 0) :
                        lev_interface.grid(level_interface::FACEDIM, iface, 1) ;
                    if (igrid >= 0) 
                    {
                        //
                        // Usual case, a fine grid extends all along the face.
                        //
                        const Box& fb = grow(fine.box(igrid), fine.nGrow());
                        assert ( is_remote(fine, igrid) || fb == fine[igrid].box());
                        // const Box& fb = fine[igrid].box();
                        task_fab* tfab = new task_fab_get(tl,dest,jgrid,fb,fine,igrid);
                        tl.add_task(new task_restriction_fill(&FORT_FANFR2,tl,dest,jgrid,cbox,tfab,rat,integrate,idim,idir));
                    }
                    else 
                    {
                        //
                        // A virtual fine grid is on the other side of the boundary.
                        //
                        Box fbox = ::refine(cbox, rat).grow(rat - 1);
                        if (geo & level_interface::LOW)
                            fbox.growHi(idim, 1 - rat[idim]);
                        else
                            fbox.growLo(idim, 1 - rat[idim]);
                        task_fab* tfab = new task_fill_patch(tl,dest,jgrid,fbox,fine,lev_interface,bdy,level_interface::FACEDIM,iface);
                        tl.add_task(new task_restriction_fill(&FORT_FANFR2,tl,dest,jgrid,cbox,tfab,rat,integrate,idim,idir));
                    }
                }
            }
        }
#if (BL_SPACEDIM == 3)
    
    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
    {
        if (lev_interface.flag(1, iedge))
            continue;
        Box cbox = lev_interface.node_box(1, iedge);
        const IntVect t = lev_interface.box(1, iedge).type();
        cbox.coarsen(rat);
        if (region.intersects(cbox)) 
        {
            //
            // This extends fine edge by one coarse cell past coarse face:
            //
            cbox &= regplus;
            cbox.grow(t - 1);
            const unsigned int geo = lev_interface.geo(1, iedge);
            if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
            {
                int igrid = lev_interface.grid(1, iedge, 0);
                for (int itmp = 1; igrid < 0; itmp++)
                    igrid = lev_interface.grid(1, iedge, itmp);
                Box fb = ::grow(fine.box(igrid), fine.nGrow());
                assert(is_remote(fine, igrid) || fb == fine[igrid].box());
                // const Box& fb = fine[igrid].box();
                task_fab* tfab = new task_fab_get(tl,dest,jgrid,fb,fine,igrid);
                tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
            }
            else 
            {
                Box fbox = ::grow(refine(cbox, rat), rat - 1);
                task_fab* tfab = new task_fill_patch(tl,dest,jgrid,fbox,fine,lev_interface,bdy,1,iedge);
                if (geo == level_interface::ALL) 
                { 
                    // fine grid on all sides
                    tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
                }
                else 
                {
                    Array<int> ga = lev_interface.geo_array(1, iedge);
                    tl.add_task(new task_restriction_fill(&FORT_FANER2,tl,dest,jgrid,cbox,tfab,rat,integrate,t,ga));
                }
            }
        }
    }
#endif
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
        if (lev_interface.flag(0, icor))
            continue;
        Box cbox = lev_interface.box(0, icor);
        cbox.coarsen(rat);
        if (region.intersects(cbox)) 
        {
            const unsigned int geo = lev_interface.geo(0, icor);
            if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
            {
                int igrid = lev_interface.grid(0, icor, 0);
                for (int itmp = 1; igrid < 0; itmp++)
                    igrid = lev_interface.grid(0, icor, itmp);
                Box fb = ::grow(fine.box(igrid), fine.nGrow());
                assert( is_remote(fine, igrid) || fb == fine[igrid].box());
                // const Box& fb = fine[igrid].box();
                task_fab* tfab = new task_fab_get(tl,dest,jgrid,fb,fine,igrid);
                tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
            }
            else 
            {
                Box fbox = ::grow(refine(cbox, rat), rat - 1);
                task_fab* tfab = new task_fill_patch(tl,dest,jgrid,fbox,fine,lev_interface,bdy,0,icor);
                if (geo == level_interface::ALL) 
                { 
                    // fine grid on all sides
                    tl.add_task(new task_restriction_fill(&FORT_FANRST2,tl,dest,jgrid,cbox,tfab,rat,integrate));
                }
                else 
                {
                    Array<int> ga = lev_interface.geo_array(0, icor);
                    tl.add_task(new task_restriction_fill(&FORT_FANCR2,tl,dest,jgrid,cbox,tfab,rat,integrate,ga));
                }
            }
        }
    }
}
tl.execute();
}
