
#include "restrictor.H"
#include "fill_patch.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define FORT_FACRST1  acrst1_
#define FORT_FANRST1  anrst1_
#define FORT_FANRST2  anrst2_
#define FORT_FANFR2   anfr2_
#define FORT_FANER2   aner2_
#define FORT_FANCR2   ancr2_
#define FORT_FANOR2   anor2_
#define FORT_FANIR2   anir2_
#define FORT_FANDR2   andr2_
#elif defined( BL_FORT_USE_UPPERCASE )
#define FORT_FACRST1  ACRST1
#define FORT_FANRST1  ANRST1
#define FORT_FANRST2  ANRST2
#define FORT_FANFR2   ANFR2
#define FORT_FANER2   ANER2
#define FORT_FANCR2   ANCR2
#define FORT_FANOR2   ANOR2
#define FORT_FANIR2   ANIR2
#define FORT_FANDR2   ANDR2
#elif defined( BL_FORT_USE_LOWERCASE )
#define FORT_FACRST1  acrst1
#define FORT_FANRST1  anrst1
#define FORT_FANRST2  anrst2
#define FORT_FANFR2   anfr2
#define FORT_FANER2   aner2
#define FORT_FANCR2   ancr2
#define FORT_FANOR2   anor2
#define FORT_FANIR2   anir2
#define FORT_FANDR2   andr2
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_FACRST1(Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANRST1(Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    //
    // Used in the parallel loops, most of these routines have bogus elements
    // in their calling sequences.
    //
    void FORT_FANRST2(Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANFR2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANER2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANCR2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
#if BL_SPACEDIM == 2
    void FORT_FANOR2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANIR2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
    void FORT_FANDR2 (Real*, intS, intS, const Real*, intS, intRS,
		      const int&, const int*, const int*, const int*);
#endif
}

class task_fab_get
    :
    public task_fab
{
public:

    task_fab_get (task_list&      tl_,
                  const MultiFab& d_,
                  int             dgrid_,
                  const Box&      bx,
                  const MultiFab& s_,
                  int             sgrid_);
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
    task_fab(tl_, d_, dgrid_, bx_, s_.nComp()),
    s(s_),
    sgrid(sgrid_),
    bx(bx_),
    tf(0)
{
    depend_on(tf = m_task_list.add_task(new task_copy_local(m_task_list,
                                                            target,
                                                            target_proc_id(),
                                                            bx,
                                                            s,
                                                            sgrid)));
}

extern "C"
{
  typedef void (*RESTFUN)(Real*, intS, intS, const Real*, intS, intRS,
			  const int&, const int*, const int*, const int*);
}

struct task_restriction_fill
    :
    public task
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
    arg1(i1_.getVect(), BL_SPACEDIM),
    arg2(i2_)
{
    depend_on(tf = m_task_list.add_task(tf_));
}

bool
task_restriction_fill::ready ()
{
    if (is_local(m, ind))
    {
        BL_ASSERT(!tf.null());
        BL_ASSERT(tf->ready());
        task_fab* tff = dynamic_cast<task_fab*>(tf.get());
        BL_ASSERT(tff != 0);
        const Box& fb = tff->fab().box();
        const Box& pb = m[ind].box();
        (*ref)(m[ind].dataPtr(), DIMLIST(pb), DIMLIST(cbox),
	       tff->fab().dataPtr(), DIMLIST(fb),
	       D_DECL(rat[0], rat[1], rat[2]), m.nComp(),
	       &integrate, arg1.dataPtr(), arg2.dataPtr());
    }
    return true;
}

amr_restrictor::~amr_restrictor () {}

Box
amr_restrictor::box (const Box&     fb,
		     const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat);
}

Box
amr_restrictor::rebox (const Box&     cb,
		       const IntVect& rat) const
{
    Box retbox(cb);
    return retbox.refine(rat);
}

void
amr_restrictor::fill_interface (MultiFab&,
				MultiFab&,
				const level_interface&,
				const amr_boundary*,
				const IntVect&) const
{
    BoxLib::Abort("I don't think I should ever get here");
}

void
cell_average_restrictor::fill (FArrayBox&       patch,
			       const Box&       region,
			       const FArrayBox& fgr,
			       const IntVect&   rat) const
{
    BL_ASSERT(patch.box().cellCentered());
    BL_ASSERT(patch.nComp() == 1);
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()),
		 DIMLIST(region),
		 fgr.dataPtr(), DIMLIST(fgr.box()),
		 D_DECL(rat[0], rat[1], rat[2]), patch.nComp(),
		 &integrate, 0, 0);
}

void
terrain_velocity_restrictor::fill (FArrayBox&       patch,
				   const Box&       region,
				   const FArrayBox& fgr,
				   const IntVect&   rat) const
{
    BL_ASSERT(patch.box().cellCentered());
    BL_ASSERT(patch.nComp() == 1);
    const int integ = 1;
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()),
		 DIMLIST(region),
		 fgr.dataPtr(), DIMLIST(fgr.box()),
		 D_DECL(rat[0], rat[1], rat[2]), 1, &integ, 0, 0);
    Real fac = 1.0 / rat[integrate];
    patch.mult(fac, region);
}

void
injection_restrictor::fill (FArrayBox&       patch,
			    const Box&       region,
			    const FArrayBox& fgr,
			    const IntVect&   rat) const
{
    BL_ASSERT(patch.box().type() == IntVect::TheNodeVector());
    BL_ASSERT(patch.nComp() == fgr.nComp());
    FORT_FANRST1(patch.dataPtr(), DIMLIST(patch.box()),
		 DIMLIST(region),
		 fgr.dataPtr(), DIMLIST(fgr.box()),
		 D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), 0, 0, 0);
}

void
default_restrictor::fill (FArrayBox&       patch,
                          const Box&       region,
                          const FArrayBox& fgr,
                          const IntVect&   rat) const
{
    BL_ASSERT(patch.box().cellCentered()
	      || patch.box().type() == IntVect::TheNodeVector());
    BL_ASSERT(patch.nComp() == fgr.nComp());
    if (patch.box().cellCentered())
    {
        cell_average_restrictor(0).fill(patch, region, fgr, rat);
    }
    else if (patch.box().type() == IntVect::TheNodeVector())
    {
        injection_restrictor().fill(patch, region, fgr, rat);
    }
}

bilinear_restrictor::bilinear_restrictor (int  i,
					  bool hg_dense)
    :
    integrate(i),
    m_hg_dense(hg_dense)
{
    BL_ASSERT(i == 0 || i == 1);
}

Box
bilinear_restrictor::box (const Box&     fb,
			  const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat).grow(-1);
}

Box
bilinear_restrictor::rebox (const Box&     cb,
			    const IntVect& rat) const
{
    Box retbox(cb);
    return retbox.refine(rat).grow(rat-1);
}

void
bilinear_restrictor::fill (FArrayBox&       patch,
			   const Box&       region,
			   const FArrayBox& fgr,
			   const IntVect&   rat) const
{
    BL_ASSERT(patch.box().type() == IntVect::TheNodeVector());
    BL_ASSERT(patch.nComp() == fgr.nComp());
    FORT_FANRST2(patch.dataPtr(), DIMLIST(patch.box()),
                 DIMLIST(region),
                 fgr.dataPtr(), DIMLIST(fgr.box()),
                 D_DECL(rat[0], rat[1], rat[2]), patch.nComp(),
		 &integrate, 0, 0);
}

void
bilinear_restrictor::fill_interface (MultiFab&                 dest,
				     MultiFab&                 fine,
				     const level_interface&    lev_interface,
				     const amr_boundary* bdy,
				     const IntVect&            rat) const
{
    BL_ASSERT(type(dest) == IntVect::TheNodeVector());
    BL_ASSERT(dest.nComp() == fine.nComp());

    int ratmax = rat[0];
    for (int i = 1; i < BL_SPACEDIM; ++i)
        ratmax = (rat[i] > ratmax) ? rat[i] : ratmax;

    if (fine.nGrow() >= ratmax - 1)
        fill_borders(fine, lev_interface, bdy, ratmax - 1, m_hg_dense);

    const BoxArray& dest_ba = dest.boxArray();
    const BoxArray& fine_ba = fine.boxArray();

    for (int jgrid = 0; jgrid < dest.size(); jgrid++)
    {
        const Box& region = dest_ba[jgrid];
        //
        // Interface restriction is sufficiently rare and specialized that
        // we will let the restrictor handle it---at least for now.
        //
        const Box regplus = BoxLib::grow(region, 1);

	task_list tl;
        for (int iface = 0;
	     iface < lev_interface.nboxes(level_interface::FACEDIM); iface++)
        {
            if (lev_interface.flag(level_interface::FACEDIM, iface))
		continue;

            Box cbox = lev_interface.node_box(level_interface::FACEDIM, iface);
            const IntVect t =
		lev_interface.box(level_interface::FACEDIM, iface).type();
            const unsigned int geo =
		lev_interface.geo(level_interface::FACEDIM, iface);
            cbox.coarsen(rat);
            if (region.intersects(cbox))
            {
                //
                // This extends fine face by one coarse cell past coarse face:
                //
                cbox &= regplus;
                int idim = lev_interface.fdim(iface);
                cbox.grow(t - 1);
                Box fb = rebox(cbox, rat);

                if (geo == level_interface::ALL)
                {
                    //
                    // Fine grid on both sides.
                    //
                    if (fine.nGrow() >= ratmax - 1)
                    {
                        int igrid =
			    lev_interface.grid(
				level_interface::FACEDIM, iface, 0);
                        if (igrid < 0)
                            igrid =
				lev_interface.grid(
				    level_interface::FACEDIM, iface, 1);
                        task_fab* tfab =
			    new task_fab_get(tl, dest, jgrid, fb, fine, igrid);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANRST2, tl,
				dest, jgrid, cbox, tfab,
				rat, integrate));
                    }
                    else
                    {
                        task_fab* tfab =
			    new task_fill_patch(
				tl, dest, jgrid, fb, fine,
				lev_interface,
				bdy,
				level_interface::FACEDIM, iface);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANRST2, tl, dest, jgrid, cbox, tfab,
				rat, integrate));
                    }
                }
                else
                {
                    //
                    // Fine grid on just one side.
                    //
                    const int idir  = (geo & level_interface::LOW) ? -1 : 1;
                    const int igrid =
			(idir < 0) ?
			lev_interface.grid(
			    level_interface::FACEDIM, iface, 0) :
                        lev_interface.grid(
			    level_interface::FACEDIM, iface, 1) ;

                    if (igrid >= 0)
                    {
                        //
                        // Usual case, a fine grid extends all along the face.
                        //
                        task_fab* tfab =
			    new task_fab_get(tl, dest, jgrid, fb, fine, igrid);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANFR2, tl, dest, jgrid, cbox, tfab,
				rat, integrate, idim, idir));
                    }
                    else
                    {
                        //
                        // A virtual fine grid is on the other side of the boundary.
                        //
                        if (geo & level_interface::LOW)
                            fb.growHi(idim, 1 - rat[idim]);
                        else
                            fb.growLo(idim, 1 - rat[idim]);
                        task_fab* tfab =
			    new task_fill_patch(
				tl, dest, jgrid, fb, fine,
				lev_interface, bdy,
				level_interface::FACEDIM, iface);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANFR2, tl, dest, jgrid, cbox, tfab,
				rat, integrate, idim, idir));
                    }
                }
            }
        }
	tl.execute("bilinear_restrictor::fill_interface(1)");

#if (BL_SPACEDIM == 3)

        for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++)
        {
            if (lev_interface.flag(1, iedge)) continue;

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
                const Box          fb  = rebox(cbox, rat);
                const unsigned int geo = lev_interface.geo(1, iedge);

                if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1)
                {
                    int igrid = lev_interface.grid(1, iedge, 0);
                    for (int itmp = 1; igrid < 0; itmp++)
                        igrid = lev_interface.grid(1, iedge, itmp);
                    task_fab* tfab =
			new task_fab_get(tl, dest, jgrid, fb, fine, igrid);
                    tl.add_task(
			new task_restriction_fill(
			    &FORT_FANRST2, tl, dest, jgrid, cbox, tfab,
			    rat, integrate));
                }
                else
                {
                    task_fab* tfab =
			new task_fill_patch(
			    tl, dest, jgrid, fb, fine, lev_interface, bdy, 1, iedge);
                    if (geo == level_interface::ALL)
                    {
                        //
                        // Fine grid on all sides.
                        //
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANRST2, tl, dest, jgrid, cbox, tfab,
				rat, integrate));
                    }
                    else
                    {
                        Array<int> ga = lev_interface.geo_array(1, iedge);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANER2, tl, dest, jgrid, cbox, tfab,
				rat, integrate, t, ga));
                    }
                }
            }
        }
	tl.execute("bilinear_restrictor::fill_interface(2)");
#endif
        for (int icor = 0; icor < lev_interface.nboxes(0); icor++)
        {
            if (lev_interface.flag(0, icor)) continue;

            Box cbox = lev_interface.box(0, icor);
            cbox.coarsen(rat);

            if (region.intersects(cbox))
            {
                const Box          fb  = rebox(cbox, rat);
                const unsigned int geo = lev_interface.geo(0, icor);

                if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1)
                {
                    int igrid = lev_interface.grid(0, icor, 0);
                    for (int itmp = 1; igrid < 0; itmp++)
                        igrid = lev_interface.grid(0, icor, itmp);
                    task_fab* tfab =
			new task_fab_get(tl, dest, jgrid, fb, fine, igrid);
                    tl.add_task(
			new task_restriction_fill(
			    &FORT_FANRST2, tl, dest, jgrid, cbox, tfab,
			    rat, integrate));
                }
                else
                {
                    task_fab* tfab =
			new task_fill_patch(tl, dest, jgrid, fb, fine,
					    lev_interface, bdy, 0, icor);
                    if (geo == level_interface::ALL)
                    {
                        //
                        // Fine grid on all sides.
                        //
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANRST2, tl, dest, jgrid, cbox, tfab,
				rat, integrate));
                    }
                    else
                    {
                        Array<int> ga = lev_interface.geo_array(0, icor);
                        tl.add_task(
			    new task_restriction_fill(
				&FORT_FANCR2, tl, dest, jgrid, cbox, tfab,
				rat, integrate, ga));
                    }
                }
            }
        }
	tl.execute("bilinear_restrictor::fill_interface(3)");
    }
}
