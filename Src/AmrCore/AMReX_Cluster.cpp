
#include <algorithm>
#include <cmath>
#include <AMReX_Cluster.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_BLProfiler.H>

namespace amrex {

namespace {
enum CutStatus { HoleCut=0, SteepCut, BisectCut, InvalidCut };
}

Cluster::Cluster () noexcept
{}

Cluster::Cluster (IntVect* a, Long len) noexcept
    :
    m_ar(a),
    m_len(len)
{
    minBox();
}

Cluster::~Cluster () {}

namespace {
//
// Predicate in call to std::partition() in Cluster::Cluster(Cluster,Box).
//
class InBox
{
public:
    explicit InBox (const Box& b) noexcept : m_box(b) {}

    // You might see compiler warning on this is never referenced.
    // The compiler is wrong.
    bool operator() (const IntVect& iv) const noexcept
    {
        return m_box.contains(iv);
    }
private:
    Box m_box;
};
}

Cluster::Cluster (Cluster&   c,
                  const Box& b) 
    :
    m_ar(0),
    m_len(0)
{
    BL_ASSERT(b.ok());
    BL_ASSERT(c.m_ar != 0 && c.m_len > 0);

    if (b.contains(c.m_bx))
    {
        m_bx    = c.m_bx;
        m_ar    = c.m_ar;
        m_len   = c.m_len;
        c.m_ar  = 0;
        c.m_len = 0;
        c.m_bx  = Box();
    }
    else
    {
        IntVect* prt_it = std::partition(c.m_ar, c.m_ar+c.m_len, InBox(b));

        if (prt_it == c.m_ar)
        {
            //
            // None of the points in c.m_ar were in b.
            //
            m_ar  = 0;
            m_len = 0;
            m_bx  = Box();
        }
        else if (prt_it == (c.m_ar+c.m_len))
        {
            //
            // All the points in c.m_ar were in b.
            //
            m_bx    = c.m_bx;
            m_ar    = c.m_ar;
            m_len   = c.m_len;
            c.m_ar  = 0;
            c.m_len = 0;
            c.m_bx  = Box();
        }
        else
        {
            m_ar    = c.m_ar;
            m_len   = prt_it - m_ar;
            c.m_ar  = prt_it;
            c.m_len = c.m_len - m_len;
            minBox();
            c.minBox();
        }
    }
}

void
Cluster::distribute (ClusterList&     clst,
                     const BoxDomain& bd)
{
    BL_ASSERT(ok());
    BL_ASSERT(bd.ok());
    BL_ASSERT(clst.length() == 0);
   
    for (BoxDomain::const_iterator bdi = bd.begin(), End = bd.end();
         bdi != End && ok();
         ++bdi)
    {
        Cluster* c = new Cluster(*this, *bdi);

        if (c->ok())
            clst.append(c);
        else
            delete c;
    }
}

Long
Cluster::numTag (const Box& b) const noexcept
{
    Long cnt = 0;
    for (int i = 0; i < m_len; i++)
    {
        if (b.contains(m_ar[i]))
            cnt++;
    }
    return cnt;
}

void
Cluster::minBox () noexcept
{
    if (m_len == 0)
    {
        m_bx = Box();
    }
    else
    {
        IntVect lo = m_ar[0], hi = lo;
        for (int i = 1; i < m_len; i++)
        {
            lo.min(m_ar[i]);
            hi.max(m_ar[i]);
        }
        m_bx = Box(lo,hi);
    }
}

//
// Finds best cut location in histogram.
//

static
int 
FindCut (const int* hist,
         int        lo,
         int        hi,
         CutStatus& status)
{
    const int MINOFF     = 2;
    const int CUT_THRESH = 2;

    status = InvalidCut;
    int len = hi - lo + 1;
    //
    // Check validity of histogram.
    //
    if (len <= 1)
        return lo;
    //
    // First find centermost point where hist == 0 (if any).
    //
    int mid = len/2;
    int cutpoint = -1;
    int i;
    for (i = 0; i < len; i++)
    {
        if (hist[i] == 0)
        {
            status = HoleCut;
            if (std::abs(cutpoint-mid) > std::abs(i-mid))
            {
                cutpoint = i;
                if (i > mid)
                    break;
            }
        }
    }
    if (status == HoleCut)
        return lo + cutpoint;
    //
    // If we got here, there was no obvious cutpoint, try
    // finding place where change in second derivative is max.
    //
    Vector<int> dhist(len,0);
    for (i = 1; i < len-1; i++)
        dhist[i] = hist[i+1] - 2*hist[i] + hist[i-1];

    int locmax = -1;
    for (i = 0+MINOFF; i < len-MINOFF; i++)
    {
        int iprev  = dhist[i-1];
        int icur   = dhist[i];
        int locdif = std::abs(iprev-icur);
        if (iprev*icur < 0 && locdif >= locmax)
        {
            if (locdif > locmax)
            {
                status   = SteepCut;
                cutpoint = i;
                locmax   = locdif;
            }
            else
            {
                //
                // Select location nearest center of range.
                //
                if (std::abs(i-mid) < std::abs(cutpoint-mid))
                    cutpoint = i;
            }
        }
    }

    if (locmax <= CUT_THRESH)
    {
        //
        // Just recommend a bisect cut.
        //
        cutpoint = mid;
        status = BisectCut;
    }

    return lo + cutpoint;
}

namespace {
//
// Predicate in call to std::partition() in Cluster::chop().
//
class Cut
{
public:
    Cut (const IntVect& cut, int dir) : m_cut(cut), m_dir(dir) {}

    // You might see compiler warning on this is never referenced.
    // The compiler is wrong.
    bool operator() (const IntVect& iv) const
    {
        return iv[m_dir] < m_cut[m_dir];
    }
private:
    IntVect m_cut;
    int     m_dir;
};
}

Cluster*
Cluster::chop ()
{
    BL_ASSERT(m_len > 1);
    BL_ASSERT(!(m_ar == 0));

    const int*    lo  = m_bx.loVect();
    const int*    hi  = m_bx.hiVect();
    const IntVect len = m_bx.size();
    //
    // Compute histogram.
    //
    Array<Vector<int>,AMREX_SPACEDIM> hist;
    AMREX_D_TERM(hist[0].resize(len[0], 0);,
                 hist[1].resize(len[1], 0);,
                 hist[2].resize(len[2], 0);)
    for (int n = 0; n < m_len; n++)
    {
        const int* p = m_ar[n].getVect();
        AMREX_D_TERM( hist[0][p[0]-lo[0]]++;,
                hist[1][p[1]-lo[1]]++;,
                hist[2][p[2]-lo[2]]++; )
     }
    //
    // Find cutpoint and cutstatus in each index direction.
    //
    CutStatus mincut = InvalidCut;
    CutStatus status[AMREX_SPACEDIM];
    IntVect cut;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        cut[n] = FindCut(hist[n].data(), lo[n], hi[n], status[n]);
        if (status[n] < mincut)
        {
            mincut = status[n];
        }
    }
    BL_ASSERT(mincut != InvalidCut);
    //
    // Select best cutpoint and direction.
    //
    int dir = -1;
    for (int n = 0, minlen = -1; n < AMREX_SPACEDIM; n++)
    {
        if (status[n] == mincut)
        {
            int mincutlen = std::min(cut[n]-lo[n],hi[n]-cut[n]);
            if (mincutlen >= minlen)
            {
                dir = n;
                minlen = mincutlen;
            }
        }
    }
    BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM);

    int nlo = 0;
    for (int i = lo[dir]; i < cut[dir]; i++) {
        nlo += hist[dir][i-lo[dir]];
    }

    BL_ASSERT(nlo > 0 && nlo < m_len);

    int nhi = m_len - nlo;

    IntVect* prt_it = std::partition(m_ar, m_ar+m_len, Cut(cut,dir));

    BL_ASSERT((prt_it-m_ar) == nlo);
    BL_ASSERT(((m_ar+m_len)-prt_it) == nhi);

    m_len = nlo;
    minBox();

    return new Cluster(prt_it, nhi);
}

Cluster*
Cluster::new_chop ()
{
    BL_ASSERT(m_len > 1);
    BL_ASSERT(!(m_ar == 0));

    const int*    lo  = m_bx.loVect();
    const int*    hi  = m_bx.hiVect();
    const IntVect len = m_bx.size();
    //
    // Compute histogram.
    //
    Array<Vector<int>,AMREX_SPACEDIM> hist;
    AMREX_D_TERM(hist[0].resize(len[0], 0);,
                 hist[1].resize(len[1], 0);,
                 hist[2].resize(len[2], 0);)
    for (int n = 0; n < m_len; n++)
    {
        const int* p = m_ar[n].getVect();
        AMREX_D_TERM( hist[0][p[0]-lo[0]]++;,
                      hist[1][p[1]-lo[1]]++;,
                      hist[2][p[2]-lo[2]]++; )
    }

    int invalid_dir = -1;
    for (int n_try = 0; n_try < 2; n_try++)
    {
       //
       // Find cutpoint and cutstatus in each index direction.
       //
       CutStatus mincut = InvalidCut;
       CutStatus status[AMREX_SPACEDIM];
       IntVect cut;
       for (int n = 0; n < AMREX_SPACEDIM; n++)
       {
           if (n != invalid_dir)
           {
              cut[n] = FindCut(hist[n].data(), lo[n], hi[n], status[n]);
              if (status[n] < mincut)
              {
                  mincut = status[n];
              }
           }
       }
       BL_ASSERT(mincut != InvalidCut);
       //
       // Select best cutpoint and direction.
       //
       int dir = -1;
       for (int n = 0, minlen = -1; n < AMREX_SPACEDIM; n++)
       {
           if (status[n] == mincut)
           {
               int mincutlen = std::min(cut[n]-lo[n],hi[n]-cut[n]);
               if (mincutlen >= minlen)
               {
                   dir = n;
                   minlen = mincutlen;
               }
           }
       }
       BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM);
   
       int nlo = 0;
       for (int i = lo[dir]; i < cut[dir]; i++) {
           nlo += hist[dir][i-lo[dir]];
       }

       if (nlo <= 0 || nlo >= m_len) return chop();

       int nhi = m_len - nlo;

       IntVect* prt_it = std::partition(m_ar, m_ar+m_len, Cut(cut,dir));

       BL_ASSERT((prt_it-m_ar) == nlo);
       BL_ASSERT(((m_ar+m_len)-prt_it) == nhi);

       // These refer to the box that was originally passed in
       Real oldeff = eff();
       int orig_mlen = m_len;

       // Define the new box "above" the cut
       std::unique_ptr<Cluster> newbox(new Cluster(prt_it, nhi));
       Real neweff = newbox->eff();

       // Replace the current box by the part of the box "below" the cut
       m_len = nlo;
       minBox();
   
       if ( (eff() > oldeff) || (neweff > oldeff) || n_try > 0)
       {
          return newbox.release();

       } else {

          // Restore the original box and try again, cutting in a different direction
          m_len = orig_mlen;
          minBox();
          invalid_dir = dir;
       }
   }

    amrex::Abort("Should never reach this point in Cluster::new_chop()");
    return this;
}

ClusterList::ClusterList ()
    :
    lst()
{}

ClusterList::ClusterList (IntVect* pts, Long len)
{
    lst.push_back(new Cluster(pts,len));
}

ClusterList::~ClusterList ()
{
    for (std::list<Cluster*>::iterator cli = lst.begin(), End = lst.end();
         cli != End;
         ++cli)
    {
        delete *cli;
    }
}

BoxArray
ClusterList::boxArray () const
{
    BoxArray ba(lst.size());

    int i = 0;

    for (std::list<Cluster*>::const_iterator cli = lst.begin(), End = lst.end();
         cli != End;
         ++cli, ++i)
    {
        ba.set(i,(*cli)->box());
    }

    return ba;   
}

void
ClusterList::boxArray (BoxArray& ba) const
{
    ba.clear();

    ba.resize(lst.size());

    int i = 0;

    for (std::list<Cluster*>::const_iterator cli = lst.begin(), End = lst.end();
         cli != End;
         ++cli, ++i)
    {
        ba.set(i,(*cli)->box());
    }
}

BoxList
ClusterList::boxList() const
{
    BoxList blst;
    blst.reserve(lst.size());
    for (std::list<Cluster*>::const_iterator cli = lst.begin(), End = lst.end();
         cli != End;
         ++cli)
    {
        blst.push_back((*cli)->box());
    }
    return blst;   
}

void
ClusterList::boxList (BoxList& blst) const
{
    blst.clear();
    blst.reserve(lst.size());
    for (std::list<Cluster*>::const_iterator cli = lst.begin(), End = lst.end();
         cli != End;
         ++cli)
    {
        blst.push_back((*cli)->box());
    }
}

void
ClusterList::chop (Real eff)
{
    BL_PROFILE("ClusterList::chop()");

    for (std::list<Cluster*>::iterator cli = lst.begin(); cli != lst.end(); )
    {
        if ((*cli)->eff() < eff)
        {
            lst.push_back((*cli)->chop());
        }
        else
        {
            ++cli;
        }
    }
}

void
ClusterList::new_chop (Real eff)
{
    BL_PROFILE("ClusterList::new_chop()");

    for (std::list<Cluster*>::iterator cli = lst.begin(); cli != lst.end(); )
    {
        if ((*cli)->eff() < eff)
        {
            lst.push_back((*cli)->new_chop());
        }
        else
        {
            ++cli;
        }
    }
}

void
ClusterList::intersect (const BoxDomain& dom)
{
    BL_PROFILE("ClusterList::intersect()");

    //
    // Make a BoxArray covering dom.
    // We'll use this to speed up the contains() test below.
    //
    BoxArray domba(dom.boxList());

    for (std::list<Cluster*>::iterator cli = lst.begin(); cli != lst.end(); )
    {
        Cluster* c = *cli;

	bool assume_disjoint_ba = true;
        if (domba.contains(c->box(),assume_disjoint_ba))
        {
            ++cli;
        }
        else
        {
            BoxDomain bxdom;

            amrex::intersect(bxdom, dom, c->box());

            if (bxdom.size() > 0)
            {
                ClusterList clst;
                c->distribute(clst,bxdom);
                lst.splice(lst.end(),clst.lst);
            }
            //
            // Must explicitly delete c.
            //
            delete c;

            lst.erase(cli++);
        }
    }
}

}
