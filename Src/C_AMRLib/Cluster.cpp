//
// $Id: Cluster.cpp,v 1.18 2001-08-09 16:20:37 marc Exp $
//

#include <winstd.H>

#include <algorithm>

#include <Cluster.H>
#include <BoxDomain.H>

enum CutStatus { HoleCut=0, SteepCut, BisectCut, InvalidCut };

Cluster::Cluster ()
    :
    m_ar(0) {}

Cluster::Cluster (IntVect* a, long len)
    :
    m_ar(a),
    m_len(len)
{
    minBox();
}

Cluster::~Cluster () {}

const Box&
Cluster::box () const
{
    return m_bx;
}

bool
Cluster::ok () const
{
    return m_ar != 0 && m_len > 0;
}

long
Cluster::numTag () const
{
    return m_len;
}

Real
Cluster::eff () const
{
    BL_ASSERT(ok());
    return ((Real) numTag())/((Real) m_bx.numPts());
}

//
// Predicate in call to std::partition() in Cluster::Cluster(Cluster,Box).
//
class InBox
{
public:
    InBox (const Box& b) : m_box(b) {}

    bool operator() (const IntVect& iv) const
    {
        return m_box.contains(iv);
    }
private:
    const Box& m_box;
};

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
            // None of the points in `c.m_ar' were in `b'.
            //
            m_ar  = 0;
            m_len = 0;
            m_bx  = Box();
        }
        else if (prt_it == (c.m_ar+c.m_len))
        {
            //
            // All the points in `c.m_ar' were in `b'.
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
   
    for (BoxDomain::const_iterator bdi = bd.begin();
         bdi != bd.end() && ok();
         ++bdi)
    {
        Cluster* c = new Cluster(*this, *bdi);

        if (c->ok())
            clst.append(c);
        else
            delete c;
    }
}

long
Cluster::numTag (const Box& b) const
{
    long cnt = 0;
    for (int i = 0; i < m_len; i++)
    {
        if (b.contains(m_ar[i]))
            cnt++;
    }
    return cnt;
}

void
Cluster::minBox ()
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
            if (abs(cutpoint-mid) > abs(i-mid))
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
    Array<int> dhist(len,0);
    for (i = 1; i < len-1; i++)
        dhist[i] = hist[i+1] - 2*hist[i] + hist[i-1];

    int locmax = -1;
    for (i = 0+MINOFF; i < len-MINOFF; i++)
    {
        int iprev  = dhist[i-1];
        int icur   = dhist[i];
        int locdif = abs(iprev-icur);
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
                if (abs(i-mid) < abs(cutpoint-mid))
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

//
// Predicate in call to std::partition() in Cluster::chop().
//
class Cut
{
public:
    Cut (const IntVect& cut, int dir) : m_cut(cut), m_dir(dir) {}

    bool operator() (const IntVect& iv) const
    {
        return iv[m_dir] < m_cut[m_dir] ? true : false;
    }
private:
    const IntVect& m_cut;
    int            m_dir;
};

Cluster*
Cluster::chop ()
{
    BL_ASSERT(m_len > 1);
    BL_ASSERT(!(m_ar == 0));

    const int* lo  = m_bx.loVect();
    const int* hi  = m_bx.hiVect();
    const int* len = m_bx.length().getVect();
    //
    // Compute histogram.
    //
    int* hist[BL_SPACEDIM];
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        hist[n] = new int[len[n]];
        for (int i = 0; i < len[n]; i++)
            hist[n][i] = 0;
    }
    for (int n = 0; n < m_len; n++)
    {
        const int* p = m_ar[n].getVect();
        D_TERM( hist[0][p[0]-lo[0]]++;,
                hist[1][p[1]-lo[1]]++;,
                hist[2][p[2]-lo[2]]++; )
     }
    //
    // Find cutpoint and cutstatus in each index direction.
    //
    CutStatus mincut = InvalidCut;
    CutStatus status[BL_SPACEDIM];
    IntVect cut;
    for (int n = 0, mincount = 0; n < BL_SPACEDIM; n++)
    {
        cut[n] = FindCut(hist[n], lo[n], hi[n], status[n]);
        if (status[n] < mincut)
        {
            mincut = status[n];
            mincount = 1;
        }
        else if (status[n] == mincut)
        {
            mincount++;
        }
    }
    BL_ASSERT(mincut != InvalidCut);
    //
    // Select best cutpoint and direction.
    //
    int dir;
    for (int n = 0, minlen = -1; n < BL_SPACEDIM; n++)
    {
        if (status[n] == mincut)
        {
            int mincutlen = std::min(cut[n]-lo[n],hi[n]-cut[n]);
            if (mincutlen > minlen)
            {
                dir = n;
                minlen = mincutlen;
            }
        }
    }
    BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM);

    int nlo = 0;
    for (int i = lo[dir]; i < cut[dir]; i++)
        nlo += hist[dir][i-lo[dir]];

    BL_ASSERT(nlo > 0 && nlo < m_len);

    int nhi = m_len - nlo;

    for (int i = 0; i < BL_SPACEDIM; i++)
        delete [] hist[i];

    IntVect* prt_it = std::partition(m_ar, m_ar+m_len, Cut(cut,dir));

    BL_ASSERT((prt_it-m_ar) == nlo);
    BL_ASSERT(((m_ar+m_len)-prt_it) == nhi);

    m_len = nlo;
    minBox();

    return new Cluster(prt_it, nhi);
}

ClusterList::ClusterList ()
    :
    lst()
{}

ClusterList::ClusterList (IntVect* pts,
                          long     len)
{
    lst.push_back(new Cluster(pts,len));
}

ClusterList::~ClusterList ()
{
    for (std::list<Cluster*>::iterator cli = lst.begin();
         cli != lst.end();
         ++cli)
    {
        delete *cli;
    }
}

int
ClusterList::length () const
{
    return lst.size();
}

void
ClusterList::append (Cluster* c)
{
    lst.push_back(c);
}

BoxArray
ClusterList::boxArray () const
{
    BoxArray ba(lst.size());

    int i = 0;

    for (std::list<Cluster*>::const_iterator cli = lst.begin();
         cli != lst.end();
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

    for (std::list<Cluster*>::const_iterator cli = lst.begin();
         cli != lst.end();
         ++cli, ++i)
    {
        ba.set(i,(*cli)->box());
    }
}

BoxList
ClusterList::boxList() const
{
    BoxList blst;
    for (std::list<Cluster*>::const_iterator cli = lst.begin();
         cli != lst.end();
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
    for (std::list<Cluster*>::const_iterator cli = lst.begin();
         cli != lst.end();
         ++cli)
    {
        blst.push_back((*cli)->box());
    }
}

void
ClusterList::chop (Real eff)
{
    for (std::list<Cluster*>::iterator cli = lst.begin();
         cli != lst.end(); )
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
ClusterList::intersect (const BoxDomain& dom)
{
    for (std::list<Cluster*>::iterator cli = lst.begin(); cli != lst.end(); )
    {
        Cluster* c = *cli;

        if (dom.contains(c->box()))
        {
            ++cli;
        }
        else
        {
            BoxDomain bxdom;

            BoxLib::intersect(bxdom, dom, c->box());

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
