//BL_COPYRIGHT_NOTICE

//
// $Id: Cluster.cpp,v 1.8 1998-01-26 22:00:01 lijewski Exp $
//

#include <Cluster.H>

enum CutStatus { HoleCut=0, SteepCut, BisectCut, InvalidCut };

Cluster::Cluster (Array<IntVect>* a) 
    :
    m_ar(a)
{
    minBox();
}

Cluster::~Cluster ()
{
    delete m_ar;
}

//
// Construct new cluster by removing all points from c that lie in box b.
//

Cluster::Cluster (Cluster&   c,
                  const Box& b) 
    :
    m_ar(0)
{
    assert(b.ok());
    assert(c.m_ar!=0 && c.m_ar->length() > 0);

    if (b.contains(c.m_bx))
    {
        m_bx = c.m_bx;
        m_ar = c.m_ar;
        c.m_ar = 0;
        c.m_bx = Box();
    }
    else
    {
        int len = c.m_ar->length();
        int* owns = new int[len];
        int nlen = 0;
        for (int i = 0; i < len; i++)
        {
            if (b.contains(c.m_ar->get(i)))
            {
                nlen++;
                owns[i] = 1;
            }
            else
            {
                owns[i] = 0;
            }
        }
        if (nlen == 0)
        {
            m_ar = 0;
            m_bx = Box();
        }
        else if (nlen == len)
        {
            m_bx = c.m_bx;
            m_ar = c.m_ar;
            c.m_ar = 0;
            c.m_bx = Box();
        }
        else
        {
            m_ar = new Array<IntVect>(nlen);
            Array<IntVect>* cm_ar = new Array<IntVect>(len-nlen);
            int i1 = 0;
            int i2 = 0;
            for (int i = 0; i < len; i++)
            {
                if (owns[i] == 1)
                {
                    m_ar->set(i1++, c.m_ar->get(i));
                }
                else
                {
                    cm_ar->set(i2++, c.m_ar->get(i));
                }
            }
            delete c.m_ar;
            c.m_ar = cm_ar;
            c.minBox();
            minBox();
        }
        delete [] owns;
    }
}

//
// Construct new cluster by removing all points from c that lie in box b.
//

void
Cluster::distribute (ClusterList&     clst,
                     const BoxDomain& bd)
{
    assert(bd.ok());
    assert(ok());
    assert(clst.length() == 0);
   
    for (BoxDomainIterator bdi(bd); bdi && ok(); ++bdi)
    {
        Cluster* c = new Cluster(*this, bdi());

        if (c->ok())
        {
            clst.append(c);
        }
        else
        {
            delete c;
        }
    }
}

long
Cluster::numTag (const Box& b) const
{
    long cnt = 0;
    for (int i = 0; i < m_ar->length(); i++)
    {
        if (b.contains((*m_ar)[i]))
            cnt++;
    }
    return cnt;
}

void
Cluster::minBox ()
{
    long len = m_ar->length();

    if (len == 0)
    {
        m_bx = Box();
    }
    else
    {
        IntVect lo = (*m_ar)[0];
        IntVect hi(lo);
        for (int i = 1; i < len; i++)
        {
            lo.min((*m_ar)[i]);
            hi.max((*m_ar)[i]);
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
    if (len <= 1) return lo;
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
        return lo+cutpoint;
    //
    // If we got here, there was no obvious cutpoint, try
    // finding place where change in second derivative is max.
    //
    int* dhist = new int[len];
    for (i = 1; i < len-1; i++)
        dhist[i] = hist[i+1] - 2*hist[i] + hist[i-1];

    int locmax = -1;
    for (i = 0+MINOFF; i < len-MINOFF; i++)
    {
        int iprev = dhist[i-1];
        int icur = dhist[i];
        int locdif = abs(iprev-icur);
        if ((iprev*icur < 0) && (locdif >= locmax))
        {
            if (locdif > locmax)
            {
                status = SteepCut;
                cutpoint = i;
                locmax = locdif;
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
    delete [] dhist;

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

Cluster* 
Cluster::chop ()
{
    int npts = m_ar->length();
    assert(npts > 1);

    const int* lo  = m_bx.loVect();
    const int* hi  = m_bx.hiVect();
    const int* len = m_bx.length().getVect();
    //
    // Compute histogram.
    //
    int* hist[BL_SPACEDIM];
    int n;
    for (n = 0; n < BL_SPACEDIM; n++)
    {
        hist[n] = new int[len[n]];
    }
    for (n = 0; n < BL_SPACEDIM; n++)
    {
        for (int i = 0; i < len[n]; i++)
            hist[n][i] = 0;
    }
    int i;
    for (i = 0; i < npts; i++)
    {
        const int* p = (*m_ar)[i].getVect();
        D_TERM( hist[0][p[0]-lo[0]]++;,
                hist[1][p[1]-lo[1]]++;,
                hist[2][p[2]-lo[2]]++; )
     }
    //
    // Find cutpoint and cutstatus in each index direction.
    //
    CutStatus mincut = InvalidCut;
    CutStatus status[BL_SPACEDIM];
    int cut[BL_SPACEDIM];
    int mincount = 0;
    for (n = 0; n < BL_SPACEDIM; n++)
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
    assert(mincut != InvalidCut);
    //
    // Select best cutpoint and direction.
    //
    int minlen = -1;
    int dir;
    for (n = 0; n < BL_SPACEDIM; n++)
    {
        if (status[n] == mincut)
        {
            int mincutlen = Min(cut[n]-lo[n],hi[n]-cut[n]);
            if (mincutlen > minlen)
            {
                dir = n;
                minlen = mincutlen;
            }
        }
    }

    int nlo = 0;
    for (i = lo[dir]; i < cut[dir]; i++)
        nlo += hist[dir][i-lo[dir]];

    assert(nlo > 0 && nlo < npts);

    int nhi = npts - nlo;

    for (i = 0; i < BL_SPACEDIM; i++)
        delete [] hist[i];
    //
    // Split intvect list.
    //
    Array<IntVect>* alo = new Array<IntVect>(nlo);
    Array<IntVect>* ahi = new Array<IntVect>(nhi);
    int ilo = 0;
    int ihi = 0;
    for (i = 0; i < npts; i++)
    {
        const IntVect& p = (*m_ar)[i];
        if (p[dir] < cut[dir])
        {
            alo->set(ilo++, p);
        }
        else
        {
            ahi->set(ihi++, p);
        }
    }
    delete m_ar;
    m_ar = alo;
    minBox();

    Cluster* result = new Cluster(ahi);

    return result;
}

ClusterList::~ClusterList ()
{
    for (ListIterator<Cluster*> cli(lst); cli; ++cli)
    {
        delete cli();
    }
}

BoxArray
ClusterList::boxArray () const
{
    long len = lst.length();
    BoxArray ba(len);
    ListIterator<Cluster*> cli(lst);
    for (int i = 0; i < len; i++)
        ba.set(i,(*cli++)->box());
    return ba;   
}

void
ClusterList::boxArray (BoxArray &ba) const
{
    ba.clear();
    long len = lst.length();
    ba.resize(len);
    ListIterator<Cluster*> cli(lst);
    for (int i = 0; i < len; i++)
        ba.set(i,(*cli++)->box());
}

BoxList
ClusterList::boxList() const
{
    BoxList blst;
    for (ListIterator<Cluster*> cli(lst); cli; ++cli)
    {
        blst.append((*cli)->box());
    }
    return blst;   
}

void
ClusterList::boxList (BoxList& blst) const
{
    blst.clear();
    for (ListIterator<Cluster*> cli(lst); cli; ++cli)
    {
        blst.append((*cli)->box());
    }
}

void
ClusterList::chop (Real eff)
{
    for (ListIterator<Cluster*> cli(lst); cli; )
    {
        if (cli()->eff() < eff)
        {
            lst.append(cli()->chop());
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
    for (ListIterator<Cluster*> cli(lst); cli; )
    {
        Cluster* c = cli();
        const Box& cbox = c->box();
        if (dom.contains(cbox))
        {
            ++cli;
        }
        else
        {
            BoxDomain bxdom;
            ::intersect(bxdom, dom, cbox);
            if (bxdom.length() > 0)
            {
                ClusterList clst;
                c->distribute(clst,bxdom);
                lst.catenate(clst.lst);
            }
            //
            // Must explicitly delete c.
            //
            delete c;
            lst.remove(cli);
        }
    }
}
