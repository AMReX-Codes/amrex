//BL_COPYRIGHT_NOTICE

//
// $Id: Cluster.cpp,v 1.7 1998-01-22 16:39:53 lijewski Exp $
//

#include <Cluster.H>

enum CutStatus{ hole_cut=0, steep_cut, bisect_cut, invalid_cut };

static int findCut(const int *hist, int lo, int hi, CutStatus &status);

// ------------------------------------------------------------
Cluster::Cluster(Array<IntVect> *a) 
    : ar(a)
{
    minBox();
}

// ------------------------------------------------------------
Cluster::~Cluster()
{
    delete ar;
}

// ------------------------------------------------------------
// construct new cluster by removing all points from c
// that lie in box b.  
Cluster::Cluster(Cluster &c,const Box& b) 
    : ar(0)
{
    assert( b.ok() );
    assert( c.ar!=0 && c.ar->length() > 0 );
    if (b.contains(c.bx)) {
        bx = c.bx;
        ar = c.ar;
        c.ar = 0;
        c.bx = Box();
        return;
    }
    int len = c.ar->length();
    int *owns = new int[len];
    int nlen = 0;
    int i;
    for (i = 0; i < len; i++) {
        if (b.contains(c.ar->get(i))) {
            nlen++;
            owns[i] = 1;
        } else {
            owns[i] = 0;
        }
    }
    if (nlen == 0) {
        ar = 0;
        bx = Box();
    } else if (nlen == len) {
        bx = c.bx;
        ar = c.ar;
        c.ar = 0;
        c.bx = Box();
    } else {
        ar = new Array<IntVect>(nlen);
        Array<IntVect>* car = new Array<IntVect>(len-nlen);
        int i1 = 0;
        int i2 = 0;
        int i;
        for (i = 0; i < len; i++) {
            const IntVect &p = (c.ar->get(i));
            if (owns[i] == 1) {
                ar->set(i1++,p);
            } else {
                car->set(i2++,p);
            }
        }
        delete c.ar;
        c.ar = car;
        c.minBox();
        minBox();
    }
    delete [] owns;
}

// ------------------------------------------------------------
// construct new cluster by removing all points from c
// that lie in box b.  
void
Cluster::distribute(ClusterList &clst, const BoxDomain &bd)
{
    assert( bd.ok() );
    assert( ok() );
    assert( clst.length() == 0 );
   
    BoxDomainIterator bdi(bd);
    while (bdi && ok()) {
        Cluster *c = new Cluster(*this,bdi());
        if (c->ok()) {
            clst.append(c);
        } else {
            delete c;
        }
        ++bdi;
    }
}

// ------------------------------------------------------------
int
Cluster::numTag(const Box& b)
{
    int cnt = 0;
    int i;
    for (i = 0; i < ar->length(); i++) {
        const IntVect &p = (*ar)[i];
        if (b.contains(p)) cnt++;
    }
    return cnt;
}

// ------------------------------------------------------------
void
Cluster::minBox()
{
    int len = ar->length();
    if (len == 0) {
        bx = Box();
        return;
    }
    IntVect lo = (*ar)[0];
    IntVect hi(lo);
    int i;
    for (i = 1; i < len; i++) {
        const IntVect &p = (*ar)[i];
        lo.min(p);
        hi.max(p);
    }
    bx = Box(lo,hi);
}

// ------------------------------------------------------------
Cluster* 
Cluster::chop()
{
    int npts = ar->length();
    assert(npts > 1);

    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    const int* len = bx.length().getVect();
    //
    // Compute histogram.
    //
    int* hist[BL_SPACEDIM];
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) {
            hist[n] = new int[len[n]];
    }
    for (n = 0; n < BL_SPACEDIM; n++) {
      for (int i = 0; i < len[n]; i++) hist[n][i] = 0;
    }
    int i;
    for (i = 0; i < npts; i++) {
        const int* p = (*ar)[i].getVect();
        D_TERM( hist[0][p[0]-lo[0]]++;,
                hist[1][p[1]-lo[1]]++;,
                hist[2][p[2]-lo[2]]++; )
            }
   
      // find cutpoint ant cutstatus in each index direction
    CutStatus mincut = invalid_cut;
    CutStatus status[BL_SPACEDIM];
    int cut[BL_SPACEDIM];
    int mincount = 0;
    for (n = 0; n < BL_SPACEDIM; n++) {
        cut[n] = findCut(hist[n], lo[n], hi[n], status[n]);
        if (status[n] < mincut) {
            mincut = status[n];
            mincount = 1;
        } else if (status[n] == mincut) {
            mincount++;
        }
    }

    assert( mincut != invalid_cut );

      // select best cutpoint and direction
    int minlen = -1;
    int dir;
    for (n = 0; n < BL_SPACEDIM; n++) {
        if (status[n] == mincut) {
            int mincutlen = Min(cut[n]-lo[n],hi[n]-cut[n]);
            if (mincutlen > minlen) {
                dir = n;
                minlen = mincutlen;
            }
        }
    }

    int nlo = 0;
    for (i = lo[dir]; i < cut[dir]; i++) nlo += hist[dir][i-lo[dir]];
    assert( nlo > 0 && nlo < npts );
    int nhi = npts - nlo;

    for (i = 0; i < BL_SPACEDIM; i++) delete [] hist[i];

      // split intvect list
    Array<IntVect> *alo = new Array<IntVect>(nlo);
    Array<IntVect> *ahi = new Array<IntVect>(nhi);
    int ilo = 0;
    int ihi = 0;
    for (i = 0; i < npts; i++) {
        const IntVect &p = (*ar)[i];
        if (p[dir] < cut[dir]) {
            alo->set(ilo++,p);
        } else {
            ahi->set(ihi++,p);
        }
    }
    delete ar;
    ar = alo;
    minBox();

    Cluster* result = new Cluster(ahi);

    return result;
}


// ------------------------------------------------------------
// function FINDCUT: finds best cut location in histogram
//
const int MINOFF = 2;
const int CUT_THRESH = 2;

static int 
findCut(const int *hist, int lo, int hi, CutStatus &status)
{
    status = invalid_cut;
    int len = hi - lo + 1;

      // check validity of histogram
    if (len <= 1) return lo;

      // first find centermost point where hist == 0 (if any)
    int mid = len/2;
    int cutpoint = -1;
    int i;
    for (i = 0; i < len; i++) {
        if (hist[i] == 0) {
            status = hole_cut;
            if (abs(cutpoint-mid) > abs(i-mid)) {
                cutpoint = i;
                if (i > mid) break;
            };
        };
    };
    if (status == hole_cut) return lo+cutpoint;

      // if we got here, there was no obvious cutpoint, try
      // finding place where change in second derivative is max
    int *dhist = new int[len];
    for (i = 1; i < len-1; i++) {
        dhist[i] = hist[i+1] - 2*hist[i] + hist[i-1];
    };

    int locmax = -1;
    for(i = 0+MINOFF; i < len-MINOFF; i++) {
        int iprev = dhist[i-1];
        int icur = dhist[i];
        int locdif = abs(iprev-icur);
        if ( (iprev*icur < 0) && (locdif >= locmax) ) {
            if (locdif > locmax) {
                status = steep_cut;
                cutpoint = i;
                locmax = locdif;
            } else {
                  // select location nearest center of range
                if (abs(i-mid) < abs(cutpoint-mid)) cutpoint = i;
            };
        };
    };
    delete [] dhist;

    if (locmax <= CUT_THRESH) {
          // just recommend a bisect cut
        cutpoint = mid;
        status = bisect_cut;
    };
    return lo + cutpoint;
}

// ------------------------------------------------------------------
// ClusterList member functions
// ------------------------------------------------------------------
ClusterList::ClusterList(Array<IntVect> *pts)
{
    Cluster *c = new Cluster(pts);
    lst.append(c);
}

// ------------------------------------------------------------------
ClusterList::~ClusterList()
{
    ListIterator<Cluster*> cli(lst);
    while (cli) {
        Cluster *c = cli();
        delete c;
        ++cli;
    }
}

// ------------------------------------------------------------------
BoxArray
ClusterList::boxArray()
{
    int len = lst.length();
    BoxArray ba(len);
    ListIterator<Cluster*> cli(lst);
    int i;
    for(i = 0; i < len; i++) {
        ba.set(i,(*cli++)->box());
    }
    return ba;   
}

// ------------------------------------------------------------------
void
ClusterList::boxArray(BoxArray &ba)
{
    ba.clear();
    int len = lst.length();
    ba.resize(len);
    ListIterator<Cluster*> cli(lst);
    int i;
    for(i = 0; i < len; i++) {
        ba.set(i,(*cli++)->box());
    }
}

// ------------------------------------------------------------------
BoxList
ClusterList::boxList()
{
    BoxList blst;
    ListIterator<Cluster*> cli(lst);
    while (cli) {
        blst.append((*cli++)->box());
    }
    return blst;   
}

// ------------------------------------------------------------------
void
ClusterList::boxList(BoxList &blst)
{
    blst.clear();
    ListIterator<Cluster*> cli(lst);
    while (cli) {
        blst.append((*cli++)->box());
    }
}

// ------------------------------------------------------------------
void
ClusterList::chop(Real eff)
{
    ListIterator<Cluster*> cli(lst);
    while (cli) {
        Cluster& c = *(cli());
        if (c.eff() < eff) {
            Cluster *tmp = c.chop();
            lst.append(tmp);
        } else {
            ++cli;
        }
    }
}

// ------------------------------------------------------------------
void
ClusterList::intersect(const BoxDomain& dom)
{
    ListIterator<Cluster*> cli(lst);
    while (cli) {
        Cluster* c = cli();
        const Box& cbox = c->box();
        if (dom.contains(cbox)) {
            ++cli;
        } else {
            BoxDomain bxdom;
            ::intersect(bxdom, dom, cbox);
            if (bxdom.length() > 0) {
                ClusterList clst;
                c->distribute(clst,bxdom);
                lst.catenate(clst.lst);
            }
              // must explicitly delete c
            delete c;
            lst.remove(cli);
        }
    }
}


