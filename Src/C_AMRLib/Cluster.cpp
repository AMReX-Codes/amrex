
//
// $Id: Cluster.cpp,v 1.1 1997-11-18 19:30:22 lijewski Exp $
//

#include <Cluster.H>

enum CutStatus{ hole_cut=0, steep_cut, bisect_cut, invalid_cut };

static int findCut(const int *hist, int lo, int hi, CutStatus &status);

// ------------------------------------------------------------
CLUSTER::CLUSTER(Array<INTVECT> *a) 
    : ar(a)
{
    minBox();
}

// ------------------------------------------------------------
CLUSTER::~CLUSTER()
{
    delete ar;
}

// ------------------------------------------------------------
// construct new cluster by removing all points from c
// that lie in box b.  
CLUSTER::CLUSTER(CLUSTER &c,const BOX& b) 
    : ar(0)
{
    assert( b.ok() );
    assert( c.ar!=0 && c.ar->length() > 0 );
    if (b.contains(c.bx)) {
	bx = c.bx;
	ar = c.ar;
	c.ar = 0;
	c.bx = BOX();
	return;
    }
    int len = c.ar->length();
    int *owns = new int[len];
    assert( owns != 0 );
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
	bx = BOX();
    } else if (nlen == len) {
	bx = c.bx;
	ar = c.ar;
	c.ar = 0;
	c.bx = BOX();
    } else {
	ar = new Array<INTVECT>(nlen);
	Array<INTVECT> *car = new Array<INTVECT>(len-nlen);
	int i1 = 0;
	int i2 = 0;
        int i;
	for (i = 0; i < len; i++) {
	    const INTVECT &p = (c.ar->get(i));
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
    delete owns;
}

// ------------------------------------------------------------
// construct new cluster by removing all points from c
// that lie in box b.  
void
CLUSTER::distribute(ClusterList &clst, const BoxDomain &bd)
{
    assert( bd.ok() );
    assert( ok() );
    assert( clst.length() == 0 );
   
    BoxDomainIterator bdi(bd);
    while (bdi && ok()) {
	CLUSTER *c = new CLUSTER(*this,bdi());
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
CLUSTER::numTag(const BOX& b)
{
    int cnt = 0;
    int i;
    for (i = 0; i < ar->length(); i++) {
	const INTVECT &p = (*ar)[i];
	if (b.contains(p)) cnt++;
    }
    return cnt;
}

// ------------------------------------------------------------
void
CLUSTER::minBox()
{
    int len = ar->length();
    if (len == 0) {
	bx = BOX();
	return;
    }
    INTVECT lo = (*ar)[0];
    INTVECT hi(lo);
    int i;
    for (i = 1; i < len; i++) {
	const INTVECT &p = (*ar)[i];
	lo.min(p);
	hi.max(p);
    }
    bx = BOX(lo,hi);
}

// ------------------------------------------------------------
CLUSTER* 
CLUSTER::chop()
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
    D_TERM( hist[0] = new int[len[0]]; ,
	    hist[1] = new int[len[1]]; ,
	    hist[2] = new int[len[2]]; )
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) {
      int i;
      for (i = 0; i < len[n]; i++) hist[n][i] = 0;
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

    for (i = 0; i < BL_SPACEDIM; i++) delete hist[i];

      // split intvect list
    Array<INTVECT> *alo = new Array<INTVECT>(nlo);
    Array<INTVECT> *ahi = new Array<INTVECT>(nhi);
    int ilo = 0;
    int ihi = 0;
    for (i = 0; i < npts; i++) {
	const INTVECT &p = (*ar)[i];
	if (p[dir] < cut[dir]) {
	    alo->set(ilo++,p);
	} else {
	    ahi->set(ihi++,p);
	}
    }
    delete ar;
    ar = alo;
    minBox();

    return new CLUSTER(ahi);
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
    assert( dhist );
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
    delete dhist;

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
ClusterList::ClusterList(Array<INTVECT> *pts)
{
    CLUSTER *c = new CLUSTER(pts);
    lst.append(c);
}

// ------------------------------------------------------------------
ClusterList::~ClusterList()
{
    ListIterator<CLUSTER*> cli(lst);
    while (cli) {
	CLUSTER *c = cli();
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
    ListIterator<CLUSTER*> cli(lst);
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
    ListIterator<CLUSTER*> cli(lst);
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
    ListIterator<CLUSTER*> cli(lst);
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
    ListIterator<CLUSTER*> cli(lst);
    while (cli) {
	blst.append((*cli++)->box());
    }
}

// ------------------------------------------------------------------
void
ClusterList::chop(REAL eff)
{
    ListIterator<CLUSTER*> cli(lst);
    while (cli) {
	CLUSTER& c = *(cli());
	if (c.eff() < eff) {
	    CLUSTER *tmp = c.chop();
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
    ListIterator<CLUSTER*> cli(lst);
    while (cli) {
	CLUSTER* c = cli();
	const BOX& cbox = c->box();
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


