//BL_COPYRIGHT_NOTICE

#ifndef _Cluster_H_
#define _Cluster_H_ 

//
// $Id: Cluster.H,v 1.3 1997-11-26 19:18:48 lijewski Exp $
//

#include <IntVect.H>
#include <Box.H>
#include <Array.H>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <REAL.H>

#if !defined(BL_SPACEDIM) || (BL_SPACEDIM != 1 && BL_SPACEDIM !=2 && BL_SPACEDIM != 3)
#error BL_SPACEDIM must be defined
#endif

#if (BL_SPACEDIM==1)
#define CLUSTER	Cluster1d
#elif (BL_SPACEDIM==2)
#define CLUSTER	Cluster2d
#elif (BL_SPACEDIM==3)
#define CLUSTER	Cluster3d
#endif

class ClusterList;

/// Utility class for tagging error cells
class CLUSTER
{

/*  This class is basically a fortran array of INTEGERSs.
The array is defined to exist for a set of points inside a box.
*/

protected:
    Box bx;
    Array<IntVect> *ar;

public:
/// administrative functions

/// construct invalid Cluster
    CLUSTER() : ar(0) {}

//  construct a cluster from an array of IntVectS
//  NOTE: Cluster object takes over memory management of array
    CLUSTER(Array<IntVect> *a);

//  construct new cluster by removing all points from c
//  that lie in box b.  Cluster c is modified and may become invalid.
    CLUSTER(CLUSTER &c,const Box& b);

//  destructor
    ~CLUSTER();

//  return minimal box containing all tagged points
    const Box& box() const {return bx;}

//  return ith tagged point in cluster
    IntVect operator[](int i) {
	assert(ar != 0);
	return (*ar)[i];
    }

//  does cluster contain any points?
    int ok() const { return (ar!=0 && ar->length() > 0); }

//  returns number of tagged points in cluster
    int numTag() const {return ar->length();}

//  return number of tagged points in intersection of cluster and Box b.
    int numTag(const Box& b);

/*  This operation splits a cluster into two pieces by selecting
    a cutting plane with certain optimal characteristics then
    dividing the tagged points into clusters on either side of the
    plane.  One cluster is returned as a new object the other
    is the modified calling object. */

    CLUSTER* chop();

/*  constructs a list of cluster objects obtained by intersecting
    this cluster with each box in bl.  The list is retruned in the
    argument clst.  For each intersection that includes tagged points,
    construct a new cluster by removing those points from this
    object.  Empty intersections or those that contain no tagged
    points will not generate a new cluster.
    NOTE that this cluster will be modified and possibly become
    invalid in the process. */

    void distribute(ClusterList &clst, const BoxDomain &bd);

//  compute ratio of tagged to total number og points in cluster
    Real eff() const {
	assert( ok() );
	return ((Real) numTag())/((Real) bx.numPts());
    }

protected:
//  compute and store minimal box containing tagged points
    void minBox();

private:
    CLUSTER(const CLUSTER&);
    CLUSTER& operator=(const CLUSTER&);

};

// --------------------------------------------------------------
/// Container class for CLUSTER
class ClusterList 
{

protected:
    List<CLUSTER*> lst;
    
public:
//  constructor
    ClusterList() : lst() {};

//  constructor
    ClusterList(Array<IntVect> *pts);

//  destructor
    ~ClusterList();

//  return number of clusters in list
    int length() const { return lst.length(); };

//  add cluster to end of list
    void append(CLUSTER* c) { lst.append(c); };

//  return array of boxes corresponding to clusters
    BoxArray boxArray();

//  same but modifies argument
    void boxArray(BoxArray &ba);

//  return list of boxes corresponding to clusters
    BoxList  boxList();

//  same but modifies argument
    void boxList(BoxList &blst);

//  chop all clusters in list that have poor efficiency
    void chop(Real eff);

//  intersect clusters with BoxDomain to insure cluster
//  boxes are interior to domain
    void intersect(const BoxDomain& dom);

};

#endif
