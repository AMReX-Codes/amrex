//BL_COPYRIGHT_NOTICE

#ifndef _Cluster_H_
#define _Cluster_H_ 

//
// $Id: Cluster.H,v 1.8 1998-01-28 04:55:10 lijewski Exp $
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
#define Cluster        Cluster1d
#elif (BL_SPACEDIM==2)
#define Cluster        Cluster2d
#elif (BL_SPACEDIM==3)
#define Cluster        Cluster3d
#endif

class ClusterList;

//
// Utility class for tagging error cells.
//
class Cluster
{
public:
    //
    // Construct invalid Cluster.
    //
    Cluster () : m_ar(0) {}
    //
    // Construct a cluster from an array of IntVects.
    // Cluster object does NOT take over memory management of array.
    //
    Cluster (IntVect* a, long len);
    //
    // Construct new cluster by removing all points from c that lie in box b.
    // Cluster c is modified and may become invalid.
    //
    Cluster (Cluster&   c,
             const Box& b);
    //
    // The destructor.
    //
    ~Cluster () {}
    //
    // Return minimal box containing all tagged points.
    //
    const Box& box () const { return m_bx; }
    //
    // Does cluster contain any points?
    //
    bool ok () const { return m_ar != 0 && m_len > 0; }
    //
    //  Returns number of tagged points in cluster.
    //
    long numTag () const { return m_len; }
    //
    // Return number of tagged points in intersection of cluster and Box b.
    //
    long numTag (const Box& b) const;
    //
    // This operation splits a cluster into two pieces by selecting
    // a cutting plane with certain optimal characteristics then
    // dividing the tagged points into clusters on either side of the
    // plane.  One cluster is returned as a new object the other
    // is the modified calling object.
    //
    Cluster* chop ();
    //
    // Constructs a list of cluster objects obtained by intersecting
    // this cluster with each box in bl.  The list is retruned in the
    // argument clst.  For each intersection that includes tagged points,
    // construct a new cluster by removing those points from this
    // object.  Empty intersections or those that contain no tagged
    // points will not generate a new cluster.
    // NOTE that this cluster will be modified and possibly become
    // invalid in the process.
    //
    void distribute (ClusterList& clst, const BoxDomain& bd);
    //
    // Compute ratio of tagged to total number og points in cluster.
    //
    Real eff () const
    {
        assert(ok());
        return ((Real) numTag())/((Real) m_bx.numPts());
    }

private:
    //
    // These are disallowed.
    //
    Cluster (const Cluster&);
    Cluster& operator= (const Cluster&);
    //
    // Compute and store minimal box containing tagged points.
    //
    void minBox ();
    //
    // The data.
    //
    Box      m_bx;
    IntVect* m_ar;
    long     m_len;
};

//
// Container class for Cluster.
//

class ClusterList 
{
public:
    //
    // Constructor.
    //
    ClusterList () : lst() {}
    //
    // Constructor.
    //
    ClusterList (IntVect* pts, long len)
    {
        lst.append(new Cluster(pts,len));
    }
    //
    // Destructor.
    //
    ~ClusterList ();
    //
    // Return number of clusters in list.
    //
    int length () const { return lst.length(); }
    //
    // Add cluster to end of list.
    //
    void append (Cluster* c) { lst.append(c); }
    //
    // Return array of boxes corresponding to clusters.
    //
    BoxArray boxArray () const;
    //
    // Same but modifies argument.
    //
    void boxArray (BoxArray& ba) const;
    //
    // Return list of boxes corresponding to clusters.
    //
    BoxList boxList() const;
    //
    // Same but modifies argument.
    //
    void boxList (BoxList& blst) const;
    //
    // Chop all clusters in list that have poor efficiency.
    //
    void chop (Real eff);
    //
    // Intersect clusters with BoxDomain to insure cluster
    // boxes are interior to domain.
    //
    void intersect (const BoxDomain& dom);

private:
    //
    // These are disallowed.
    //
    ClusterList (const ClusterList&);
    ClusterList& operator= (const ClusterList&);
    //
    // The data.
    //    
    List<Cluster*> lst;
};

#endif
