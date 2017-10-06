
#ifndef _Cluster_H_
#define _Cluster_H_ 

#include <list>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>

namespace amrex {

class BoxDomain;
class ClusterList;

//
// A cluster of tagged cells.
//
// Utility class for tagging error cells.
//

class Cluster
{
public:
    //
    // The default constructor -- builds invalid Cluster.
    //
    Cluster ();
    //
    // Construct a cluster from an array of IntVects.
    // The Cluster object does NOT take over memory management
    // of array; i.e. that's the user's ultimate responsibility.
    //
    Cluster (IntVect* a,
             long     len);
    //
    // Construct new cluster by removing all points from c that lie
    // in box b.  Cluster c is modified and may become invalid.
    //
    Cluster (Cluster&   c,
             const Box& b);
    //
    // The destructor.  Does NOT destroy the initial IntVect array.
    //
    ~Cluster ();
    //
    // Return minimal box containing all tagged points.
    //
    const Box& box () const { return m_bx; }
    //
    // Does cluster contain any points?
    //
    bool ok () const { return m_ar != 0 && m_len > 0; }
    //
    // Returns number of tagged points in cluster.
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
    // argument clst.  For each intersection that includes tagged
    // points, construct a new cluster by removing those points from
    // this object.  Empty intersections or those that contain no
    // tagged points will not generate a new cluster.
    // Note that this cluster will be modified and possibly become
    // invalid in the process.
    //
    void distribute (ClusterList&     clst,
                     const BoxDomain& bd);
    //
    // Compute ratio of tagged to total number of points in cluster.
    //
    Real eff () const { BL_ASSERT(ok()); return numTag()/m_bx.d_numPts(); }

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
// A list of Cluster objects
//
// A container class for Cluster.
//

class ClusterList 
{
public:
    //
    // The default constructor.
    //
    ClusterList ();
    //
    // Construct a list containing Cluster(pts,len).
    //
    ClusterList (IntVect* pts,
                 long     len);
    //
    // The destructor.
    //
    ~ClusterList ();
    //
    // Return number of clusters in list.
    //
    int length () const { return lst.size(); }
    //
    // Add cluster to end of list.
    //
    void append (Cluster* c) { lst.push_back(c); }
    //
    // Return array of boxes corresponding to clusters.
    //
    BoxArray boxArray () const;
    //
    // Return array of boxes corresponding to clusters in argument.
    //
    void boxArray (BoxArray& ba) const;
    //
    // Return list of boxes corresponding to clusters.
    //
    BoxList boxList() const;
    //
    // Return list of boxes corresponding to clusters in argument.
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
    std::list<Cluster*> lst;
};

}

#endif /*_Cluster_H_*/
