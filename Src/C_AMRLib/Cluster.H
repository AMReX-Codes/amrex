
#ifndef _Cluster_H_
#define _Cluster_H_ 

//
// $Id: Cluster.H,v 1.12 2000-10-02 20:48:41 lijewski Exp $
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
//@Man:
//@Memo: A cluster of tagged cells.
/*@Doc:

  Utility class for tagging error cells.
*/

class Cluster
{
public:
    //
    //@ManDoc: The default constructor -- builds invalid Cluster.
    //
    Cluster ()
        :
        m_ar(0) {}

    /*@ManDoc: Construct a cluster from an array of IntVects.
               The Cluster object does NOT take over memory management
               of array; i.e. that's the user's ultimate responsibility.
    */
    Cluster (IntVect* a,
             long     len);

    /*@ManDoc: Construct new cluster by removing all points from c that lie
               in box b.  Cluster c is modified and may become invalid.
    */
    Cluster (Cluster&   c,
             const Box& b);
    //
    //@ManDoc: The destructor.  Does NOT destroy the initial IntVect array.
    //
    ~Cluster () {}
    //
    //@ManDoc: Return minimal box containing all tagged points.
    //
    const Box& box () const { return m_bx; }
    //
    //@ManDoc: Does cluster contain any points?
    //
    bool ok () const { return m_ar != 0 && m_len > 0; }
    //
    //@ManDoc: Returns number of tagged points in cluster.
    //
    long numTag () const { return m_len; }
    //
    //@ManDoc: Return number of tagged points in intersection of cluster and Box b.
    //
    long numTag (const Box& b) const;

    /*@ManDoc: This operation splits a cluster into two pieces by selecting
                a cutting plane with certain optimal characteristics then
                dividing the tagged points into clusters on either side of the
                plane.  One cluster is returned as a new object the other
                is the modified calling object.
    */
    Cluster* chop ();

    /*@ManDoc: Constructs a list of cluster objects obtained by intersecting
               this cluster with each box in bl.  The list is retruned in the
               argument clst.  For each intersection that includes tagged
               points, construct a new cluster by removing those points from
               this object.  Empty intersections or those that contain no
               tagged points will not generate a new cluster.
               Note that this cluster will be modified and possibly become
               invalid in the process.
    */
    void distribute (ClusterList&     clst,
                     const BoxDomain& bd);
    //
    //@ManDoc: Compute ratio of tagged to total number of points in cluster.
    //
    Real eff () const
    {
        BL_ASSERT(ok());
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
//@Man:
//@Memo: A list of Cluster objects
/*@Doc:

  A container class for Cluster.
*/

class ClusterList 
{
public:
    //
    //@ManDoc: The default constructor.
    //
    ClusterList ()
        :
        lst() {}
    //
    //@ManDoc: Construct a list containing Cluster(pts,len).
    //
    ClusterList (IntVect* pts,
                 long     len)
    {
        lst.append(new Cluster(pts,len));
    }
    //
    //@ManDoc: The destructor.
    //
    ~ClusterList ();
    //
    //@ManDoc: Return number of clusters in list.
    //
    int length () const { return lst.length(); }
    //
    //@ManDoc: Add cluster to end of list.
    //
    void append (Cluster* c) { lst.append(c); }
    //
    //@ManDoc: Return array of boxes corresponding to clusters.
    //
    BoxArray boxArray () const;
    //
    //@ManDoc: Return array of boxes corresponding to clusters in argument.
    //
    void boxArray (BoxArray& ba) const;
    //
    //@ManDoc: Return list of boxes corresponding to clusters.
    //
    BoxList boxList() const;
    //
    //@ManDoc: Return list of boxes corresponding to clusters in argument.
    //
    void boxList (BoxList& blst) const;
    //
    //@ManDoc: Chop all clusters in list that have poor efficiency.
    //
    void chop (Real eff);

    /*@ManDoc: Intersect clusters with BoxDomain to insure cluster
               boxes are interior to domain.
    */
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

#endif /*_Cluster_H_*/
