
#ifndef AMREX_Cluster_H_
#define AMREX_Cluster_H_
#include <AMReX_Config.H>

#include <AMReX_BoxList.H>
#include <AMReX_REAL.H>

#include <list>

namespace amrex {

class BoxDomain;
class BoxArray;
class ClusterList;


/**
* \brief A cluster of tagged cells.
*
* Utility class for tagging error cells.
*/

class Cluster
{
public:

    /**
    * \brief The default constructor -- builds invalid Cluster.
    */
    Cluster () noexcept = default;

    /**
    * \brief Construct a cluster from an array of IntVects.
    * The Cluster object does NOT take over memory management
    * of array; i.e. that's the user's ultimate responsibility.
    *
    * \param a
    * \param len
    */
    Cluster (IntVect* a, Long len) noexcept;

    /**
    * \brief Construct new cluster by removing all points from c that lie
    * in box b.  Cluster c is modified and may become invalid.
    *
    * \param c
    * \param b
    */
    Cluster (Cluster& c, const Box& b);

    ~Cluster () = default;
    Cluster (const Cluster&) = delete;
    Cluster (Cluster&&) = delete;
    Cluster& operator= (const Cluster&) = delete;
    Cluster& operator= (Cluster&&) = delete;

    /**
    * \brief Return minimal box containing all tagged points.
    */
    [[nodiscard]] const Box& box () const noexcept { return m_bx; }

    /**
    * \brief Does cluster contain any points?
    */
    [[nodiscard]] bool ok () const noexcept { return m_ar != nullptr && m_len > 0; }

    /**
    * \brief Returns number of tagged points in cluster.
    */
    [[nodiscard]] Long numTag () const noexcept { return m_len; }

    /**
    * \brief Return number of tagged points in intersection of cluster and Box b.
    *
    * \param b
    */
    [[nodiscard]] Long numTag (const Box& b) const noexcept;

    /**
    * \brief This operation splits a cluster into two pieces by selecting
    * a cutting plane with certain optimal characteristics then
    * dividing the tagged points into clusters on either side of the
    * plane.  One cluster is returned as a new object the other
    * is the modified calling object.  This is called by chop(eff)
    */
    Cluster* chop ();


    /**
    * \brief This version of chop has slightly different logic - in this case if a
    * cut results in two boxes with the same grid efficiency as the original
    * box then the cut is reverted and a cut in a different direction is chosen
    * This is called by new_chop(eff)
    */
    Cluster* new_chop ();

    /**
    * \brief Constructs a list of cluster objects obtained by intersecting
    * this cluster with each box in bl.  The list is returned in the
    * argument clst.  For each intersection that includes tagged
    * points, construct a new cluster by removing those points from
    * this object.  Empty intersections or those that contain no
    * tagged points will not generate a new cluster.
    * Note that this cluster will be modified and possibly become
    * invalid in the process.
    *
    * \param clst
    * \param bd
    */
    void distribute (ClusterList&     clst,
                     const BoxDomain& bd);

    /**
    * \brief Compute ratio of tagged to total number of points in cluster.
    */
    [[nodiscard]] Real eff () const noexcept {
        BL_ASSERT(ok());
        return static_cast<Real>(double(numTag()) / m_bx.d_numPts());
    }

private:

    /**
    * \brief Compute and store minimal box containing tagged points.
    */
    void minBox () noexcept;

    //! The data.
    Box      m_bx;
    IntVect* m_ar = nullptr;
    Long     m_len = 0;
};


/**
* \brief A list of Cluster objects
*
* A container class for Cluster.
*/

class ClusterList
{
public:

    /**
    * \brief The default constructor.
    */
    ClusterList () = default;

    /**
    * \brief Construct a list containing Cluster(pts,len).
    *
    * \param pts
    * \param len
    */
    ClusterList (IntVect* pts, Long len);

    /**
    * \brief The destructor.
    */
    ~ClusterList ();

    ClusterList (const ClusterList&) = delete;
    ClusterList (ClusterList&&) = delete;
    ClusterList& operator= (const ClusterList&) = delete;
    ClusterList& operator= (ClusterList&&) = delete;

    /**
    * \brief Return number of clusters in list.
    */
    [[nodiscard]] int length () const { return static_cast<int>(lst.size()); }

    /**
    * \brief Add cluster to end of list.
    *
    * \param c
    */
    void append (Cluster* c) { lst.push_back(c); }

    /**
    * \brief Return array of boxes corresponding to clusters.
    */
    [[nodiscard]] BoxArray boxArray () const;

    /**
    * \brief Return array of boxes corresponding to clusters in argument.
    *
    * \param ba
    */
    void boxArray (BoxArray& ba) const;

    /**
    * \brief Return list of boxes corresponding to clusters.
    */
    [[nodiscard]] BoxList boxList() const;

    /**
    * \brief Return list of boxes corresponding to clusters in argument.
    *
    * \param blst
    */
    void boxList (BoxList& blst) const;

    /**
    * \brief Chop all clusters in list that have poor efficiency.
    *
    * \param eff
    */
    void chop (Real eff);

    /**
    * \brief Chop all clusters in list that have poor efficiency.
    * This version calls new_chop() instead of chop()
    *
    * \param eff
    */
    void new_chop (Real eff);

    /**
    * \brief Intersect clusters with BoxArray to insure cluster
    * boxes are interior to the domain of BoxArray.  Note that
    * ba is modified during the process.
    *
    * \param ba
    */
    void intersect (BoxArray& ba);

private:

    //! The data.
    std::list<Cluster*> lst;
};

}

#endif /*_Cluster_H_*/
