#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

#include <map>
#include <limits>
#include <memory>
#include <cstddef>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {

class BoxArray;
class MultiFab;
template <typename T> class FabArray;
class FabArrayBase;

/**
* \brief Calculates the distribution of FABs to MPI processes.
*
*  This class calculates the distribution of FABs to MPI processes in a
*  FabArray in a multi-processor environment.  By distribution is meant what
*  MPI process in the multi-processor environment owns what FAB.  Only the BoxArray
*  on which the FabArray is built is used in determining the distribution.
*  The three types of distributions supported are round-robin, knapsack, and SFC.
*  In the round-robin distribution FAB i is owned by CPU i%N where N is total
*  number of CPUs.  In the knapsack distribution the FABs are partitioned
*  across CPUs such that the total volume of the Boxes in the underlying
*  BoxArray are as equal across CPUs as is possible.  The SFC distribution is
*  based on a space filling curve.
*/

class DistributionMapping
{
  public:

    template <typename T> friend class FabArray;
    friend class FabArrayBase;

    //! The distribution strategies
    enum Strategy { UNDEFINED = -1, ROUNDROBIN, KNAPSACK, SFC, RRSFC };

    //! The default constructor.
    DistributionMapping ();

    //! The copy constructor.
    DistributionMapping (const DistributionMapping& rhs) = default;

    //! The copy assignment operator.
    DistributionMapping& operator= (const DistributionMapping& rhs) = default;

    //! The move constructor.
    DistributionMapping (DistributionMapping&& rhs) noexcept = default;

    //! The move assignment operator.
    DistributionMapping& operator= (DistributionMapping&& rhs) noexcept = default;

   //! The destructor.
    ~DistributionMapping() noexcept = default;

    /**
    * \brief Create an object with the specified mapping.
    */
    explicit DistributionMapping (const Vector<int>& pmap);
    explicit DistributionMapping (Vector<int>&& pmap) noexcept;
    //! Build mapping out of BoxArray over nprocs processors.
    explicit DistributionMapping (const BoxArray& boxes,
				  int nprocs = ParallelDescriptor::NProcs());
    /**
    * \brief This is a very specialized distribution map.
    * Do NOT use it unless you really understand what it does.
    */
    DistributionMapping (const DistributionMapping& d1,
                         const DistributionMapping& d2);

    /**
    * \brief Build mapping out of BoxArray over nprocs processors.
    * You need to call this if you built your DistributionMapping
    * with the default constructor.
    */
    void define (const BoxArray& boxes, int nprocs = ParallelDescriptor::NProcs());
    /**
    * \brief Build mapping out of an Array of ints. You need to call this if you
    * built your DistributionMapping with the default constructor.
    */
    void define (const Vector<int>& pmap);
    void define (Vector<int>&& pmap) noexcept;
    /**
    * \brief Returns a constant reference to the mapping of boxes in the
    * underlying BoxArray to the CPU that holds the FAB on that Box.
    * ProcessorMap()[i] is an integer in the interval [0, NCPU) where
    * NCPU is the number of CPUs being used.
    */
    const Vector<int>& ProcessorMap () const noexcept;

    //! Length of the underlying processor map.
    long size () const noexcept { return m_ref->m_pmap.size(); }
    long capacity () const noexcept { return m_ref->m_pmap.capacity(); }
    bool empty () const noexcept { return m_ref->m_pmap.empty(); }

    //! Number of references to this DistributionMapping
    long linkCount () const noexcept { return m_ref.use_count(); }

    //! Equivalent to ProcessorMap()[index].
    int operator[] (int index) const noexcept { return m_ref->m_pmap[index]; }

    //! Set/get the distribution strategy.
    static void strategy (Strategy how);

    static Strategy strategy ();

    //! Set/get the space filling curve threshold.
    static void SFC_Threshold (int n);

    static int SFC_Threshold ();

    //! Are the distributions equal?
    bool operator== (const DistributionMapping& rhs) const noexcept;

    //! Are the distributions different?
    bool operator!= (const DistributionMapping& rhs) const noexcept;

    void SFCProcessorMap(const BoxArray& boxes, const std::vector<long>& wgts, int nprocs,
                         bool sort=true);
    void KnapSackProcessorMap(const std::vector<long>& wgts, int nprocs,
                              Real* efficiency = 0,
			      bool do_full_knapsack = true,
			      int nmax = std::numeric_limits<int>::max());
    void RoundRobinProcessorMap(int nboxes, int nprocs);
    void RoundRobinProcessorMap(const std::vector<long>& wgts, int nprocs);

    /**
    * \brief Initializes distribution strategy from ParmParse.
    *
    * ParmParse options are:
    *
    *   DistributionMapping.strategy = ROUNDROBIN
    *   DistributionMapping.strategy = KNAPSACK
    *   DistributionMapping.strategy = SFC
    *   DistributionMapping.strategy = RRFC
    */
    static void Initialize ();

    static void Finalize ();

    static bool SameRefs (const DistributionMapping& lhs,
                          const DistributionMapping& rhs)
		  { return lhs.m_ref == rhs.m_ref; }

    static DistributionMapping makeKnapSack   (const MultiFab& weight,
                                               int nmax=std::numeric_limits<int>::max());
    static DistributionMapping makeKnapSack   (const Vector<Real>& rcost);

    static DistributionMapping makeRoundRobin (const MultiFab& weight);
    static DistributionMapping makeSFC        (const MultiFab& weight, bool sort=true);

    /**
    * if use_box_vol is true, weight boxes by their volume in Distribute
    * otherwise, all boxes will be treated with equal weight
    */
    static std::vector<std::vector<int> > makeSFC (const BoxArray& ba, bool use_box_vol=true);

private:

    const Vector<int>& getIndexArray ();
    const std::vector<bool>& getOwnerShip ();

    //! Ways to create the processor map.
    void RoundRobinProcessorMap (const BoxArray& boxes, int nprocs);
    void KnapSackProcessorMap   (const BoxArray& boxes, int nprocs);
    void SFCProcessorMap        (const BoxArray& boxes, int nprocs);
    void RRSFCProcessorMap      (const BoxArray& boxes, int nprocs);

    using LIpair = std::pair<long,int>;

    struct LIpairLT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const noexcept
            {
                return lhs.first < rhs.first;
            }
    };

    struct LIpairGT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const noexcept
            {
                return lhs.first > rhs.first;
            }
    };

    static void Sort (std::vector<LIpair>& vec, bool reverse);

    void RoundRobinDoIt (int                  nboxes,
                         int                  nprocs,
                         std::vector<LIpair>* LIpairV = 0);

    void KnapSackDoIt (const std::vector<long>& wgts,
                       int                      nprocs,
                       Real&                    efficiency,
                       bool                     do_full_knapsack,
		       int                      nmax=std::numeric_limits<int>::max());

    void SFCProcessorMapDoIt (const BoxArray&          boxes,
                              const std::vector<long>& wgts,
                              int                      nprocs,
                              bool                     sort=true);

    void RRSFCDoIt           (const BoxArray&          boxes,
                              int                      nprocs);

    //! Least used ordering of CPUs (by # of bytes of FAB data).
    void LeastUsedCPUs (int nprocs, Vector<int>& result);
    /**
    * \brief rteam: Least used ordering of Teams
    * rworker[i]: Least used ordering of team workers for Team i
    */
    void LeastUsedTeams (Vector<int>& rteam, Vector<Vector<int> >& rworker, int nteams, int nworkers);

    //! A useful typedef.
    typedef void (DistributionMapping::*PVMF)(const BoxArray &, int);

    //! Everyone uses the same Strategy -- defaults to SFC.
    static Strategy m_Strategy;
    /**
    * \brief Pointer to one of the CreateProcessorMap() functions.
    * Corresponds to the one specified by m_Strategy.
    */
    static PVMF m_BuildMap;

    struct Ref
    {
        friend class DistributionMapping;

        //! Constructors to match those in DistributionMapping ....
        Ref () noexcept {}

        explicit Ref (int len) : m_pmap(len) {}

        explicit Ref (const Vector<int>& pmap) : m_pmap(pmap) {}

        explicit Ref (Vector<int>&& pmap) noexcept : m_pmap(std::move(pmap)) {}

	//! dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

        void clear () { m_pmap.clear();  m_index_array.clear();   m_ownership.clear(); }

        Vector<int> m_pmap; //!< index array for all boxes
        Vector<int> m_index_array;  //!< index array for local boxes owned by the team
        std::vector<bool> m_ownership; //!< true ownership
    };
    //
    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<Ref> m_ref;

public:
    struct RefID {
        RefID () noexcept : data(nullptr) {}
        explicit RefID (Ref* data_) noexcept : data(data_) {}
        bool operator<  (const RefID& rhs) const noexcept { return std::less<Ref*>()(data,rhs.data); }
        bool operator== (const RefID& rhs) const noexcept { return data == rhs.data; }
        bool operator!= (const RefID& rhs) const noexcept { return data != rhs.data; }
	const Ref *dataPtr() const noexcept { return data; }
	void PrintPtr(std::ostream &os) const { os << data << '\n'; }
        friend std::ostream& operator<< (std::ostream& os, const RefID& id);
   private:
        Ref* data;
    };

    //! This gives a unique ID of the reference, which is different from dmID above.
    RefID getRefID () const noexcept { return RefID { m_ref.get() }; }
};

//! Our output operator.
std::ostream& operator<< (std::ostream& os, const DistributionMapping& pmap);

std::ostream& operator<< (std::ostream& os, const DistributionMapping::RefID& id);

}

#endif /*BL_DISTRIBUTIONMAPPING_H*/
