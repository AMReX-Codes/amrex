#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

#include <map>
#include <limits>
#include <memory>
#include <cstddef>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {

class BoxArray;
class MultiFab;
template <typename T> class FabArray;

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

    //! The distribution strategies
    enum Strategy { UNDEFINED = -1, ROUNDROBIN, KNAPSACK, SFC, PFC, RRSFC };

    //! The default constructor.
    DistributionMapping ();

    /**
    * \brief Create an object with the specified mapping.
    */
    explicit DistributionMapping (const Array<int>& pmap,
				  ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //! Build mapping out of BoxArray over nprocs processors.
    explicit DistributionMapping (const BoxArray& boxes,
				  int nprocs = ParallelDescriptor::NProcs(), 
				  ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    /**
    * \brief This is a very specialized distribution map.
    * Do NOT use it unless you really understand what it does.
    */
    DistributionMapping (const DistributionMapping& d1,
                         const DistributionMapping& d2);

    //! The destructor.
    ~DistributionMapping ();

    //! Copy constructor.
    DistributionMapping (const DistributionMapping& rhs);

    //! Assignment operator.
    DistributionMapping& operator= (const DistributionMapping& rhs);

    //! Move constructor.
    DistributionMapping (DistributionMapping&& rhs) noexcept;

    /**
    * \brief Build mapping out of BoxArray over nprocs processors.
    * You need to call this if you built your DistributionMapping
    * with the default constructor.
    */
    void define (const BoxArray& boxes, int nprocs = ParallelDescriptor::NProcs(),
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    /**
    * \brief Build mapping out of an Array of ints. You need to call this if you
    * built your DistributionMapping with the default constructor.
    */
    void define (const Array<int>& pmap);
    /**
    * \brief Returns a constant reference to the mapping of boxes in the
    * underlying BoxArray to the CPU that holds the FAB on that Box.
    * ProcessorMap()[i] is an integer in the interval [0, NCPU) where
    * NCPU is the number of CPUs being used.
    */
    const Array<int>& ProcessorMap () const;
    //
    ParallelDescriptor::Color color() const { return m_color; }

    //! Length of the underlying processor map.
    long size () const { return m_ref->m_pmap.size(); }
    long capacity () const { return m_ref->m_pmap.capacity(); }
    bool empty () const { return m_ref->m_pmap.empty(); }

    //! Number of references to this DistributionMapping
    long linkCount () const { return m_ref.use_count(); }

    //! Equivalent to ProcessorMap()[index].
    int operator[] (int index) const { return m_ref->m_pmap[index]; }

    // Replace the cached processor map.  This is to support FabArray::MoveAllFabs
    //   All FabArrays using the cached map must move their fabs before
    //   calling this function.
    void ReplaceCachedProcessorMap(const Array<int> &newProcMap);

    //! Set/get the distribution strategy.
    static void strategy (Strategy how);

    static Strategy strategy ();

    //! Set/get the space filling curve threshold.
    static void SFC_Threshold (int n);

    static int SFC_Threshold ();
 
    //! Are the distributions equal?
    bool operator== (const DistributionMapping& rhs) const;

    //! Are the distributions different?
    bool operator!= (const DistributionMapping& rhs) const;

    void SFCProcessorMap(const BoxArray& boxes, const std::vector<long>& wgts,
                         int nprocs);
    void PFCProcessorMap(const BoxArray& boxes, const std::vector<long>& wgts,
                         int nprocs);
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
    *   DistributionMapping.strategy = PFC
    *   DistributionMapping.strategy = RRFC
    */
    static void Initialize ();

    static void Finalize ();

    static bool SameRefs (const DistributionMapping& lhs,
                          const DistributionMapping& rhs)
		  { return lhs.m_ref == rhs.m_ref; }

    static void PrintDiagnostics(const std::string &filename);

    //! Initialize the topological proximity map
    static void InitProximityMap(bool makeMap = false, bool reinit = false);
    static int NHops(const Box &tbox, const IntVect &ivfrom, const IntVect &ivto);
    static int ProcNumberFromRank(const int rank);
    static std::vector<int> RanksFromProcNumber(const int procnum);
    static IntVect TopIVFromProcNumber(const int procnum);
    static std::vector<int> ProcNumbersFromTopIV(const IntVect &iv);
    static IntVect TopIVFromRank(const int rank);
    static std::vector<int> RanksFromTopIV(const IntVect &iv);
    static std::string GetProcName();
    static int GetProcNumber();
    static int ProximityMap(const int rank)   { return proximityMap[rank];   }
    static int ProximityOrder(const int rank) { return proximityOrder[rank]; }

#if !defined(BL_NO_FORT)
    static void ReadCheckPointHeader(const std::string &filename,
				     Array<IntVect>  &refRatio,
                                     Array<BoxArray> &allBoxes);
#endif
    static void PFCMultiLevelMap(const Array<IntVect>  &refRatio,
                                 const Array<BoxArray> &allBoxes);

    static Array<Array<int> > MultiLevelMapPFC(const Array<IntVect>  &refRatio,
                                               const Array<BoxArray> &allBoxes,
					       int maxgrid);
    static Array<Array<int> > MultiLevelMapRandom(const Array<IntVect>  &refRatio,
                                                  const Array<BoxArray> &allBoxes,
						  int maxgrid,
						  int maxRank = -1, int minRank = 0);
    static Array<Array<int> > MultiLevelMapKnapSack(const Array<IntVect>  &refRatio,
                                                    const Array<BoxArray> &allBoxes,
						    int maxgrid);

    static int NDistMaps() { return nDistMaps; }
    static void SetNDistMaps(int ndm) { nDistMaps = ndm; }
    int DistMapID() const { return dmID; }
    void SetDistMapID(int dmid) { dmID = dmid; }

    bool Check() const;

#ifdef BL_USE_MPI
    static Array<int> TranslateProcMap(const Array<int> &pm_old, const MPI_Group group_new, const MPI_Group group_old);
#endif

    static DistributionMapping makeKnapSack   (const MultiFab& weight,
                                               int nmax=std::numeric_limits<int>::max());
    static DistributionMapping makeRoundRobin (const MultiFab& weight);
    static DistributionMapping makeSFC        (const MultiFab& weight, const BoxArray& boxes);

private:

    //! Ways to create the processor map.
    void RoundRobinProcessorMap (const BoxArray& boxes, int nprocs);
    void KnapSackProcessorMap   (const BoxArray& boxes, int nprocs);
    void SFCProcessorMap        (const BoxArray& boxes, int nprocs);
    void PFCProcessorMap        (const BoxArray& boxes, int nprocs);
    void RRSFCProcessorMap      (const BoxArray& boxes, int nprocs);

    using LIpair = std::pair<long,int>;

    struct LIpairLT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const
            {
                return lhs.first < rhs.first;
            }
    };

    struct LIpairGT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const
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
                              int                      nprocs);

    void PFCProcessorMapDoIt (const BoxArray&          boxes,
                              const std::vector<long>& wgts,
                              int                      nprocs);

    void RRSFCDoIt           (const BoxArray&          boxes,
                              int                      nprocs);

    //! Current # of bytes of FAB data.
    static void CurrentBytesUsed (int nprocs, Array<long>& result);
    static void CurrentCellsUsed (int nprocs, Array<long>& result);

    //! Least used ordering of CPUs (by # of bytes of FAB data).
    void LeastUsedCPUs (int nprocs, Array<int>& result);
    /**
    * \brief rteam: Least used ordering of Teams
    * rworker[i]: Least used ordering of team workers for Team i
    */
    void LeastUsedTeams (Array<int>& rteam, Array<Array<int> >& rworker, int nteams, int nworkers);

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
        Ref () {}

        explicit Ref (int len) : m_pmap(len) {}

        explicit Ref (const Array<int>& pmap) : m_pmap(pmap) {}

	//! dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

	//! Local data -- our processor map.
        Array<int> m_pmap;
    };
    //
    // The data -- a reference-counted pointer to a Ref.
    //
    std::shared_ptr<Ref> m_ref;
    ParallelDescriptor::Color m_color;
    int dmID;

    //! Topological proximity map
    static long totalCells;
    static Real bytesPerCell;
    static Array<int> proximityMap;    // i == rank, pMap[i]   == proximity mapped rank
    static Array<int> proximityOrder;  // i == rank, pOrder[i] == proximity mapped order
    static Array<long> totalBoxPoints;  // i == rank

    static int nDistMaps;

public:
    struct RefID {
        RefID () : data(nullptr) {}
        explicit RefID (Ref* data_) : data(data_) {}
        bool operator<  (const RefID& rhs) const { return std::less<Ref*>()(data,rhs.data); }
        bool operator== (const RefID& rhs) const { return data == rhs.data; }
        bool operator!= (const RefID& rhs) const { return data != rhs.data; }
	const Ref *dataPtr() const { return data; }
	void PrintPtr(std::ostream &os) const { os << data << '\n'; }
   private:
        Ref* data;
    };

    //! This gives a unique ID of the reference, which is different from dmID above.
    RefID getRefID () const { return RefID { m_ref.get() }; }
};

//! Our output operator.
std::ostream& operator<< (std::ostream& os, const DistributionMapping& pmap);

}

#endif /*BL_DISTRIBUTIONMAPPING_H*/
