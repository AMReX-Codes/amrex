
#ifndef AMREX_StateData_H_
#define AMREX_StateData_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFCopyDescriptor.H>
#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_StateDescriptor.H>

/**
 * \file AMReX_StateData.H
 * \brief Manages the old/new time slices of a state descriptor on a single AMR level.
 */

#include <memory>

namespace amrex {

class StateDataPhysBCFunct;


/**
* \brief Current and previous level-time data.
*
* StateData holds state data on a level for the current and previous time step.
*/
class StateData
{
    friend class StateDataPhysBCFunct;

public:

    /**
    * \brief The default constructor.
    */
    StateData ();

    /**
     * \brief Construct and allocate state storage for a level.
     *
     * \param p_domain  Domain box of the AMR problem (used for BCs).
     * \param grds      Cell-centered BoxArray for this level.
     * \param dm        Distribution mapping describing ownership of grids.
     * \param d         State descriptor providing metadata (components, BCs, etc.).
     * \param cur_time  Current simulation time.
     * \param dt        Initial time step (dt_old = dt_new = dt).
     * \param factory   Fab factory controlling allocation (EB-aware, etc.).
     */
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const DistributionMapping& dm,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt,
               const FabFactory<FArrayBox>& factory);

    /**
    * \brief The destructor.
    */
    ~StateData ();

    //! Move constructor; transfers ownership of stored MultiFabs and metadata.
    StateData (StateData&& rhs) noexcept;
    //! Copy-assignment; deep-copies metadata and allocates new FAB storage.
    StateData& operator= (StateData const& rhs);

    //! Copy construction is disabled to avoid ambiguous ownership of state buffers.
    StateData (StateData const& rhs) = delete;
    //! Move assignment is disabled (use define()/swap helpers instead).
    StateData& operator= (StateData && rhs) = delete;

    /**
     * \brief Initialize the object after default construction.
     *
     * \param p_domain  Domain box of the AMR problem (used for BCs).
     * \param grds      Cell-centered BoxArray for this level.
     * \param dm        Distribution mapping for that grid layout.
     * \param d         State descriptor describing component layout.
     * \param cur_time  Current simulation time.
     * \param dt        Initial dt (applied to both old/new slots).
     * \param factory   Fab factory used for allocation.
     */
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const DistributionMapping& dm,
                 const StateDescriptor& d,
                 Real                   cur_time,
                 Real                   dt,
                 const FabFactory<FArrayBox>& factory);

    /**
     * \brief Copy the "old" time slice from \p state, allocating storage if needed.
     *
     * \param state Source StateData whose old data/time level are cloned.
     */
    void copyOld (const StateData& state);

    /**
     * \brief Copy the "new" time slice from \p state, allocating storage if needed.
     *
     * \param state Source StateData whose new data/time level are cloned.
     */
    void copyNew (const StateData& state);

    /**
    * \brief Allocates space for old timestep data.
    */
    void allocOldData ();

    /**
    * \brief Deletes the space used by the old timestep data.
    */
    void removeOldData () { old_data.reset(); }

    /**
    * \brief Reverts back to initial state.
    */
    void reset ();

    /**
    * \brief Old data becomes new data and new time is incremented by \p dt.
    *
    * \param dt Time step separating the old and new slots.
    */
    void swapTimeLevels (Real dt);

    /**
    * \brief Swaps old data with a new MultiFab.
    * Deletes the previous old data.
    *
    * \param mf MultiFab that will replace the stored old data (ownership transferred).
    */
    void replaceOldData ( MultiFab&& mf );

    /**
    * \brief Swaps old data with another StateData.
    * Does not delete the previous old data.
    *
    * \param s Peer whose old data will be exchanged with this instance.
    */
    void replaceOldData ( StateData& s );

    /**
    * \brief Swaps new data with a new MultiFab.
    * Deletes the previous new data.
    *
    * \param mf MultiFab that will replace the stored new data (ownership transferred).
    */
    void replaceNewData ( MultiFab&& mf );

    /**
    * \brief Swaps new data with another StateData.
    * Does not delete the previous new data.
    *
    * \param s Peer whose new data will be exchanged with this instance.
    */
    void replaceNewData ( StateData& s );


    /**
    * \brief Sets time of old and new data.
    *
    * \param time    New-time value.
    * \param dt_old  Dt associated with the old time interval.
    * \param dt_new  Dt associated with the new time interval.
    */
    void setTimeLevel (Real time,
                       Real dt_old,
                       Real dt_new);

    /**
    * \brief Sets time of old data.
    *
    * \param time Old-time value to record.
    */
    void setOldTimeLevel (Real time);

    /**
    * \brief Sets time of new data.
    *
    * \param time New-time value to record.
    */
    void setNewTimeLevel (Real time);

    //! Set both start/stop endpoints of the new-time interval to \p time.
    void syncNewTimeLevel (Real time);

    /**
     * \brief Register the state data with a MultiFabCopyDescriptor for later interpolation/fills.
     *
     * \param multiFabCopyDesc Descriptor that tracks MultiFabs participating in copy/interp operations.
     * \param mfid             Output vector receiving identifiers for the registered MultiFabs.
     */
    void RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                       Vector<MultiFabId>&      mfid);

    /**
     * \brief Request interpolation/copied data covering \p subbox.
     *
     * \param multiFabCopyDesc Descriptor returned by RegisterData().
     * \param mfid             Identifiers previously produced by RegisterData().
     * \param returnedUnfillableBoxes Optional output list of boxes that could not be filled.
     * \param returnedFillBoxIds      Output fill-box identifiers referencing queued requests.
     * \param subbox           Region that needs to be populated.
     * \param time             Physical time at which data should be sampled.
     * \param src_comp         Source component index in the state data.
     * \param dest_comp        Destination component index inside the caller's buffers.
     * \param num_comp         Number of components requested.
     * \param extrap           If true, allow extrapolation beyond old/new time intervals.
     */
    void InterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
                       Vector<MultiFabId>&      mfid,
                       BoxList*                returnedUnfillableBoxes,
                       Vector<FillBoxId>&       returnedFillBoxIds,
                       const Box&              subbox,
                       Real                    time,
                       int                     src_comp,
                       int                     dest_comp,
                       int                     num_comp,
                       bool                    extrap = false);

    /**
     * \brief Complete outstanding interpolation requests and fill \p dest.
     *
     * \param fabCopyDesc Descriptor produced by RegisterData().
     * \param mfid        Identifiers previously produced by RegisterData().
     * \param fillBoxIds  Identifiers returned by InterpAddBox().
     * \param dest        Destination FAB.
     * \param time        Physical time at which data should be sampled.
     * \param src_comp    Source component offset in the state data.
     * \param dest_comp   Destination component offset inside \p dest.
     * \param num_comp    Number of components to fill.
     * \param extrap      If true, allow extrapolation beyond old/new time intervals.
     */
    void InterpFillFab (MultiFabCopyDescriptor&  fabCopyDesc,
                        const Vector<MultiFabId>& mfid,
                        const Vector<FillBoxId>&  fillBoxIds,
                        FArrayBox&               dest,
                        Real                     time,
                        int                      src_comp,
                        int                      dest_comp,
                        int                      num_comp,
                        bool                     extrap = false);


    /**
    * \brief Set physical boundary values for a FAB extracted from the MultiFab.
    *
    * \param dest         Destination FAB to modify.
    * \param time         Physical time at which BCs are evaluated.
    * \param dx           Cell spacing array.
    * \param prob_domain  Real-space domain extents.
    * \param dest_comp    Destination component offset in \p dest.
    * \param src_comp     Source component offset in the state data.
    * \param num_comp     Number of components to update (default 1).
    */
    void FillBoundary (FArrayBox&     dest,
                       Real           time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dest_comp,
                       int            src_comp,
                       int            num_comp = 1);

    //! Variant that fills data on region \p bx into \p dest at components starting \p dest_comp,
    //! sampling state components beginning at \p src_comp for \p num_comp values using geometry \p geom at \p time.
    void FillBoundary (Box const&      bx,
                       FArrayBox&      dest,
                       Real            time,
                       Geometry const& geom,
                       int             dest_comp,
                       int             src_comp,
                       int             num_comp);

    /**
    * \brief Write the state data to a checkpoint file.
    *
    * \param name         State name (prefix for component folders).
    * \param fullpathname Directory that will hold the data.
    * \param os           Output stream for metadata.
    * \param how          VisMF write mode.
    * \param dump_old     Whether to include the old-time slot.
    */
    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os,
                     VisMF::How         how,
                     bool               dump_old = true);

    /**
    * \brief Restart with domain box, grids, and dmap provided.
    *
    * \param is        Input stream tied to the checkpoint data.
    * \param p_domain  Domain box for the level.
    * \param grds      Cell-centered grid layout.
    * \param dm        Distribution mapping describing ownership of grids.
    * \param factory   Fab factory controlling allocation strategy.
    * \param d         State descriptor describing components/BCs.
    * \param chkfile   Name of the checkpoint file (for logging).
    */
    void restart (std::istream&          is,
                  const Box&             p_domain,
                  const BoxArray&        grds,
                  const DistributionMapping& dm,
                  const FabFactory<FArrayBox>& factory,
                  const StateDescriptor& d,
                  const std::string&     chkfile);

    /**
    * \brief Initialize from another StateData sharing descriptor \p d.
    *
    * \param d   State descriptor describing the components.
    * \param rhs Source data that will be cloned.
    */
    void restart (const StateDescriptor& d,
                  const StateData& rhs);

    /**
    * \brief Returns the StateDescriptor.
    */
    const StateDescriptor* descriptor () const noexcept { return desc; }

    /**
    * \brief Returns the valid domain.
    */
    const Box& getDomain () const noexcept { return domain; }

    /**
    * \brief Returns the BoxArray.
    */
    const BoxArray& boxArray () const noexcept { return grids; }

    //! Replace the stored BoxArray layout with \p a_grids (no redistribution occurs).
    void setBoxArray (BoxArray const& a_grids) noexcept { grids = a_grids; }

    const DistributionMapping& DistributionMap () const noexcept { return dmap; }

    //! Replace the stored DistributionMapping with \p new_dmap (no redistribution occurs).
    void setDistributionMap ( DistributionMapping& new_dmap ) noexcept { dmap = new_dmap; }

    const FabFactory<FArrayBox>& Factory () const noexcept { return *m_factory; }

    //! Replace the internal factory with a clone of \p a_factory.
    void setFactory (FabFactory<FArrayBox> const& a_factory) {
        m_factory.reset(a_factory.clone());
    }

    /**
    * \brief Returns the current time.
    */
    Real curTime () const noexcept {
        return (desc->timeType() == StateDescriptor::Point) ?
            new_time.stop : 0.5_rt*(new_time.start + new_time.stop);
    }

    /**
    * \brief Returns the previous time.
    */
    Real prevTime () const noexcept {
        return (desc->timeType() == StateDescriptor::Point) ?
            old_time.stop : 0.5_rt*(old_time.start + old_time.stop);
    }

    /**
    * \brief Returns the new data.
    */
    MultiFab& newData () noexcept { BL_ASSERT(new_data != nullptr); return *new_data; }

    /**
    * \brief Returns the new data.
    */
    const MultiFab& newData () const noexcept { BL_ASSERT(new_data != nullptr); return *new_data; }

    /**
    * \brief Returns the old data.
    */
    MultiFab& oldData () noexcept { BL_ASSERT(old_data != nullptr); return *old_data; }

    /**
    * \brief Returns the old data.
    */
    const MultiFab& oldData () const noexcept { BL_ASSERT(old_data != nullptr); return *old_data; }

    //! Returns the FAB of new data at grid index \p i.
    FArrayBox& newGrid (int i) noexcept { BL_ASSERT(new_data != nullptr); return (*new_data)[i]; }

    //! Returns the FAB of old data at grid index \p i.
    FArrayBox& oldGrid (int i) noexcept { BL_ASSERT(old_data != nullptr); return (*old_data)[i]; }

    /**
    * \brief Returns boundary conditions of specified component on the specified grid.
    *
    * \param comp Component index.
    * \param i    Grid index.
    */
    BCRec getBC (int comp, int i) const noexcept;

    /**
    * \brief Prints out the time interval.
    *
    * \param os Output stream receiving a textual summary of old/new times.
    */
    void printTimeInterval (std::ostream& os) const;

    /**
    * \brief True if there is any old data available.
    */
    bool hasOldData () const noexcept { return old_data != nullptr; }

    /**
    * \brief True if there is any new data available.
    */
    bool hasNewData () const noexcept { return new_data != nullptr; }

    /**
     * \brief Gather data pointers and their time stamps corresponding to \p time.
     *
     * \param data      Output list receiving pointers to MultiFabs.
     * \param datatime  Output list receiving the associated times.
     * \param time      Requested interpolation time.
     */
    void getData (Vector<MultiFab*>& data,
                  Vector<Real>& datatime,
                  Real time) const;

    /**
    * \brief Get the Arena used.
    */
    Arena* getArena () const noexcept { return arena; }

    //! Set the Arena used for subsequent allocations to \p ar.
    void setArena (Arena* ar) noexcept { arena = ar; }

    /**
    * \brief These facilitate prereading FabArray headers to avoid
    * synchronization when reading multiple FabArrays
    */
    //! View of FabArray header filenames recorded during the last checkpoint write.
    static const Vector<std::string> &FabArrayHeaderNames() { return fabArrayHeaderNames; }
    //! Clear the cached list of FabArray header names.
    static void ClearFabArrayHeaderNames() { fabArrayHeaderNames.clear(); }

    //! Install a pointer \p fahmp to the global map storing pre-read FabArray headers.
    static void SetFAHeaderMapPtr(std::map<std::string, Vector<char> > *fahmp) { faHeaderMap = fahmp; }


private:

    std::unique_ptr<FabFactory<FArrayBox> > m_factory;

    struct TimeInterval
    {
        Real start;
        Real stop;
    };

    //! Pointer to data descriptor.
    const StateDescriptor* desc{nullptr};

    //! Problem domain.
    Box domain;

    //! Grids defined at this level.
    BoxArray grids;
    //
    DistributionMapping dmap;

    //! Time variable assoc with new data.
    TimeInterval new_time;

    //! Time variable assoc with old data.
    TimeInterval old_time;

    //! Pointer to new-time data.
    std::unique_ptr<MultiFab> new_data;

    //! Pointer to previous time data.
    std::unique_ptr<MultiFab> old_data;

    //! Arena we should use for allocating the data.
    Arena* arena{nullptr};

    /**
    * \brief This is used as a temporary collection of FabArray header
    * names written during a checkpoint
    */
    static Vector<std::string> fabArrayHeaderNames;

    //! This is used to store preread FabArray headers
    static std::map<std::string, Vector<char> > *faHeaderMap;  // ---- [faheader name, the header]

    //! Internal helper that reads FAB headers/data from a checkpoint stream.
    void restartDoit (std::istream& is, const std::string& chkfile);
};

class StateDataPhysBCFunct
{
public:

    /**
     * \brief Helper that applies physical BCs for a particular StateData component range.
     *
     * \param sd     StateData whose data will be modified.
     * \param sc     Starting component index within \p sd.
     * \param geom_  Geometry describing the level.
     */
    StateDataPhysBCFunct (StateData& sd, int sc, const Geometry& geom_);

    /**
     * \brief Apply BCs to \p mf for components [\p dest_comp, dest_comp+\p num_comp).
     *
     * \param mf        MultiFab to modify.
     * \param dest_comp Destination component offset.
     * \param num_comp  Number of components to update.
     * \param nghost    How many ghost cells should be filled.
     * \param time      Physical time at which BCs are evaluated.
     * \param bccomp    Component offset into the descriptor's BC array.
     */
    void operator() (MultiFab& mf, int dest_comp, int num_comp, IntVect const& nghost,
                     Real time, int bccomp);

    void FillBoundary (MultiFab& mf, int dest_comp, int num_comp, IntVect const& nghost,
                       Real time, int bccomp) {
        this->operator()(mf,dest_comp,num_comp,nghost,time,bccomp); // NOLINT(readability-suspicious-call-argument)
    }
private:
    StateData* statedata;
    int src_comp;
    const Geometry* geom;

};

}

#endif /*_StateData_H_*/
