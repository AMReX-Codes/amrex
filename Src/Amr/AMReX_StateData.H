
#ifndef AMREX_StateData_H_
#define AMREX_StateData_H_

#include <memory>

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
    * \brief Constructor that properly initializes data members.
    *
    * \param p_domain
    * \param grds
    * \param dm
    * \param d
    * \param cur_time
    * \param dt
    * \param factory
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

    StateData (StateData&& rhs) noexcept;

    StateData (StateData const& rhs) = delete;
    void operator= (StateData && rhs) = delete;
    void operator= (StateData const& rhs);


    /**
    * \brief Initializes data members if you used default constructor.
    *
    * \param p_domain
    * \param grds
    * \param dm
    * \param d
    * \param cur_time
    * \param dt
    * \param factory
    */
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
		 const DistributionMapping& dm,
                 const StateDescriptor& d,
                 Real                   cur_time,
                 Real                   dt,
                 const FabFactory<FArrayBox>& factory);

    /**
    * \brief Copies old data from another StateData object and sets the same time level.
    * If old data is uninitialized, allocates it with same properties as the input data.
    *
    * \param state
    */
    void copyOld (const StateData& state);

    /**
    * \brief Copies new data from another StateData object and sets the same time level.
    * If new data is uninitialized, allocates it with the same properties as the input data.
    *
    * \param state
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
    * \brief Old data becomes new data and new time is incremented by dt.
    *
    * \param dt
    */
    void swapTimeLevels (Real dt);

    /**
    * \brief Swaps old data with a new MultiFab.
    * Deletes the previous old data.
    *
    * \param mf
    */
    void replaceOldData ( MultiFab&& mf );

    /**
    * \brief Swaps old data with another StateData.
    * Does not delete the previous old data.
    *
    * \param s
    */
    void replaceOldData ( StateData& s );

    /**
    * \brief Swaps new data with a new MultiFab.
    * Deletes the previous new data.
    *
    * \param mf
    */
    void replaceNewData ( MultiFab&& mf );

    /**
    * \brief Swaps new data with another StateData.
    * Does not delete the previous new data.
    *
    * \param s
    */
    void replaceNewData ( StateData& s );


    /**
    * \brief Sets time of old and new data.
    *
    * \param t_new
    * \param dt_old
    * \param dt_new
    */
    void setTimeLevel (Real t_new,
                       Real dt_old,
                       Real dt_new);

    /**
    * \brief Sets time of old data.
    *
    * \param t_old
    */
    void setOldTimeLevel (Real t_old);

    /**
    * \brief Sets time of new data.
    *
    * \param t_new
    */
    void setNewTimeLevel (Real t_new);

    void syncNewTimeLevel (Real t_new);

    void RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                       Vector<MultiFabId>&      mfid);

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
    * \brief Set physical bndry values
    *
    * \param dest
    * \param time
    * \param dx
    * \param prob_domain
    * \param dest_comp
    * \param src_comp
    * \param num_comp
    */
    void FillBoundary (FArrayBox&     dest,
                       Real           time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dest_comp,
                       int            src_comp,
                       int            num_comp = 1);

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
    * \param name
    * \param fullpathname
    * \param os
    * \param how
    * \param dump_old
    */
    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os,
                     VisMF::How         how,
                     bool               dump_old = true);

    /**
    * \brief Restart with domain box, grids, and dmap provided
    *
    * \param is
    * \param p_domain
    * \param grds
    * \param dm
    * \param factroy
    * \param d
    * \param restart_file
    */
    void restart (std::istream&          is,
		  const Box&             p_domain,
		  const BoxArray&        grds,
		  const DistributionMapping& dm,
                  const FabFactory<FArrayBox>& factroy,
                  const StateDescriptor& d,
                  const std::string&     restart_file);

    /**
    * \brief or from another similar state
    *
    * \param d
    * \param rhs
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

    const DistributionMapping& DistributionMap () const noexcept { return dmap; }

    /**
    *
    * \param new_dmap
    */
    void setDistributionMap ( DistributionMapping& new_dmap ) noexcept { dmap = new_dmap; }

    const FabFactory<FArrayBox>& Factory () const noexcept { return *m_factory; }

    /**
    * \brief Returns the current time.
    */
    Real curTime () const noexcept {
	return (desc->timeType() == StateDescriptor::Point) ?
	    new_time.stop : 0.5*(new_time.start + new_time.stop);
    }

    /**
    * \brief Returns the previous time.
    */
    Real prevTime () const noexcept {
	return (desc->timeType() == StateDescriptor::Point) ?
	    old_time.stop : 0.5*(old_time.start + old_time.stop);
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

    /**
    * \brief Returns the FAB of new data at grid index `i'.
    *
    * \param i
    */
    FArrayBox& newGrid (int i) noexcept { BL_ASSERT(new_data != nullptr); return (*new_data)[i]; }

    /**
    * \brief Returns the FAB of old data at grid index `i'.
    *
    * \param i
    */
    FArrayBox& oldGrid (int i) noexcept { BL_ASSERT(old_data != nullptr); return (*old_data)[i]; }

    /**
    * \brief Returns boundary conditions of specified component on the specified grid.
    *
    * \param comp
    * \param i
    */
    BCRec getBC (int comp, int i) const noexcept;

    /**
    * \brief Prints out the time interval.
    *
    * \param os
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

    void getData (Vector<MultiFab*>& data,
		  Vector<Real>& datatime,
		  Real time) const;


    /**
    * \brief These facilitate prereading FabArray headers to avoid
    * synchronization when reading multiple FabArrays
    */
    static const Vector<std::string> &FabArrayHeaderNames() { return fabArrayHeaderNames; }
    static void ClearFabArrayHeaderNames() { fabArrayHeaderNames.clear(); }

    static void SetFAHeaderMapPtr(std::map<std::string, Vector<char> > *fahmp) { faHeaderMap = fahmp; }


private:

    std::unique_ptr<FabFactory<FArrayBox> > m_factory;

    struct TimeInterval
    {
        Real start;
        Real stop;
    };

    //! Pointer to data descriptor.
    const StateDescriptor* desc;

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

    /**
    * \brief This is used as a temporary collection of FabArray header
    * names written during a checkpoint
    */
    static Vector<std::string> fabArrayHeaderNames;

    //! This is used to store preread FabArray headers
    static std::map<std::string, Vector<char> > *faHeaderMap;  // ---- [faheader name, the header]

    void restartDoit (std::istream& is, const std::string& restart_file);
};

class StateDataPhysBCFunct
    : public PhysBCFunctBase
{
public:

    StateDataPhysBCFunct (StateData& sd, int sc, const Geometry& geom_);
    
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time, int bccomp) override;
    using PhysBCFunctBase::FillBoundary;
private:
    StateData* statedata;
    int src_comp;
    const Geometry& geom;

};

}

#endif /*_StateData_H_*/
