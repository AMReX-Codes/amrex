#ifndef _AmrTask_H_
#define _AmrTask_H_ 

#include <fstream>
#include <memory>
#include <list>

#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
#include <AMReX_FabArrayBase.H>
#include <AMReX_FabArray.H>
#include <AMReX_IntVect.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>
#include <AMReX_AmrCore.H>
bool verbose=false;

namespace amrex {
    // These are all ParmParse'd in.  Set defaults in Initialize()!!!
    extern int  plot_nfiles;
    extern int  mffile_nstreams;
    extern int  probinit_natonce;
    extern bool plot_files_output;
    extern int  checkpoint_nfiles;
    extern int  regrid_on_restart;
    extern int  use_efficient_regrid;
    extern int  plotfile_on_restart;
    extern int  checkpoint_on_restart;
    extern bool checkpoint_files_output;
    extern int  compute_new_dt_on_regrid;
    extern bool precreateDirectories;
    extern bool prereadFAHeaders;

    class AmrLevelTask;
    enum task_states{
	_coarseTimeStep=0,
	_post_coarseTimeStep,
	_timeStep,
	_timeStep_advance,
	_post_timeStep,
	_final,
    };

    /**** State transition diagrams of tasks on the coarse level (lower part) and fine level (upper part)
     *
     *							   	     _timeStep -> _timeStep_advance  -> _post_timestep	
     *							  	      ^\            ^\              |                 |   \\
     * 									\             \_____________|                 |    \\
     *								    //	 \                                            |	    \\ (signal)
     *							    	   // 	  \____________(subcycling)___________________|      \\
     *								  //							      \\
     *						        	 // 	                     			                            			        fine task
     *			-------------------------------------------------------------------------------------------------------------------------------------------------------------
     *							       // (signal)	
     *	_coarseTimeStep -> _timeStep  -> _timeStep_advance     ---> _post_timestep    			----------------------->          _post_coarseTimeStep 		coarsest task
     *	 |^\		  	        ^\		  |                      				  					      |
     *	 |  \		     	          \_______________|		         				  					      |
     *	 |   \ 		     						         				  					      |
     *	 |    \	              											  					      |
     *	 |     \									 				  				      |
     * _final   \_____________________________________________________________________________________________________________________________________________|
     */

    class AmrTask: public Task
    {
	LocalConnection l_con;
	RemoteConnection r_con;
	struct TileArray
	{
	    Vector<int> numLocalTiles;
	    Vector<int> localTileIndexMap;
	    Vector<Box> tileArray;
	} ta;
	bool _do_tiling;
	int current_tile;
	Amr* _amr;
	int max_step;
	Real stop_time;

	public:
	AmrTask(){}
	void InitAmrTask (const void* amr, int m_step, Real st_time){
	    assert(_amr);
	    TaskName name= MyName();
	    _taskState= name[0]==0?_coarseTimeStep:_timeStep;
	    _amr= (Amr*)amr;
	    max_step= m_step;
	    stop_time= st_time;
	    cumtime= _amr->cumtime;
	    dt_level.resize(_amr->finest_level+1);
	    dt_min.resize(_amr->finest_level+1);
	    level_steps.resize(_amr->finest_level+1);
	    level_count.resize(_amr->finest_level+1);
	    for (int i = 0; i <= _amr->finest_level; i++){
		dt_level[i] = _amr->dtLevel(i);
		dt_min[i] = _amr->dt_min[i];
		level_steps[i] = _amr->level_steps[i];
		level_count[i] = _amr->level_count[i];
	    }
	    which_level_being_advanced= _amr->which_level_being_advanced;
	    sub_cycle= _amr->sub_cycle;
	    subcycling_iteration =0;
            parent_subcycling_iteration=0;
	    if(name[0]<_amr->finest_level){
		MultiFab& mfSrc = _amr->amr_level[name[0]]->get_new_data(0);
		MultiFab& mfDst = _amr->amr_level[name[0]+1]->get_new_data(0);
		const BoxArray& fine_BA = mfDst.boxArray();
		BoxArray crse_S_fine_BA = fine_BA;
		crse_S_fine_BA.coarsen(_amr->refRatio(1));
		MultiFab *crse_S_fine = new MultiFab(crse_S_fine_BA, mfDst.DistributionMap(), mfDst.nComp(),0);
		TheCPC_sendup= (amrex::FabArrayBase::CPC*)&(crse_S_fine->getCPC(0, mfSrc, 0, Periodicity::NonPeriodic()));
	    }
	    if(name[0]>0){
		MultiFab& mfSrc1 = _amr->amr_level[name[0]-1]->get_new_data(0);
		MultiFab& mfDst1 = _amr->amr_level[name[0]]->get_new_data(0);

		const BoxArray& fine_BA = mfDst1.boxArray();
		BoxArray crse_S_fine_BA = fine_BA;
		crse_S_fine_BA.coarsen(_amr->refRatio(1));
		MultiFab *crse_S_fine = new MultiFab(crse_S_fine_BA, mfDst1.DistributionMap(), mfDst1.nComp(),0);
		TheCPC_pullup= (amrex::FabArrayBase::CPC*)&(crse_S_fine->getCPC(0, mfSrc1, 0, Periodicity::NonPeriodic()));
	    }
	}

	void Signal(TaskName dest, int tag){
	    Data* data= new Data(_id, dest, 1);
	    data->SetTag(tag);
	    _outputs.push(data);
	}

	bool DependSignal(TaskName src, int tag){
	    bool arrived= Depend_on(src, tag);
	    if(!arrived) return false;
	    return true;
	}
	void DependSignal_reset(TaskName src, int tag){
	    Data* msg= _neighbors_in.pop_front(src, tag);
	    msg->Free();
	}

	void Job(){
	    TaskName name= MyName();
	    if(name[0]==0){//tasks on the coarsest level
		switch(_taskState){
		    case _coarseTimeStep: 
			if ( _amr->okToContinue() &&
				(levelSteps(0) < max_step || max_step < 0) &&
				(cumtime < stop_time || stop_time < 0.0) )

			{
			    coarseTimeStepTask_init(stop_time);
			    _taskState= _timeStep;
			}else{ 
			    //_final:
			    coarseTimeStepTask_end(stop_time);
			    SelfDestroy();
			    break;
			}
		    case _timeStep: 
			timeStepTask_init(0,cumtime,1,1,stop_time);
			_taskState= _timeStep_advance;

			//advanceObj= new advanceTask(cumtime, dt_level[0], 1, 1);
		    case _timeStep_advance:
			dt_new= dt_level[0];//advanceObj->Job();
			//if(advance->state()!=advance->finalState()){
			//    delete advanceObj;
			//    break; //keep doing advance
			//}
			//else continue the main control flow
			if (_amr->finest_level >0)
			{
			    const int nloc_cpAsc = TheCPC_sendup->m_LocTags->size();
			    const int nsnds_cpAsc = TheCPC_sendup->m_SndTags->size();
			    const int nrcvs_cpAsc = TheCPC_sendup->m_RcvTags->size();
			    for(int i=0; i<nloc_cpAsc; i++)
			    {
				const amrex::FabArrayBase::CopyComTag& tag = (*TheCPC_sendup->m_LocTags)[i];
				if(name[2] == tag.srcIndex){
				    Signal(TaskName(name[0]+1, name[1], tag.dstIndex), subcycling_iteration);
				}
			    }
			}
			_taskState= _post_timeStep;
			break;
		    case _post_timeStep: 
			dt_level[0] = dt_new;
			level_steps[0]++;
			level_count[0]++;
			post_timestepTask(1);
			_taskState= _post_coarseTimeStep;
			break;
		    case _post_coarseTimeStep: 
			//_amr->amr_level[0]->postCoarseTimeStep(cumtime);
			cumtime += dt_level[0];
			_taskState= _coarseTimeStep;
			break;
		    default:
			cout<<"Error state"<<endl;
			exit(0);
		}
	    }else if(name[0]<=_amr->finest_level){
		int lev= name[0];
		switch(_taskState){
		    case _timeStep:
			timeStepTask_init(lev, cumtime+ (subcycling_iteration-1)*dt_level[lev], subcycling_iteration, _amr->n_cycle[lev], stop_time);
			_taskState= _timeStep_advance;
			//advanceObj= new advanceTask(cumtime+ (subcycling_iteration-1)*dt_level[lev], dt_level[lev], subcycling_iteration, _amr->n_cycle[lev]);
		    case _timeStep_advance:
			dt_new= dt_level[lev];//advanceObj->Job();
			//if(advance->state()!=advance->finalState()){
			//    delete advanceObj;
			//    break; //keep doing advance
			//}
			if (lev < _amr->finest_level)
			{
			    const int nloc_cpAsc = TheCPC_sendup->m_LocTags->size();
			    const int nsnds_cpAsc = TheCPC_sendup->m_SndTags->size();
			    const int nrcvs_cpAsc = TheCPC_sendup->m_RcvTags->size();
			    for(int i=0; i<nloc_cpAsc; i++)
			    {
				const amrex::FabArrayBase::CopyComTag& tag = (*TheCPC_sendup->m_LocTags)[i];
				if(name[2] == tag.srcIndex){
				    Signal(TaskName(name[0]+1, name[1], tag.dstIndex), subcycling_iteration);
				}
			    }
			}
			_taskState= _post_timeStep;
			break;
		    case _post_timeStep:
			cumtime += dt_level[lev];
			subcycling_iteration++;
			if(subcycling_iteration < _amr->n_cycle[lev]){
			    _taskState= _timeStep;
			}else{
			    parent_subcycling_iteration++;
			    subcycling_iteration=0;
			    if(parent_subcycling_iteration== _amr->n_cycle[lev-1]){
			         SelfDestroy();
			    }else _taskState= _timeStep;
			}
			//_amr->amr_level[lev]->postCoarseTimeStep(cumtime);
			break;
		    default:
			cout<<"Error state"<<endl;
			exit(0);
		}
	    } 
	}

	/*
	   if(name[0]>0 && name[0]<=_amr->finest_level){ //execute advance
	   TaskName parentName(name);
	   int tag=0;
	   Data* msg= _neighbors_in.pop_front(TaskName(parentName), tag);
	   char* pos;
	   Real time= ((Real*)msg)[0];
	   Real dt= ((Real*)msg)[1];
	   pos= sizeof(Real)*2; 
	   int  iteration= ((int*)pos)[0];
	   int  ncycle= ((int*)pos)[1];
	   advanceTask(time, dt, iteration, ncycle);
	   cout<<"task "<<name[0]<<" "<<name[1]<<" "<<name[2]<<" subcycling"<<endl;
	   SelfDestroy();//just for now
	   }
	 */

	bool Dependency(){
	    TaskName name= MyName();
	    switch(_taskState){
		case _coarseTimeStep:
		    assert(name[0]==0);
		    return true; 
		case _timeStep:
		    if(name[0]==0) return true; //at this state, coarse tasks just need to perform the initialization of the timestep
		    else { //waiting for signals from tasks at the coarser level
		        if(subcycling_iteration>0) return true;
			const int nloc_cpAsc = TheCPC_pullup->m_LocTags->size();
			const int nsnds_cpAsc = TheCPC_pullup->m_SndTags->size();
			const int nrcvs_cpAsc = TheCPC_pullup->m_RcvTags->size();
			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const amrex::FabArrayBase::CopyComTag& tag = (*TheCPC_pullup->m_LocTags)[i];
			    if(name[2] == tag.dstIndex){
				if(DependSignal(TaskName(name[0]-1, name[1], tag.srcIndex), parent_subcycling_iteration)==false) return false;
			    }
			}

			//now all signals have arrived
			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const amrex::FabArrayBase::CopyComTag& tag = (*TheCPC_pullup->m_LocTags)[i];
			    if(name[2] == tag.dstIndex){
				DependSignal_reset(TaskName(name[0]-1, name[1], tag.srcIndex), parent_subcycling_iteration);
			    }
			}
			return true;
		    }
		case _timeStep_advance:
		    return true;//advanceObj->Dependency();
		case _post_timeStep:
		    if(name[0]==0) return true; //just for now
		    else{
			return true; //just for now
		    }
		case _post_coarseTimeStep:
		    return true;
	    }
	}

	void PostCompletion(){
	}

	void SetFab(FArrayBox* fab){_fab= fab;}

	void SetMF(const FabArray<FArrayBox> &mf){
	    _mf= (FabArray<FArrayBox>*)&mf;
	}

	void SetIdx(int idx){
	    _idx=idx;
	}

	void SetLocalIdx(int lIdx){_lIdx= lIdx;}

	void enable_tiling(){_do_tiling=true;}

	Box validbox() const{
	    return _mf->box(_idx);
	}

	void generateTileArray(const IntVect& tileSize){
	    IntVect nt_in_fab, tsize, nleft;
	    int ntiles = 1;
	    const Box& bx = validbox();

	    for (int d=0; d<BL_SPACEDIM; d++) {
		int ncells = bx.length(d);
		nt_in_fab[d] = std::max(ncells/tileSize[d], 1);
		tsize    [d] = ncells/nt_in_fab[d];
		nleft    [d] = ncells - nt_in_fab[d]*tsize[d];
		ntiles *= nt_in_fab[d];
	    }

	    IntVect small, big, ijk;  // note that the initial values are all zero.
	    ijk[0] = -1;
	    for (int t = 0; t < ntiles; ++t) {
		ta.localTileIndexMap.push_back(t);
		ta.numLocalTiles.push_back(ntiles);

		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d]<nt_in_fab[d]-1) {
			ijk[d]++;
			break;
		    } else {
			ijk[d] = 0;
		    }
		}

		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d] < nleft[d]) {
			small[d] = ijk[d]*(tsize[d]+1);
			big[d] = small[d] + tsize[d];
		    } else {
			small[d] = ijk[d]*tsize[d] + nleft[d];
			big[d] = small[d] + tsize[d] - 1;
		    }
		}

		Box tbx(small, big, IndexType::TheCellType());
		tbx.shift(bx.smallEnd());

		ta.tileArray.push_back(tbx);
	    }
	}

	FabArray<FArrayBox> *_mf;
	FArrayBox* _fab;
	int _idx;
	int _lIdx;

	//! Set the timestep on each level.
	void setDtLevelTask (const Vector<Real>& dt_lev){
	    for (int i = 0; i <= _amr->finest_level; i++)
		dt_level[i] = dt_lev[i];
	}

	//! Set the timestep at one level.
	void setDtLevelTask (Real dt, int lev){
	    dt_level[lev] = dt;
	}

	//! Set the dtmin on each level.
	void setDtMinTask (const Vector<Real>& dt_min_in){
	    for (int i = 0; i <= _amr->finest_level; i++)
		dt_min[i] = dt_min_in[i];
	}

	/**  
	 * \brief What is "level" in Amr::timeStep?  This is only relevant if we are still in Amr::timeStep;
	 *      it is set back to -1 on leaving Amr::timeStep.
	 */  
	int level_being_advancedTask () const { return which_level_being_advanced; }

	int levelSteps (int lev) const { return level_steps[lev]; }

	//! Physical time.
	Real cumTimeTask () const { return cumtime; }

	void setCumTimeTask (Real t) {cumtime = t;}

	//! Time step at specified level.
	Real dtLevelTask (int level) const { return dt_level[level]; }

	//! Max time step (typically based on physics) at specified level
	Real dtMinTask (int level) const { return dt_min[level]; }

	//! Array of time steps at all levels.
	const Vector<Real>& dtLevelTask () const { return dt_level; }

	void coarseTimeStepTask_init (Real stop_time){
	    if (levelSteps(0) > 0)
	    {
		int post_regrid_flag = 0;
		_amr->amr_level[0]->computeNewDt(_amr->finest_level,
			sub_cycle,
			_amr->n_cycle,
			_amr->ref_ratio,
			_amr->dt_min,
			dt_level,
			stop_time,
			post_regrid_flag);
	    }
	    else
	    {
		_amr->amr_level[0]->computeInitialDt(_amr->finest_level,
			sub_cycle,
			_amr->n_cycle,
			_amr->ref_ratio,
			dt_level,
			stop_time);
	    }
	}

	//postCoarseTimeStepTask(cumtime);

	void coarseTimeStepTask_end(Real stop_time){
	    if (verbose > 0)
	    {
		if(isMasterTask()){
		    amrex::Print()
			<< "\nSTEP = " << level_steps[0]
			<< " TIME = "  << cumtime
			<< " DT = "    << dt_level[0] << "\n\n";
		}
	    }
	    if (_amr->record_run_info && ParallelDescriptor::IOProcessor())
	    {
		if(isMasterTask()){
		    _amr->runlog << "STEP = "  << level_steps[0]
			<< " TIME = " << cumtime
			<< " DT = "   << dt_level[0] << '\n';
		}
	    }
	    if (_amr->record_run_info_terse && ParallelDescriptor::IOProcessor())
		if(isMasterTask())
		    _amr->runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';

	    int check_test = 0;

	    if (_amr->check_per > 0.0)
	    {
		const int num_per_old = (cumtime-_amr->dt_level[0]) / _amr->check_per;
		const int num_per_new = (cumtime            ) / _amr->check_per;

		if (num_per_old != num_per_new)
		{
		    check_test = 1;
		}
	    }

	    int to_stop       = 0;
	    int to_checkpoint = 0;
	    int to_plot       = 0;
	    if (_amr->message_int > 0 && level_steps[0] % _amr->message_int == 0) {
		if(isMasterTask()){
		    if (ParallelDescriptor::IOProcessor())
		    {
			FILE *fp;
			if ((fp=fopen("dump_and_continue","r")) != 0)
			{
			    remove("dump_and_continue");
			    to_checkpoint = 1;
			    fclose(fp);
			}
			else if ((fp=fopen("stop_run","r")) != 0)
			{
			    remove("stop_run");
			    to_stop = 1;
			    fclose(fp);
			}
			else if ((fp=fopen("dump_and_stop","r")) != 0)
			{
			    remove("dump_and_stop");
			    to_checkpoint = 1;
			    to_stop = 1;
			    fclose(fp);
			}

			if ((fp=fopen("plot_and_continue","r")) != 0)
			{
			    remove("plot_and_continue");
			    to_plot = 1;
			    fclose(fp);
			}
		    }
		    int packed_data[2];
		    packed_data[0] = to_stop;
		    packed_data[1] = to_checkpoint;
		    ParallelDescriptor::Bcast(packed_data, 2, ParallelDescriptor::IOProcessorNumber());
		    to_stop = packed_data[0];
		    to_checkpoint = packed_data[1];
		    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
			_amr->last_checkpoint = level_steps[0];
			_amr->last_plotfile   = level_steps[0];
		    }

		    if (to_checkpoint && _amr->write_plotfile_with_checkpoint)
			to_plot = 1;

		    if ((_amr->check_int > 0 && level_steps[0] % _amr->check_int == 0) || check_test == 1
			    || to_checkpoint)
		    {
			_amr->checkPoint();
		    }

		    if (_amr->writePlotNow() || to_plot)
		    {
			_amr->writePlotFile();
		    }

		    if (_amr->writeSmallPlotNow())
		    {
			_amr->writeSmallPlotFile();
		    }
		    _amr->bUserStopRequest = to_stop;
		    if (to_stop)
		    {
			ParallelDescriptor::Barrier("Amr::coarseTimeStep::to_stop");
			if(ParallelDescriptor::IOProcessor()) {
			    if (to_checkpoint)
			    {
				std::cerr << "Stopped by user w/ checkpoint" << std::endl;
			    }
			    else
			    {
				std::cerr << "Stopped by user w/o checkpoint" << std::endl;
			    }
			}
		    }
		}
	    }
	}

	virtual Real advanceTask (Real time, Real dt, int iteration, int ncycle) = 0;

	void timeStepTask_init (int level, Real time, int iteration, int niter, Real stop_time){
	    which_level_being_advanced = level;	    
	    int lev_top = std::min(_amr->finest_level, _amr->max_level-1);

	    for (int i(level); i <= lev_top; ++i)
	    {
		const int old_finest = _amr->finest_level;

		if (_amr->okToRegrid(i))
		{
		    _amr->regrid(i,time);
		    if (compute_new_dt_on_regrid && (i == 0) )
		    {
			int post_regrid_flag = 1;
			_amr->amr_level[0]->computeNewDt(_amr->finest_level,
				sub_cycle,
				_amr->n_cycle,
				_amr->ref_ratio,
				dt_min,
				dt_level,
				stop_time,
				post_regrid_flag);
		    }
		    for (int k(i); k <= _amr->finest_level; ++k) {
			level_count[k] = 0;
		    }
		    if (old_finest < _amr->finest_level)
		    {
			// The new levels will not have valid time steps
			// and iteration counts.
			for (int k(old_finest + 1); k <= _amr->finest_level; ++k)
			{
			    dt_level[k]    = dt_level[k-1]/_amr->n_cycle[k];
			}
		    }
		}
		if (old_finest > _amr->finest_level) {
		    lev_top = std::min(_amr->finest_level, _amr->max_level - 1);
		}
	    }
	}

	//	Real dt_new = advanceTask(time,dt_level[level],iteration,niter);

#if 0
	void _postTimeStepTask(int level, Real time, int iteration, int niter, Real stop_time){
	    dt_min[level] = iteration == 1 ? dt_new : std::min(dt_min[level],dt_new);
	    level_steps[level]++;
	    level_count[level]++;
	    if (verbose > 0)
	    {
		amrex::Print() << "[Level " << level << " step " << level_steps[level] << "] "
		    << "Advanced " << _amr->amr_level[level]->countCells() << " cells\n";
	    }

	    // If the level signified that it wants a regrid after the advance has
	    // occurred, do that now.
	    if (_amr->amr_level[level]->postStepRegrid()) {

		int old_finest = _amr->finest_level;

		_amr->regrid(level, time);

		if (old_finest < _amr->finest_level)
		{
		    // The new levels will not have valid time steps.
		    for (int k = old_finest + 1; k <= _amr->finest_level; ++k)
		    {
			dt_level[k] = dt_level[k-1] / _amr->n_cycle[k];
		    }
		}
	    }

	    //signal timeStepTask at the next level
	    if (level < _amr->finest_level)
	    {
		const int lev_fine = level+1;

		if (_amr->sub_cycle)
		{
		    const int ncycle = _amr->n_cycle[lev_fine];

		    for (int i = 1; i <= ncycle; i++)
			timeStepTask(lev_fine,time+(i-1)*dt_level[lev_fine],i,ncycle,stop_time);
		}
		else
		{
		    timeStepTask(lev_fine,time,1,1,stop_time);
		}
	    }

	    //_amr->amr_level[level]->post_timestep(iteration);

	    // Set this back to negative so we know whether we are in fact in this routine
	    which_level_being_advanced = -1;
	}//end timeStep
#endif

	virtual void post_timestepTask(int  iteration)=0;
	private:
	amrex::FabArrayBase::CPC *TheCPC_sendup;
	amrex::FabArrayBase::CPC *TheCPC_pullup;
	amrex::FabArrayBase::CPC *TheCPC_senddown;
	amrex::FabArrayBase::CPC *TheCPC_pulldown;
	protected:
	int subcycling_iteration;
        int parent_subcycling_iteration;
	int _taskState;
	// The data that need to be privatized
	Real             cumtime;      // Physical time variable.
	Real 	         dt_new;
	Vector<Real>     dt_level;     // Timestep at this level.
	Vector<Real>     dt_min;
	Array<int>       level_count;
	Array<int>       level_steps;  // Number of time steps at this level.
	int              which_level_being_advanced; // Only >=0 if we are in Amr::timeStep(level,...)
	int              sub_cycle;
	//advanceTask 	 *advanceObj;
    };
}//end namespace
#endif /*_Amr_H_*/
