#ifndef AMREX_ASYNC_MFITER
#define AMREX_ASYNC_MFITER
//Question? email tannguyen@lbl.gov
//Created 09-07-2017
//Last modification 09-07-2017

#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
#include <AMReX_FabArrayBase.H>
#include <AMReX_IntVect.H>
#include <AMReX_FArrayBox.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

    class Action :public Task{
	protected:
	    int _iter, _nIters;
	public:
	    Action(){
		_iter=-1;
	    }
	    void SetSteps(int nIters){
		_nIters= nIters;
	    }
	    virtual void Compute()=0;
	    void Job(){
		if(_iter<_nIters){
		    //FillBoundary_Pull();
		    //visit tiles in current Fab
		    Compute();
		    //if(_iter<nStep-1) FillBoundary_Push();
		}else SelfDestroy();
		_iter++;
	    }
	    bool Dependency(){
		if(_iter==-1) return true;
		else{
		    //use metadata to check for msgs
		} 
	    }
	    void PostCompletion(){
	    }
	    bool isValid(){
	    }
	    Box validbox(){
	    }
    };

    template <class T>
	class MFGraph: public AbstractTaskGraph<T>{
	    protected:
		string _graphName;
	    public:
		MFGraph(const FabArrayBase &mf, int nSteps, int rank, int nProcs){
		    AbstractTaskGraph<T>::_nProcs= nProcs;
		    AbstractTaskGraph<T>::_rank= rank;
		    for(int i; i<=mf.local_size(); i++){
			TaskName name(mf.IndexArray()[i]);
			T *t= new T();
			t->SetSteps(nSteps);
			t->SetName(name);
			AbstractTaskGraph<T>::_initialTasks.push_back(t);
			AbstractTaskGraph<T>::_taskPool[name]= t;
		    }
		    AbstractTaskGraph<T>::_begin= *(AbstractTaskGraph<T>::_initialTasks.begin());
		    AbstractTaskGraph<T>::_end= *(AbstractTaskGraph<T>::_initialTasks.end());
		    AbstractTaskGraph<T>::_currIt= AbstractTaskGraph<T>::_initialTasks.begin();
		    AbstractTaskGraph<T>::_current= *(AbstractTaskGraph<T>::_currIt);
		    AbstractTaskGraph<T>::_mode= _Push;
		}
		int FindProcessAssociation(TaskName name){
		    return 0;
		}
	};

    template <class A>
	class AMFIter{
	    private:
		RTS rts;
		MFGraph<A> *graph;
	    public:
		AMFIter(const FabArrayBase &mf, int nSteps){
		    int rank= ParallelDescriptor::MyProc();
		    int nProcs= ParallelDescriptor::NProcs();
		    rts.Init(rank, nProcs);
		    graph= new MFGraph<A>(mf, nSteps, rank, nProcs);
		}
		void Iterate(){
		    rts.Iterate(graph);
		    rts.Finalize();
		}
	};
} //end namespace
#endif
