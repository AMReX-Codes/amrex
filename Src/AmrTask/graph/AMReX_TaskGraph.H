#ifndef AMREX_ABSTRACT_TASKGRAPH
#define AMREX_ABSTRACT_TASKGRAPH
//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-24-2017

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <sstream>
#include <cassert>
#include "AMReX_AbstractTask.H"
using namespace std;

namespace amrex{

    template <int D=1>
	class BlockMapping{
	    private:
		PointVect<D> _first, _last, _range, _size;
		size_t _linearSize;
		size_t _totalSize;

	    public:
		BlockMapping(PointVect<D> gSize, PointVect<D> myProc, PointVect<D> nProcs){
		    for(int i=0; i<D; i++){
			size_t block;
			size_t remainder;
			block= gSize[i]/nProcs[i]; 
			remainder= gSize[i]%nProcs[i];
			_first[i]= myProc[i]*block + (myProc[i]<remainder?myProc[i]:remainder);  
			_last[i]= _first[i]+block -1 + (myProc[i]<remainder?1:0);  
			_range[i]= _last[i]-_first[i] +1;
			_size[i]= _range[i];
		    }
		}

		BlockMapping(PointVect<D> gSize, int myProc, int nProcs){
		    size_t block;
		    size_t remainder;
		    size_t totalSize=gSize[0];
		    for(int i=1; i<D; i++) totalSize=totalSize*gSize[i];
		    _totalSize= totalSize;
		    block= totalSize/nProcs;
		    remainder= totalSize%nProcs;
		    int first_linear, last_linear;
		    first_linear= myProc*block + (myProc<remainder?myProc:remainder);
		    last_linear= first_linear+block -1 + (myProc<remainder?1:0);
		    _linearSize= last_linear-first_linear+1;
		    size_t subSize=totalSize;
		    for(int i=D-1; i>=0; i--){
			subSize= subSize/gSize[i];
			_last[i]=last_linear/subSize;
			_first[i]=first_linear/subSize;
			_range[i]= _last[i]-_first[i] +1;
			_size[i]= _range[i];
			last_linear= last_linear%subSize;
			first_linear= first_linear%subSize;
		    }
		}

		PointVect<D> first(){return _first;}
		PointVect<D> last(){return _last;}
		PointVect<D> range(){return _range;}
		PointVect<D> size(){return _size;}
		size_t linearSize(){return _linearSize;}
		size_t totalSize(){return _totalSize;}
	};

    enum RunningMode{
	_Push=0,//The task will not be scheduled until all data dependencies are satified. Data are pushed by other tasks without its consent.
	_Pull//The task will be fetched and it will pull dependencies from other tasks 
    };

    class SelfAssociate{
	public:
	    TaskName TaskAssociate(TaskName name){
		return name;
	    }
    };

    //!This class defines data and task spaces, as well as options to decompose these spaces
    template <class T>
	class AbstractTaskGraph{
	    struct TaskComp {
		bool operator() (const TaskName& lhs, const TaskName& rhs) const
		{
		    return lhs<rhs;
		}
	    };
	    protected:
	    int _nProcs, _rank;
	    //! Graph name, e.g. "3DStencil", "SparseLU", etc. The default value is empty string
	    string _graphName; 
	    std::vector<Task*> _initialTasks; 
	    std::map<TaskName, Task*, TaskComp> _taskPool;
	    //! A process can iterate over tasks that it owns. However, how application tasks are mapped to processes will be defined at higher level classes.
	    Task *_begin, *_end, *_current;
	    typename std::vector<Task*>::iterator _currIt;
	    size_t _nLocalTasks;
	    RunningMode _mode;

	    public:
	    AbstractTaskGraph(string graphName=""): _graphName(graphName), _rank(0), _nProcs(1){
		_begin=NULL;
		_end=NULL;
		_current=NULL;
		_nLocalTasks=0;
		_mode= _Push;
	    }
	    void DestroyGraph(){
		    for(typename std::vector<T*>::iterator it= _initialTasks.begin(); it!= _initialTasks.end(); it++){
			delete (*it);
		    }
		    _initialTasks.clear();
		    _taskPool.clear();
	    }
	    void DestroyTask(Task* t){
		_taskPool.erase(t->MyName());
	    }
	    int MyProc(){return _rank;}  
	    int ProcCount(){return _nProcs;}  
	    RunningMode GetRunningMode(){return _mode;}
	    string GetTaskName(){return _graphName;}
	    void GraphSynchronize();
	    std::map<TaskName, Task*, TaskComp> &GetTaskPool(){return _taskPool;}
	    Task* LocateTask(TaskName name){
		if(_taskPool.find(name)!= _taskPool.end())
		    return _taskPool[name];
		return NULL;
	    }
	    virtual int FindProcessAssociation(TaskName name){ //maps task name to process rank
	    }
	    //!First element stored in the process
	    Task* Begin(){
		return _begin;
	    }
	    //!Last element stored in the process
	    Task* End(){
		return _end;
	    }
	    //! The next element
	    Task* Next(){
		_currIt++;
		_current= *_currIt;
		return _current;
	    }
	    //! The current element
	    Task* Current(){
		return _current;
	    }
	};


    /**
     * \brief This class is useful when we need a task graph to partition an irregular geometry. 
     * Tasks are created from a list of names (for example, box/data tile names), which can have multiple dimensions.
     * These tasks are stored in a distributed hash map and can be accessed in constant time.
     * A task can create another task dynamically and set it location (i.e., which process owns the newly created task). 
     * The default location is that of the the parent task.
     */
    template <class T, int D=1, class A=SelfAssociate>
	class ArrayGraph: public AbstractTaskGraph<T>{
	    protected:
		string _graphName;
		BlockMapping<D> *_taskMap;
		PointVect<D> _graphSize;

	    public:
		//! Create a 1D Task Graph
		ArrayGraph(string graphName, int graphSize, int rank, int nProcs): _graphName(graphName), _graphSize(graphSize){
		    assert(D==1);
		    AbstractTaskGraph<T>::_nProcs= nProcs;
		    AbstractTaskGraph<T>::_rank= rank;
		    _taskMap= new BlockMapping<1>(PointVect<1>(graphSize), PointVect<1>(rank), PointVect<1>(nProcs));
		    for(int i=_taskMap->first()[0], idx=0; i<=_taskMap->last()[0]; i++, idx++){
			TaskName name(i);
			T *t= new T();
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
		//! Create a multidimensional graph and LINEARLY map it to processors
		ArrayGraph(string graphName, PointVect<D> graphSize, int rank, int nProcs): _graphName(graphName){
		    _graphSize= graphSize;
		    AbstractTaskGraph<T>::_nProcs= nProcs;
		    AbstractTaskGraph<T>::_rank= rank;
		    _taskMap= new BlockMapping<D>(graphSize, rank, nProcs);
		    PointVect<D> p=_taskMap->first();
		    for(int i=0; i<_taskMap->linearSize(); i++){
			TaskName name;
			name.SetSize(D);
			for(int j=0; j<D; j++){
			    name[j]= p[j];
			}
			T *t= new T();
			t->SetName(name);
			AbstractTaskGraph<T>::_initialTasks.push_back(t);
			AbstractTaskGraph<T>::_taskPool[name]= t;
			for(int d=0; d<D; d++)
			    if(p[d]+1 < graphSize[d]){
				p[d]++; 
				//reset lower dimensions
				for(int dd=0; dd<d; dd++) p[dd]=0;
				break;
			    }
		    }
		    AbstractTaskGraph<T>::_begin= *(AbstractTaskGraph<T>::_initialTasks.begin());
		    AbstractTaskGraph<T>::_end= *(AbstractTaskGraph<T>::_initialTasks.end());
		    AbstractTaskGraph<T>::_currIt= AbstractTaskGraph<T>::_initialTasks.begin();
		    AbstractTaskGraph<T>::_current= *(AbstractTaskGraph<T>::_currIt);
		    AbstractTaskGraph<T>::_mode= _Push;
		}
		int FindProcessAssociation(TaskName name){
		    A associate;
		    TaskName n= associate.TaskAssociate(name);
		    size_t block_size= _taskMap->totalSize()/AbstractTaskGraph<T>::ProcCount();
		    size_t val=n[0];
		    size_t stride=1;
		    for(int d=1; d<D; d++){
			stride *= _graphSize[d-1];
			val+= n[d]*stride;
		    }
		    return val/block_size;
		}

		BlockMapping<D>* GetTaskMap(){return &_taskMap;}
		void Destroy(){
		    for(typename std::vector<T*>::iterator it= AbstractTaskGraph<T>::_initialTasks.begin(); it!= AbstractTaskGraph<T>::_initialTasks.end(); it++){
			delete (*it);
		    }
		    AbstractTaskGraph<T>::_initialTasks.clear();
		    AbstractTaskGraph<T>::_taskPool.clear();
		    delete _taskMap;
		}
	};

#if 0
    //! This class should be used when we want to create a rectangular task space (data and/or time).
    template <class T, int D=1>
	class CartesianGraph: public AbstractTaskGraph<T>{
	    protected:
		string _graphName;
		BlockMapping<T, D>* _taskMap;

	    public:
		//create a multi-dimensional taskGraph using P processes
		CartesianGraph(string graphName="", PointVect<D> graphSize, PointVect<D> Prank, PointVect<D> Psize){
		    _taskMap= new BlockMapping<T, D>(graphSize, Prank, Psize);
		    if(_taskMap.size() !=0){ 
			TaskName name(_taskMap.first()[0]);
			size_t nTasks= _taskMap.size()[0];
			for (int d=1; d<D; d++){
			    name[d]=_taskMap.first()[d];
			    nTasks= nTasks*_taskMap.size()[d];
			}
			for (int d=0; d<D; d++){
			    for(int i=_taskMap->first()[d], idx=0; i<=_taskMap->last()[d]; i++, idx++){
				T *t= new T();
				t->SetName(name);
				AbstractTaskGraph<T>::_initialTasks.push_back(t);
				AbstractTaskGraph<T>::_taskPool[name]= t;
			    }
			}
		    }


		}
		BlockMapping<T, D>* GetTaskMap(){return &_taskMap;}
	};
#endif


#if 0

    //Data are distributed linearly, though a multi-dimensional key (e.g. task name) is required to locate them
    template <class T, int D=1>
	class LinearMapping::public AbstractMapping{
	    private:
		std::unorderedMap<PointVect<D>, std::list<T*>, PointVect::shift_hasher> _data; 
	    public:
		LinearMapping(std::vector< PointVect<D> > nameVect){
		    _begin[0]= nameVect[0];
		    _end[0]= nameVect[nameVect.size()-1];
		    _size[0]= nameVect.size();
		}
	};

    //

    template <class T, int D=1>
	class CyclicMapping: public AbstractMapping{
	    private:

	    public: 
		CyclicMapping(PointVect<D> gSize, PointVect<D> myProc, PointVect<D> nProcs){
		    for(int i=0; i<D; i++){
			size_t block;
			size_t remainder;
			block= gSize[i]/nProcs[i]; 
			remainder= gSize[i]%nProcs[i];
			_begin[i]= myProc[i]<gSize[i]?myProc[i]:gSize[i];
			_size[i]= _begin[i]<gSize[i]?block+myProc[i]<remainder?1:0; 
			_end[i]= _begin[i]+ (size[i]-1)*nProcs[i]; 
		    }
		}
	};

    template <class T, int D=1>
	class BlockyCylicMapping{
	    private:
		CyclicMapping(PointVect<D> gSize, PointVect<D> blockSize, PointVect<D> myProc, PointVect<D> nProcs){
		    for(int i=0; i<D; i++){
			size_t block;
			size_t block_block;
			size_t remainder;
			block= gSize[i]/blockSize[i]+gSize[i]%blockSize[i]==0?0:1;
			block_block= block/nProcs[i];
			remainder= block%nProcs[i];

			_begin[i]= myProc[i]<block?myProc[i]:gSize[i];
			_size[i]= _begin[i]<gSize[i]?block_block+myProc[i]<remainder?1:0; //number of blocks
			_end[i]= _begin[i]+ (size[i]-1)*nProcs[i];
		    }
		}
	};
#endif

}

#endif
