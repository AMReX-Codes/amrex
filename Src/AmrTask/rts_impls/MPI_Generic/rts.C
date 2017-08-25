//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 08-14-2017
#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
#include <mpi.h>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#include "sysInfo.H"
#include "mylock.h"
#include <pthread.h>

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

#include "dl_malloc.h"
namespace amrex{
    //we don't use template for task and message queuese since in the future we may implement them in different ways
    class _TaskQueue {
	private:
	std::queue<Task*> _queue;
	MyLock _lock;
	public:
	void push(Task* t){
	    _lock.lock();
	    _queue.push(t);
	    _lock.unlock();
	}
	Task* pop(){
	    _lock.lock();
	    if(_queue.size()>0) {
		Task*t = _queue.front();
		_queue.pop();
	        _lock.unlock();
		return t;
            }
	    _lock.unlock();
	    return NULL;
	}
	size_t size(){ return _queue.size();}
    };

    class _MessageQueue{
	private:
	std::queue<Data*> _queue;
	MyLock _lock;
	public:
	void push(Data* &d){
	    _lock.lock();
	    _queue.push(d);
	    _lock.unlock();
	}
        Data* pop(){
            _lock.lock();
            if(_queue.size()>0) {
                Data*d = _queue.front();
                _queue.pop();
	        _lock.unlock();
                return d;
            }
            _lock.unlock();
            return NULL;
        }
	size_t size(){ return _queue.size();}
    };

    struct RtsDomain{
	_TaskQueue _WaitingQueue; 
	_TaskQueue _DataFetchingQueue;  //used in Pull model
	_TaskQueue _ReadyQueue; 
	_TaskQueue _RunningQueue; 
	_TaskQueue _ToCreateTaskQueue; 
	_TaskQueue _ToDestroyTaskQueue; 
	_MessageQueue _MsgQueue; 
	pthread_t *_threads;
	bool _doneSignal;
        int _size;
        RtsDomain(){_threads=NULL; _doneSignal=false; _size=0;};
    };
    mspace *mspaces;
    void ** segment_ptrs;
    int numa_nodes;
    RtsDomain *dom;
    AbstractTaskGraph<Task>* graph;


    void NUMA_malloc(int numaID, void** ptr, size_t size){
	*ptr = mspace_memalign(mspaces[numaID], 64, size);
    }
    void NUMA_free(int numaID, void* ptr){
	mspace_free(mspaces[numaID], ptr);
    }


    int RTS::RTS_ProcCount(){
	int nProcs=1;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	return nProcs;
    }

    int RTS::RTS_MyProc(){
	int myRank=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	return myRank;
    }

    int RTS::RTS_WorkerThreadCount(){
	return _nWrks;
    }

    int RTS::RTS_MyWorkerThread(){
	return 0;
    }

    struct argT {
	int numaID;
	int tid;
    };
    void run(void* threadInfo){
	argT *args= (argT*)threadInfo;
	int numaID= args->numaID;
	int tid= args->tid;
        while(true){
                if(dom[numaID]._ReadyQueue.size()){
                    Task* t= dom[numaID]._ReadyQueue.pop();
                    if(t){
                      t->RunJob();
                      t->RunPostCompletion();
                      //Flush all outputs
                      while(t->GetOutputs().size()>0){
                        Data* outdata= t->GetOutputs().front();
                        t->GetOutputs().pop();
                        if(t){
                           TaskName dst= outdata->GetRecipient();
                           int tag= outdata->GetTag();
                           if(graph->LocateTask(dst)){
                             graph->LocateTask(dst)->GetInputs().push_back(outdata->GetSource(), outdata, tag);
                           }else dom[numaID]._MsgQueue.push(outdata);
                        }
                      }
                      //process newly created tasks
                      while(t->GetNewTasks().size()>0){
                        Task* nt= t->GetNewTasks().front();
                        t->GetNewTasks().pop();
                        dom[numaID]._ToCreateTaskQueue.push(nt);
                      }
                      //keep or destroy current task
                      if(t->isPersistent()){
                        if(t->Dependency()){
                            dom[numaID]._ReadyQueue.push(t);
                        }else{
                            dom[numaID]._WaitingQueue.push(t);
                        }
                      }else{
                        dom[numaID]._ToDestroyTaskQueue.push(t);
                      }
                   }
                }
	    if(dom[numaID]._doneSignal) break;
        }
	free(args);
    }

    void RTS::RTS_Init(){
	int provided;
	MPI_Init_thread(0, 0, MPI_THREAD_FUNNELED, &provided);
	if(provided == MPI_THREAD_SINGLE){//with this MPI, process can't spawn threads
	    if(_nWrks>1){
		cerr << "Spawning threads is not allowed by the MPI implementation" << std::endl;;
	    }    
	}
	NodeHardware hw = query_node_hardware();
	assert(_nWrks>0 && _nWrks <= hw.core_per_numa * hw.numa_per_node);

	bool numaAware=true;
	char* env= getenv("ENABLE_NUMA_AWARE");
	numaAware= (env!=NULL);
	if(numaAware){ //the process covers multiple NUMA nodes
	    numa_nodes= hw.numa_per_node;
	    mspaces= new mspace[numa_nodes];
	    segment_ptrs= new void*[numa_nodes];
	    int worker_per_numa = _nWrks / numa_nodes;
	    int remainder= _nWrks % numa_nodes;
	    int r=0;
	    int base=0; 
	    int localID=-1;
	    //create a list of persistent threads for each NUMA node
	    cpu_set_t cpuset;
	    pthread_attr_t attr;
	    pthread_attr_init(&attr);
	    dom= new RtsDomain[numa_nodes];
	    for(int i=0; i<numa_nodes; i++)
		dom[i]._threads= new pthread_t[worker_per_numa+1];
	    for(int i=0, domNo=-1; i<_nWrks; i++){
		localID++;
		if(localID==0){
		    domNo++;
		    segment_ptrs[domNo] = new char[SEGMENT_SIZE];
		    mspaces[domNo] = create_mspace_with_base(segment_ptrs[domNo], SEGMENT_SIZE, 1);
		    mspace_set_footprint_limit(mspaces[domNo], SEGMENT_SIZE);
		}
		CPU_ZERO(&cpuset);
		CPU_SET(base+localID, &cpuset);
                if(! (localID==0 && domNo==0)){
		  pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
		  argT* arg= new argT;
		  arg->numaID= domNo;
		  arg->tid= localID;
		  int err = pthread_create(&(dom[domNo]._threads[localID]), &attr, (void*(*)(void*))run, arg);
                }else dom[domNo]._threads[localID]= pthread_self();// master thread
                dom[domNo]._size++;
		if(r<remainder && localID == worker_per_numa){
		    localID=-1;
		    base+= hw.core_per_numa;
		    r++;
		}else if(r==remainder && localID == (worker_per_numa-1)){
		    localID=-1;
		    base+= hw.core_per_numa;
		}
	    }
	}else{
	    numa_nodes= 1;
	    dom= new RtsDomain[1];
	    mspaces= new mspace[1];
	    segment_ptrs= new void*[1];
	    segment_ptrs[0]= new char[SEGMENT_SIZE];
	    mspaces[0] = create_mspace_with_base(segment_ptrs[0], SEGMENT_SIZE, 1);
	    mspace_set_footprint_limit(mspaces[0], SEGMENT_SIZE);

	    //create a single list of persistent threads and set the thread affinity 
	    cpu_set_t cpuset;
	    cpu_set_t mycpuset;
	    pthread_attr_t attr;
	    pthread_attr_init(&attr);
	    dom[0]._threads= new pthread_t[_nWrks];
	    pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	    for(int i=0, j=0; i<CPU_SETSIZE && j<_nWrks; i++) {
		if(CPU_ISSET(i, &cpuset)){ 
		    CPU_ZERO(&mycpuset);
		    CPU_SET(i, &mycpuset);
                    if(j!=0){
		      argT* arg= new argT;
		      arg->numaID= 0;
		      arg->tid= j;
		      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &mycpuset);
		      int err = pthread_create(&(dom[0]._threads[j]), &attr, (void*(*)(void*))run, arg);
                    }else dom[0]._threads[j]= pthread_self();// master thread
                    dom[0]._size++;
		    j++;
		}
	    }
	}
    }

    void RTS::RTS_Init(int *rank, int *nProcs){
	RTS_Init();
	MPI_Comm_rank(MPI_COMM_WORLD, rank);
	MPI_Comm_size(MPI_COMM_WORLD, nProcs);
    }

    void RTS::RTS_Finalize(){
	for(int i=0; i<numa_nodes; i++)
	{
	    //Now, no task should be alive. Thus, this routine check the content of all task queues.
	    assert(dom[i]._WaitingQueue.size()==0);
	    assert(dom[i]._DataFetchingQueue.size()==0);
	    assert(dom[i]._ReadyQueue.size()==0);
	    assert(dom[i]._RunningQueue.size()==0);
	    free(dom[i]._threads);
	    free(segment_ptrs[i]);
	}
	//free(dom); //check why this causes a segfault
	free(mspaces);
	free(segment_ptrs);
    }

    void RTS::RTS_Run(void* taskgraph){
	//the master thread distributes tasks to workers
	graph= (AbstractTaskGraph<Task>*)taskgraph;
	//visit all initial tasks 
	{
	    Task* t= graph->Begin();
	    int numaID=0;
	    while(t != graph->End()){
		if(graph->GetRunningMode()== _Push)
		{
		    if(t->Dependency()){//all data have arrived
			dom[numaID]._ReadyQueue.push(t);
		    }else{
			dom[numaID]._WaitingQueue.push(t);
		    }
		}else{//Pull mode
		    dom[numaID]._DataFetchingQueue.push(t);
		}
		t = graph->Next();
		numaID= (numaID+1)%numa_nodes; //just use a simple round robin distribution for now
	    }
	}
	bool keepRunning=true;
        //allocate a static buffer for incoming messages
        size_t max_buf_size=SEGMENT_SIZE/16;
        char* env= getenv("MAX_MSG_SIZE");
        if(env) max_buf_size= atoi(env);
        char* _recvbuffer= new char[max_buf_size];
	while (keepRunning){
	    //Handle communication
	    {
		if(graph->GetRunningMode()== _Push)
		{
		    //Process outgoing messages for all domains
		    for(int d=0; d<numa_nodes; d++){
   		       int nMsgs= dom[d]._MsgQueue.size();
		       for(int i=0; i<nMsgs; i++){
		  	   Data* msg= dom[0]._MsgQueue.pop();
  			   if(msg){
			     TaskName name= msg->GetRecipient();
			     if(graph->LocateTask(name)){
			        Task* t= graph->LocateTask(name);
			        t->GetInputs().push_back(msg->GetSource(), msg, msg->GetTag());
			     }
			     else dom[0]._MsgQueue.push(msg); //Recipient has not been created
                           }
                       }
		    }
		}else{
		}
	    }
	    //visit waiting tasks in domain 0 
	    if(graph->GetRunningMode()== _Push)
	    { //no else
		int nWaitingTasks= dom[0]._WaitingQueue.size();
		for(int i=0; i<nWaitingTasks; i++){
		    Task* t= dom[0]._WaitingQueue.pop();
		    if(t->Dependency()){ 
			dom[0]._ReadyQueue.push(t);
		    }else{
			dom[0]._WaitingQueue.push(t);
		    }
		}
	    }
	    //Execute only one task at a time
	    {
		if(dom[0]._ReadyQueue.size()){
		    Task* t= dom[0]._ReadyQueue.pop();
		    if(t){
		      t->RunJob(); 
		      t->RunPostCompletion(); 
		      //Flush all outputs
		      while(t->GetOutputs().size()>0){
			Data* outdata= t->GetOutputs().front();
			t->GetOutputs().pop();
                        if(t){
			   TaskName dst= outdata->GetRecipient();
			   int tag= outdata->GetTag();
			   if(graph->LocateTask(dst)){
			     graph->LocateTask(dst)->GetInputs().push_back(outdata->GetSource(), outdata, tag);
			   }else dom[0]._MsgQueue.push(outdata); 
                        }
		      }
		      //process newly created tasks for domain 0
		      while(t->GetNewTasks().size()>0){
			Task* nt= t->GetNewTasks().front();
			t->GetNewTasks().pop();
			graph->GetTaskPool()[nt->MyName()]=nt;
			if(nt->Dependency()){//all data have arrived
			    dom[0]._ReadyQueue.push(nt);
			}else{
			    dom[0]._WaitingQueue.push(nt);
			}
		      } 
		      //keep or destroy task for domain 0
		      if(t->isPersistent()){
			if(t->Dependency()){
			    dom[0]._ReadyQueue.push(t);
			}else{
			    dom[0]._WaitingQueue.push(t);
			}
		      }else{
			//remove task from the task pool and delete it
			graph->DestroyTask(t);
		      }
                   }
		}
	    }
            //service new task creation and destroy for other workers
            for(int d=0; d<numa_nodes; d++){
		if(dom[d]._ToCreateTaskQueue.size()){
		   Task* nt= dom[d]._ToCreateTaskQueue.pop();
		   if(nt){
                        graph->GetTaskPool()[nt->MyName()]=nt;
                        if(nt->Dependency()){//all data have arrived
                            dom[d]._ReadyQueue.push(nt);
                        }else{
                            dom[d]._WaitingQueue.push(nt);
                        }
		   }
		}
	    }
            for(int d=0; d<numa_nodes; d++){
                if(dom[d]._ToDestroyTaskQueue.size()){
		   Task* ot= dom[d]._ToDestroyTaskQueue.pop();
		   if(ot){
			graph->DestroyTask(ot);
		   }
		}
            }

	    keepRunning=false;
	    for(int d=0; d<numa_nodes; d++)
       	        if(dom[d]._WaitingQueue.size()>0 || dom[d]._DataFetchingQueue.size()>0|| dom[d]._ReadyQueue.size()>0|| dom[d]._RunningQueue.size()>0|| dom[d]._MsgQueue.size()>0 || graph->GetTaskPool().size())
	            keepRunning=true;
	}
        for(int d=0; d<numa_nodes; d++) dom[d]._doneSignal=true;
        free(_recvbuffer);
    }

    const double kMicro = 1.0e-6;
    double RTS::Time()
    {
	struct timeval TV;

	const int RC = gettimeofday(&TV, NULL);
	if(RC == -1)
	{
	    printf("ERROR: Bad call to gettimeofday\n");
	    return(-1);
	}
	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
    } 

    void RTS::Barrier(){
	//nothing
    }

}//end namespace

