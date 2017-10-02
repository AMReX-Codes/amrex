#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-21-2017

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

namespace amrex{
    typedef std::queue<Task*> _TaskQueue;
    typedef std::queue<Data*> _MessageQueue;
    _TaskQueue _WaitingQueue;
    _TaskQueue _DataFetchingQueue;  //used in Pull model
    _TaskQueue _ReadyQueue;
    _TaskQueue _RunningQueue;
    _MessageQueue _MsgQueue;

    int RTS::ProcCount(){
	return 1;
    }

    int RTS::MyProc(){
	return 0;
    }

    int RTS::WorkerThreadCount(){
	return 1;
    }

    int RTS::MyWorkerThread(){
	return 0;
    }

    void RTS::Init(){ 
    }

    void RTS::Init(int rank, int nProcs){
	_rank=0;
	_nProcs=1;
    }

    void RTS::Finalize(){
	//Now, no task should be alive. Thus, this routine check the content of all task queues.
	assert(_WaitingQueue.size()==0);
	assert(_DataFetchingQueue.size()==0);
	assert(_ReadyQueue.size()==0);
	assert(_RunningQueue.size()==0);
    }

    void RTS::Iterate(void* taskgraph){
	AbstractTaskGraph<Task>* graph= (AbstractTaskGraph<Task>*)taskgraph;
	//visit all initial tasks 
	{
	    Task* t= graph->Begin();
	    while(t != graph->End()){
		if(graph->GetRunningMode()== _Push)
		{
		    if(t->Dependency()){//all data have arrived
			_ReadyQueue.push(t);
		    }else{
			_WaitingQueue.push(t);
		    }
		}else{//Pull mode
		    _DataFetchingQueue.push(t);
		}
		t = graph->Next();
	    }
	}
	bool keepRunning=true;
	while (keepRunning){
	    //Handle communication
	    {
		if(graph->GetRunningMode()== _Push)
		{
		    //Process messages
		    int nMsgs= _MsgQueue.size();
		    for(int i=0; i<nMsgs; i++){
			Data* msg= _MsgQueue.front();
			_MsgQueue.pop();
			TaskName name= msg->GetRecipient();
			if(graph->LocateTask(name)){
			    Task* t= graph->LocateTask(name);
			    t->GetInputs().push_back(msg->GetSource(), msg, msg->GetTag());
			}
			else _MsgQueue.push(msg); //Recipient has not been created
		    }
		}else{
		    while(_DataFetchingQueue.size()){
			Task* t= _DataFetchingQueue.front();
			_DataFetchingQueue.pop();
			t->Dependency();//send active messages to pull data from source tasks    
		    }
		}
	    }
	    //visit waiting tasks (only in push mode)
	    if(graph->GetRunningMode()== _Push)
	    { //no else
		int nWaitingTasks= _WaitingQueue.size();
		for(int i=0; i<nWaitingTasks; i++){
		    Task* t= _WaitingQueue.front();
		    _WaitingQueue.pop();
		    if(t->Dependency()){ 
			_ReadyQueue.push(t);
		    }else{
			_WaitingQueue.push(t);
		    }
		}
	    }
	    //Execute ready tasks
	    {
		while(_ReadyQueue.size()){
		    Task* t= _ReadyQueue.front();
		    _ReadyQueue.pop();
		    t->RunJob(); 
		    t->RunPostCompletion(); 
		    //Flush all outputs
		    while(t->GetOutputs().size()>0){
			Data* outdata= t->GetOutputs().front();
			t->GetOutputs().pop();
			TaskName dst= outdata->GetRecipient();
			int tag= outdata->GetTag();
			if(graph->LocateTask(dst)){
			    graph->LocateTask(dst)->GetInputs().push_back(outdata->GetSource(), outdata, tag);
			}else _MsgQueue.push(outdata); 
		    }
		    //process newly created tasks
		    while(t->GetNewTasks().size()>0){
			Task* nt= t->GetNewTasks().front();
			t->GetNewTasks().pop();
			graph->GetTaskPool()[nt->MyName()]=nt;
			if(nt->Dependency()){//all data have arrived
			    _ReadyQueue.push(nt);
			}else{
			    _WaitingQueue.push(nt);
			}
		    }
		    //keep or destroy task
		    if(t->isPersistent()){
			if(t->Dependency()){
			    _ReadyQueue.push(t);
			}else{
			    _WaitingQueue.push(t);
			}
		    }else{
			//remove task from the task pool and delete it
			graph->DestroyTask(t);
		    }
		}
	    }
	    keepRunning= _WaitingQueue.size()>0 || _DataFetchingQueue.size()>0|| _ReadyQueue.size()>0|| _RunningQueue.size()>0|| _MsgQueue.size()>0 || graph->GetTaskPool().size()>0;
	}
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

