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

    int RTS::RTS_ProcCount(){
	return 1;
    }

    int RTS::RTS_MyProc(){
	return 0;
    }

    int RTS::RTS_WorkerThreadCount(){
	return 1;
    }

    int RTS::RTS_MyWorkerThread(){
	return 0;
    }

    void RTS::RTS_Init(){
    }

    void RTS::RTS_Init(int *rank, int *nProcs){
	*rank=0;
	*nProcs=1;
    }

    void RTS::RTS_Finalize(){
	//Now, no task should be alive. Thus, this routine check the content of all task queues.
	assert(WaitingQueue.size()==0);
	assert(DataFetchingQueue.size()==0);
	assert(ReadyQueue.size()==0);
	assert(RunningQueue.size()==0);
    }

    void RTS::RTS_Run(void* taskgraph, bool keepInitialGraph=false){
	AbstractTaskGraph<Task>* graph= (AbstractTaskGraph<Task>*)taskgraph;
	//visit all initial tasks 
	{
	    Task* t= graph->Begin();
	    while(t != graph->End()){
		if(graph->GetRunningMode()== _Push)
		{
		    if(t->Dependency()){//all data have arrived
			ReadyQueue.push(t);
		    }else{
			WaitingQueue.push(t);
		    }
		}else{//Pull mode
		    DataFetchingQueue.push(t);
		}
		t = graph->Next();
	    }
	}
#if 1
	bool keepRunning=true;
	while (keepRunning){
	    //Handle communication
	    {
		if(graph->GetRunningMode()== _Push)
		{
#if 0
		    //Process all messages
		    while(MsgQueue.size()){
			Data* msg= MsgQueue.front();
			MsgQueue.pop();
			TaskName name= msg->GetRecipient();
			Task* t= graph->LocateTask(name);
			t->Pull(msg);
		    }
#endif
		}else{
		    while(DataFetchingQueue.size()){
			Task* t= DataFetchingQueue.front();
			DataFetchingQueue.pop();
			t->Dependency();//send active messages to pull data from source tasks    
		    }
		}
	    }
	    //visit waiting tasks (only in push mode)
	    if(graph->GetRunningMode()== _Push)
	    { //no else
		int nWaitingTasks= WaitingQueue.size();
		for(int i=0; i<nWaitingTasks; i++){
		    Task* t= WaitingQueue.front();
		    WaitingQueue.pop();
		    if(t->Dependency()){ 
			ReadyQueue.push(t);
		    }else{
			WaitingQueue.push(t);
		    }
		}
	    }
	    //Execute ready tasks
	    {
		while(ReadyQueue.size()){
		    Task* t= ReadyQueue.front();
		    ReadyQueue.pop();
		    t->RunJob(); 
		    t->RunPostCompletion(); 
		    //Flush all outputs
		    std::queue<Data*> taskOutputs= t->GetOutputs();
		    while(taskOutputs.size()>0){
			Data* outdata= taskOutputs.front();
			taskOutputs.pop();
			TaskName dst= outdata->GetRecipient();
			int tag= outdata->GetTag();
			graph->LocateTask(dst)->GetInputs().push_back(t->MyName(), outdata, tag);
		    }
		}
	    }
	    keepRunning= WaitingQueue.size()>0 || DataFetchingQueue.size()>0|| ReadyQueue.size()>0|| RunningQueue.size()>0|| MsgQueue.size()>0;
	}
	if(!keepInitialGraph){
	    //graph->clear();
	}
#endif
    }

    void Barrier(){
	//Since we run the graph in sequential mode, this routine can be just a noop
    }

}//end namespace

