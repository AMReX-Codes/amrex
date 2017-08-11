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
		    while(t->GetOutputs().size()>0){
			Data* outdata= t->GetOutputs().front();
			t->GetOutputs().pop();
			TaskName dst= outdata->GetRecipient();
			int tag= outdata->GetTag();
			graph->LocateTask(dst)->GetInputs().push_back(t->MyName(), outdata, tag);
		    }
		    //keep or destroy task
		    if(t->isPersistent()){
			if(t->Dependency()){
			    ReadyQueue.push(t);
			}else{
			    WaitingQueue.push(t);
			}
		    }else{
			//later

		    }
		}
	    }
	    keepRunning= WaitingQueue.size()>0 || DataFetchingQueue.size()>0|| ReadyQueue.size()>0|| RunningQueue.size()>0|| MsgQueue.size()>0;
	}
	if(!keepInitialGraph){
	    //graph->clear();
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

