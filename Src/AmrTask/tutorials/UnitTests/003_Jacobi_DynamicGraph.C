#include <iostream>
#include <string.h>

#include "AMReX_AbstractTask.H"
#include "AMReX_AbstractTaskGraph.H"
#include "RTS.H"

using namespace amrex;


class TokenRingTask :public Task{
    private:
	int _token;
    public:
	void Job(){
	    int taskID= MyName().GetID(0);
	    int leftNeighbor= taskID>0?taskID-1:-1; 
	    int rightNeighbor= taskID< _nTasks?taskID+1:-1; 
	    if(taskID==0){
		_token=1;
	    }else{
		if(leftNeighbor!=-1){
		    Pull(leftNeighbor, (char*)&_token);
		}
	    }
	    if(rightNeighbor!=-1){
		Push(rightNeighbor, (char*)&_token);
	    }
	    else cout<< "Task" << taskID << "has received the token with value" << _token;
	}
	void dependency(){
	    int taskID= MyName().GetID(0);
	    int leftNeighbor= taskID>0?taskID-1:-1; 
	    if(leftNeighbor!=-1) Depend_on(leftNeighbor);
	}
	void PostCompletion(){
	    //Do nothing. This task will be destroyed and no further action will be made.
	}
	static size_t _nTasks;
};


int main(int argc,char *argv[])
{
    int argCount = 0;
    int t, verbose=0;
    int rank, nProcs;
    /* Argument list
       -t: number of tasks
       -v: print out task graph information
    */ 
    while(++argCount <argc) {
	if(!strcmp(argv[argCount], "-t")) t = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-v")) verbose = atoi(argv[++argCount]);
    }
    RTS rts;
    rts.RTS_Init(&rank, &nProcs);
    TokenRingTask::_nTasks= t;
    string graphName= "TokenRing";
    if(verbose && rank==0){
	printf("Creating a 1D Token Ring Graph with %d tasks\n", t);
	printf("Running the graph with %d processes\n", rts.RTS_ProcCount());
    }
    ArrayGraph<TokenRingTask> *TokenRingGraph= new ArrayGraph<TokenRingTask>(graphName, t);
    rts.RTS_Run(TokenRingGraph, false);
    if(verbose && rank==0) printf("Graph execution completed\n");
    rts.RTS_Finalize();
};
