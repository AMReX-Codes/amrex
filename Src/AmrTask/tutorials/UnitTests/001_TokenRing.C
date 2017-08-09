#include <iostream>
#include <string.h>

#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"

using namespace amrex;


class TokenRingTask :public Task{
    private:
	int _token;
    public:
	static int _nTasks;
	void Job(){
	    int taskID= MyName()[0];
	    int leftNeighbor= taskID>0?taskID-1:-1; 
	    int rightNeighbor= taskID< _nTasks-1?taskID+1:-1; 
            TaskName left(leftNeighbor), right(rightNeighbor);
	    if(taskID==0){
		_token=1;
	    }else{
		if(leftNeighbor!=-1){
		    Pull(left, (char*)&_token, sizeof(int));
	            cout<< endl << "Task " << taskID << " has received token " << _token <<" from Task "<<leftNeighbor <<endl;
		}
	    }
	    if(rightNeighbor!=-1){
		Push(right, (char*)&_token, sizeof(int));
	    }
	}
	bool Dependency(){
	    int taskID= MyName()[0];
	    int leftNeighbor= taskID>0?taskID-1:-1; 
            TaskName left(leftNeighbor);
	    if(leftNeighbor!=-1) return Depend_on(left);
            else return true;
	}
	void PostCompletion(){
	    //Do nothing. This task will be destroyed and no further action will be made.
	}
};

int TokenRingTask::_nTasks;

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
	if(!(strcmp(argv[argCount], "-t"))) t = atoi(argv[++argCount]);
	if(!(strcmp(argv[argCount], "-v"))) verbose = true;
    }
    RTS rts;
    rts.RTS_Init(&rank, &nProcs);
    TokenRingTask::_nTasks= t;
    string graphName= "TokenRing";
    if(verbose && rank==0){
	printf("Creating a 1D Token Ring Graph with %d tasks\n", t);
	printf("Running the graph with %d processes\n", rts.RTS_ProcCount());
    }
    ArrayGraph<TokenRingTask> *TokenRingGraph= new ArrayGraph<TokenRingTask>(graphName, t, rank, nProcs);
    rts.RTS_Run(TokenRingGraph, false);
    if(verbose && rank==0) printf("Graph execution completed\n");
    delete TokenRingGraph;
    rts.RTS_Finalize();
};
