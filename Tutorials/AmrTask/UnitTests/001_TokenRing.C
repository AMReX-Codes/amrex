//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 08-14-2017
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
	    SelfDestroy();
	}
};

int TokenRingTask::_nTasks;


/* Example commands
   -Serial run:
   ./001_TokenRing -t 8 -v
 */

int main(int argc,char *argv[])
{
    int argCount = 0;
    int t=1, verbose=0;
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
    rts.Init();
    rank= rts.MyProc();
    nProcs= rts.ProcCount();
    TokenRingTask::_nTasks= t;
    string graphName= "TokenRing";
    if(verbose && rank==0){
	cout<<"Creating a 1D Token Ring Graph with "<<t <<" tasks";
	cout<<"Running the graph with "<< nProcs<<" processes";
    }
    double time = -rts.Time();
    rts.Barrier();
    ArrayGraph<TokenRingTask> *TokenRingGraph= new ArrayGraph<TokenRingTask>(graphName, t, rank, nProcs);
    rts.Iterate(TokenRingGraph);
    rts.Barrier();
    time += rts.Time();
    if(verbose && rank==0) cout<<"Graph execution takes "<< time <<" seconds"<<endl;
    delete TokenRingGraph;
    rts.Finalize();
};
