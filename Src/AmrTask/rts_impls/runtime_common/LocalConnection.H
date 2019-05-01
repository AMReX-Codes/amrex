#ifndef P_LOCALCONNECTION_H
#define P_LOCALCONNECTION_H

#include <PackageQueue.H>
#include <Barrier.H>
#include <AMReX_Box.H>
#include <pthread.h>

using namespace perilla;
#ifdef PERILLA_DEBUG
#include "PerillaMemCheck.H"
extern PerillaMemCheck memcheck;
#endif

using namespace perilla;
using namespace amrex;

class LocalCopyDescriptor
{
    public:
	int ns;   //Source box in layout
	int nd;   //Destination box in layout
	int sz;
	Box sbx;     //Sub-box for this copy
	Box dbx;     //Sub-box for this copy
	PackageQueue pQueue;  //store incoming or outgoing messages, both fab and the runtime can access this queue
	PackageQueue recycleQueue; //just for now, I'll replace this with a NUMA aware package allocator
	int sPartner, dPartner;
	int dcpyCnt,scpyCnt;
	pthread_mutex_t ghostLock;
	LocalCopyDescriptor() : ns(-1), nd(-1), scpyCnt(0), dcpyCnt(0), sz(0), sPartner(-1), dPartner(-1), ghostLock(PTHREAD_MUTEX_INITIALIZER)
	{
#ifdef PERILLA_DEBUG
//                memcheck.add(memcheck.genKey(this), (void*)this, "LocalCopyDescriptor");
#endif
	}
	~LocalCopyDescriptor(){
#ifdef PERILLA_DEBUG
//                memcheck.remove(memcheck.genKey(this));
#endif
	}
};

class LocalConnection
{
    public:
	int nscpy;   //Number of cpy chunks
	int ndcpy;   //Number of cpy chunks
	pthread_mutex_t sLock, dLock, ghostLock;
	int firingRuleCnt;
	int scpyCnt, dcpyCnt;
	Barrier *localBarrier;
	LocalCopyDescriptor *scpy;
	LocalCopyDescriptor *dcpy;
	LocalConnection() : nscpy(0), ndcpy(0), firingRuleCnt(0), scpy(NULL), dcpy(NULL), scpyCnt(0), dcpyCnt(0), localBarrier(NULL),sLock(PTHREAD_MUTEX_INITIALIZER),dLock(PTHREAD_MUTEX_INITIALIZER),ghostLock(PTHREAD_MUTEX_INITIALIZER){
//                memcheck.add(memcheck.genKey(this), (void*)this, "LocalCopyDescriptor");
}

	~LocalConnection()
	{
	    if(localBarrier) free(localBarrier);
	    if(scpy) delete [] scpy;
	    if(dcpy) delete [] dcpy;
//                memcheck.remove(memcheck.genKey(this));
	}
};
#endif
