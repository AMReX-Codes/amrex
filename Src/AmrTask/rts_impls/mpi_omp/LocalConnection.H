#ifndef P_LOCALCONNECTION_H
#define P_LOCALCONNECTION_H

#include <PackageQueue.H>
#include <Barrier.H>
#include <AMReX_Box.H>
#include <omp.h>

using namespace perilla;
#ifdef PERILLA_DEBUG
#include "PerillaMemCheck.H"
extern PerillaMemCheck memcheck;
#endif
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
	omp_lock_t ghostLock;
	LocalCopyDescriptor() : ns(-1), nd(-1), scpyCnt(0), dcpyCnt(0), sz(0), sPartner(-1), dPartner(-1) {
	    omp_init_lock(&ghostLock);
	}
};

class LocalConnection
{
    public:
	int nscpy;   //Number of cpy chunks
	int ndcpy;   //Number of cpy chunks
	omp_lock_t sLock, dLock, ghostLock;
	int firingRuleCnt;
	int scpyCnt, dcpyCnt;
	Barrier *localBarrier;
	LocalCopyDescriptor *scpy;
	LocalCopyDescriptor *dcpy;
	LocalConnection() : nscpy(0), ndcpy(0), firingRuleCnt(0), scpy(NULL), dcpy(NULL), scpyCnt(0), dcpyCnt(0), localBarrier(NULL){
	    omp_init_lock(&sLock);
	    omp_init_lock(&dLock);
	    omp_init_lock(&ghostLock);
	}

	~LocalConnection()
	{
	    if(localBarrier) free(localBarrier);
	    if(scpy) free(scpy);
	    if(dcpy) free(dcpy);
	}
};
#endif
