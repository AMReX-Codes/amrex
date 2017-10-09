#ifndef AMREX_CONNECTIONS_H
#define AMREX_CONNECTION_H

#include <AMReX_FArrayBox.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_IntVect.H>

namespace amrex {

    struct LocalCopyDescriptor
    {
	int ns;   // ! Source box in layout
	int nd;   //! Destination box in layout
	Box sbx;     // ! Sub-box for this copy
	Box dbx;     // ! Sub-box for this copy
        size_t sz;
	int sPartner, dPartner;
	int dcpyCnt, scpyCnt;
    };

    struct LocalConnection
    {
	int nscpy;   //! Number of cpy chunks
	int ndcpy;   //! Number of cpy chunks
        int firingRuleCnt;
	LocalCopyDescriptor *scpy;
	LocalCopyDescriptor *dcpy;
    };

    struct RemoteCommDescriptor
    {
	int ns, lns;   // ! Source box in layout
	int nd, lnd;   //! Destination box in layout
	size_t sz;
        int pr;
	Box sbx;     // ! Sub-box for this copy
	Box dbx;     // ! Sub-box for this copy
	int cnt;
    };

    struct TransDescriptor
    {
	int sz, pv, pr;
    };

    struct RemoteConnection
    {
	int nsnd;
	int nrcv;
	int nrp, nsp;
        int firingRuleCnt;
	RemoteCommDescriptor *snd;
	RemoteCommDescriptor *rcv;
	TransDescriptor *str;
	TransDescriptor *rtr;
    };
}//end namespace
#endif

