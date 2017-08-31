#ifndef MYATOMICS
#define MYATOMICS

//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-21-2017

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>
#include "mylock.h"

namespace amrex{
    extern MyLock _l;
    template<typename T> void LocalAtomicAdd_impl(T *addr, T val){
        _l.lock();
        *addr+= val;
        _l.unlock();
    }
    template<typename T> void GlobalAtomicAdd_impl(T *addr, T val){
	assert(false);//not defined
    }

}//end namespace

#endif
