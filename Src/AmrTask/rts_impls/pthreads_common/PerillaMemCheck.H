#ifndef P_MEMCHECK_H
#define P_MEMCHECK_H
#include <string>
#include <sstream>
#include <map>
#include "mylock.h"
using namespace std;

namespace perilla{

    class PerillaMemCheck
    {
	public:
	void add(string key, void* obj, string classname);
	void remove(string key);
	string genKey(void *obj){
	    std::stringstream ost;
	    ost<<obj;
	    return ost.str();
	}
	PerillaMemCheck():addCnt(0), rmCnt(0){}
	void report();
	private:
	MyLock lock;
	std::map<string, void* > objMap; 
	int addCnt, rmCnt;
    };
}//end namespace

#endif
