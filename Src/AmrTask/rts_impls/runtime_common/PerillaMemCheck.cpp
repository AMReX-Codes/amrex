#include "PerillaMemCheck.H"
using namespace perilla;

void PerillaMemCheck::add(string key, void* obj, string classname)
{
    lock.lock();
    if(objMap.find(key) == objMap.end())
    {
        objMap[key]= obj;
	addCnt++;
    }
    else{
        printf("MemCheck Error: Reinsert an object\n");
        exit(0);
    }
    lock.unlock();
}


void PerillaMemCheck::remove(string key){
    lock.lock();
    if(objMap.find(key) != objMap.end())
    {
        objMap.erase(key);
	rmCnt++;
    }
    else{
        printf("MemCheck Error: Object not found (%d Allocated vs %d Deleted)\n", addCnt, rmCnt);
        exit(0);
    }

    lock.unlock();
}
void PerillaMemCheck::report(){
    if(objMap.size()) {
        printf("Memory leak found: %d objects (%d Allocated vs %d Deleted)\n", objMap.size(), addCnt, rmCnt);
    }else printf("All allocated objects have been deallocated (%d Allocated vs %d Deleted)\n", addCnt, rmCnt);
}


