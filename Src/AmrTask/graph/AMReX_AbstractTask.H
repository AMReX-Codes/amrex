#ifndef AMREX_ABSTRACT_TASK
#define AMREX_ABSTRACT_TASK
//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 08-07-2017

#include "AMReX_DataTypes.H"
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <cassert>
#include "rts_taskimpl.H"
using namespace std;

namespace amrex{

    class TaskName{
	protected:
	    std::vector<int> _ids;
	public:
	    TaskName(){}
	    TaskName(int id){_ids.push_back(id);}
	    TaskName(int id0, int id1){_ids.push_back(id0); _ids.push_back(id1);}
	    TaskName(int id0, int id1, int id2){_ids.push_back(id0); _ids.push_back(id1); _ids.push_back(id2);}
	    int& operator[](int i){return _ids[i];}
	    const int& operator[](int i) const{return _ids[i];}
	    TaskName& operator=(const TaskName &rhs){
		assert(_ids.size()==0 || _ids.size()== rhs.Dim());
		if(_ids.size()==0)
		    for(int i=0; i< rhs.Dim(); i++)
			_ids.push_back(rhs[i]);
		else //already have the same dim
		    for(int i=0; i< _ids.size(); i++)
			_ids[i] = rhs[i];
		return *this;
	    }
	    bool operator== (const TaskName &rhs) const{
		if(_ids.size() != rhs.Dim()) return false;
		for(int i=0; i< _ids.size(); i++)
		    if(_ids[i] != rhs[i]) return false;
		return true;
	    }
	    bool operator!= (const TaskName &rhs) const{
		if(_ids.size() != rhs.Dim()) return true;
		for(int i=0; i< _ids.size(); i++)
		    if(_ids[i] != rhs[i]) return true;
		return false;
	    }
	    bool operator< (const TaskName &rhs) const{
		if(_ids.size() < rhs.Dim()) return true;
		if(_ids.size() > rhs.Dim()) return false;
		for(int i= _ids.size()-1; i>=0; i--){
		    if(_ids[i] > rhs[i]) return false;
		    if(_ids[i] < rhs[i]) return true;
                }
		return false;
	    }
            void SetSize(int dim){
		for(int i=0; i<dim; i++) _ids.push_back(0);
            }
	    int Dim() const{return _ids.size();}
    };

    //! This class can be used to express data dependencies. The runtime may serialize the data or augment it with other control information
    class Data{
	private:
	    size_t _size; //in bytes
	    TaskName _src, _recipient;
	    int _tag;
	    char* _buffer;
	public:
	    //!Create a message with empty load
	    Data(){_buffer=NULL; _size=0; _tag=0;}
	    //!Create a message with a pre-exsiting buffer
	    Data(void* buffer, size_t size){_buffer= (char*)buffer; _size=size; _tag=0;}
	    //!Create a new message with size byte load
	    Data(size_t size) {_buffer= new char[size]; _size=size; _tag=0;}

	    void SetSource(TaskName src){_src= src;}
	    void SetRecipient(TaskName recipient){_recipient= recipient;}
	    TaskName GetRecipient(){return _recipient;}
	    TaskName GetSource(){return _src;}
	    void* GetBuffer(){return (void*)_buffer;}  
	    void SetSize(size_t size){_size= size;}
	    void SetTag(size_t tag){_tag= tag;}
	    int GetTag(){return _tag;}
	    void Free(){delete _buffer;}
    };

    class DependencyMap{
	struct TaskComp {
	    bool operator() (const TaskName& lhs, const TaskName& rhs) const
	    {
		return lhs<rhs;
	    }
	};
	private:
	std::map<TaskName, map<int, std::queue< Data* > > > dataMap;
	size_t _size; 

	public:
	//return the total number of inputs/outputs of a task
	size_t size(){
	    size_t s=0;
	    std::map<TaskName, std::map<int, std::queue< Data* > > >::iterator it= dataMap.begin();
	    while (it!= dataMap.end()){
		std::map<int, std::queue< Data* > >::iterator tagIt= ((*it).second).begin(); 
		while (tagIt!= ((*it).second).end()){
		    s+= (*tagIt).second.size();
		    tagIt++;
		}
		it++;
	    }
	    return s;
	}
	bool isSatisfied(TaskName name, int tag=0){
	    if(dataMap.find(name) != dataMap.end()){
		if(dataMap[name].find(tag) != dataMap[name].end()){
		    return dataMap[name][tag].size() >0;
                }
	    }
	    return false;
	}
	bool empty(){return _size==0;}
	void push_back(TaskName name, Data* d, int tag=0){
	    dataMap[name][tag].push(d);
	    _size++;
	}
	Data* pop_front(TaskName name, int tag=0){
	    if(dataMap[name][tag].size()==0) return NULL;
	    Data* d= dataMap[name][tag].front();
	    dataMap[name][tag].pop();
	    _size--;
	    return d;
	}
	std::queue< Data* >& GetDependencies(TaskName name, int tag){return dataMap[name][tag];}
    };

    //! The most abstract task
    class Task{
	protected:
	    TaskName _id;
	    DependencyMap neighbors_in; 
	    std::queue<Data*> outputs;
            bool _isPersistent;
	public:
            Task():_isPersistent(false){}
	    //Describe Data Dependency
	    virtual bool Dependency()=0;
	    //! What the task is supposed to do
	    virtual void Job()=0;
	    //! Once the task finished its computation, any actions should be taken (like create new taks)?
	    virtual void PostCompletion()=0;
	    TaskName MyName(){return _id;}
	    void SetName(TaskName id){ _id=id;}

	    bool TestDependencies(){return Dependency();}
	    void RunJob(){Job();}
	    void RunPostCompletion(){PostCompletion();} 
	    void Pull(TaskName src, char* d, size_t size, int tag=0);
	    void Push(TaskName dest, char* d, size_t size, int tag=0);
	    bool Depend_on(TaskName src, int tag=0){
		return neighbors_in.isSatisfied(src, tag);
	    }
	    std::queue<Data*>& GetOutputs(){return outputs;}
	    DependencyMap& GetInputs(){return neighbors_in;}
            void KeepTaskAlive(){_isPersistent=true;}
            void DestroyTask(){_isPersistent=false;}
            bool isPersistent(){return _isPersistent;}

            template<typename T> void LocalAtomicAdd(T *addr, int val){
		LocalAtomicAdd_impl(addr, val);
            }
            template<typename T> void GlobalAtomicAdd(T *addr, int val){
		GlobalAtomicAdd_impl(addr, val);
	    }
    };

#if 0
    //! This task is created to do LOCAL jobs (thus it is not migratable), when all data dependencies have been satisfied
    class NonMigratableTask :public Task{
	protected:
	    TaskState st;
	    //! A task may depend on multiple pieces of data from another task
	    //! A task may send multiple pieces of data to another task
	    //DependencyMap neighbors_out;
	public:
    };


    //! This task can be migrated to a remote process (e.g. to be closer to data). 
    class MigratableTask :public Task{
	protected:
	    int origin;//original process
	public:
	    virtual void pullData();  
    };




    enum state_t{FRESH /*just created*/, PENDING/*waiting for data*/, READY/*data dependencies satisfied*/, RUNNING /*scheduled to run*/, FINISHED /*no more work to do*/, ERROR};
    class TaskState{
	private:
	    state_t _st;

	public:
	    TaskState():_st(PENDING){}
	    TaskState(state_t st):_st(st){}
	    void shift();
	    void shift(state_t newState);
    };
#endif




}//end namespace

#endif
