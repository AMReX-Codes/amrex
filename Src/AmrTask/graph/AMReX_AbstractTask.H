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
	    TaskName(int id0, int id1, int id2, int id3){_ids.push_back(id0); _ids.push_back(id1); _ids.push_back(id2);_ids.push_back(id3);}
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
	    int _serializedDescSize;
	    void* SetBuffer(char* buffer){_buffer= buffer;}  
	public:
	    //!Create a message with empty load
	    Data(TaskName src, TaskName recipient, size_t size) {
		_src= src; _recipient= recipient;
		_size=size; _tag=0;
		_serializedDescSize= sizeof(size_t) /*_size*/+ sizeof(int) /*_tag*/ + 2*sizeof(int) /*_src and _recipient length*/ + (src.Dim()+recipient.Dim()) * sizeof(size_t);
		_buffer= new char[size+_serializedDescSize];
	    }

	    Data (char* a){ //deserialize data
		_buffer= a;
		_size= *((size_t*) a);
		_tag= *((int*) (a+sizeof(size_t)));
		int srclen= *((int*) (a+sizeof(size_t)+ sizeof(int)));
		int recipientlen= *((int*) (a+sizeof(size_t)+ 2*sizeof(int)));
		_src.SetSize(srclen); _recipient.SetSize(recipientlen);
		size_t *entries= (size_t*) (a+sizeof(size_t)+ 3*sizeof(int));
		for(int i=0; i<srclen+recipientlen; i++){
		    if(i<srclen) _src[i]= entries[i];
		    else _recipient[i-srclen]= entries[i];
		}
		_serializedDescSize= sizeof(size_t) /*_size*/+ sizeof(int) /*_tag*/ + 2*sizeof(int) /*_src and _recipient length*/ + (_src.Dim()+_recipient.Dim()) * sizeof(size_t);
	    }

	    TaskName GetRecipient(){return _recipient;}
	    TaskName GetSource(){return _src;}
	    void* GetBuffer(){return (void*)(_buffer+_serializedDescSize);}  
	    size_t GetSize(){return _size;}
	    void SetTag(size_t tag){_tag= tag;}
	    int GetTag(){return _tag;}
	    void Free(){
		delete _buffer;
		delete this;
	    }

	    char* SerializeData(){
		char* p=_buffer;
		memcpy(p, &_size, sizeof(size_t));
		p+= sizeof(size_t);
		memcpy(p, &_tag, sizeof(int));
		p+= sizeof(int);
		int srclen= _src.Dim(), recipientlen= _recipient.Dim();
		memcpy(p, &srclen, sizeof(int));
		p+= sizeof(int);
		memcpy(p, &recipientlen, sizeof(int));
		p+= sizeof(int);
		size_t *entries= (size_t*) p;
		for(int i=0; i<srclen+recipientlen; i++){ 
		    if(i<srclen) entries[i]= _src[i];
		    else entries[i]= _recipient[i-srclen];
		}
		return _buffer;
	    }
	    size_t GetSerializedSize(){return _size+ _serializedDescSize;}
    };

    class DependencyMap{
	private:
	    std::map<TaskName, map<int, std::queue< Data* > > > _dataMap;
	    size_t _size; 

	public:
	    //return the total number of inputs/outputs of a task
	    size_t size(){
		size_t s=0;
		std::map<TaskName, std::map<int, std::queue< Data* > > >::iterator it= _dataMap.begin();
		while (it!= _dataMap.end()){
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
		if(_dataMap.find(name) != _dataMap.end()){
		    if(_dataMap[name].find(tag) != _dataMap[name].end()){
			return _dataMap[name][tag].size() >0;
		    }
		}
		return false;
	    }
	    bool empty(){return _size==0;}
	    void push_back(TaskName name, Data* d, int tag=0){
		_dataMap[name][tag].push(d);
		_size++;
	    }
	    Data* pop_front(TaskName name, int tag=0){
		if(_dataMap[name][tag].size()==0) return NULL;
		Data* d= _dataMap[name][tag].front();
		_dataMap[name][tag].pop();
		_size--;
		return d;
	    }
	    std::queue< Data* >& GetDependencies(TaskName name, int tag=0){return _dataMap[name][tag];}
    };

    //! The most abstract task
    class Task{
	protected:
	    TaskName _id;
	    DependencyMap _neighbors_in; 
	    std::queue<Data*> _outputs;
	    std::queue<Task*> _newTasks;
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
		return _neighbors_in.isSatisfied(src, tag);
	    }
	    std::queue<Data*>& GetOutputs(){return _outputs;}
	    std::queue<Task*>& GetNewTasks(){return _newTasks;}
	    DependencyMap& GetInputs(){return _neighbors_in;}
	    void KeepTaskAlive(){_isPersistent=true;}
	    void SelfDestroy(){_isPersistent=false;}
	    bool isPersistent(){return _isPersistent;}

	    template<typename T> void LocalAtomicAdd(T *addr, T val){
		LocalAtomicAdd_impl(addr, val);
	    }
	    template<typename T> void GlobalAtomicAdd(T *addr, T val){
		GlobalAtomicAdd_impl(addr, val);
	    }
	    void RegisterTask(Task* t){
		_newTasks.push(t);
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
