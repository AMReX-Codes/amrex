#include<AMReX_AbstractTask.H>
//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-24-2017

namespace amrex{
    void Task::Pull(TaskName src, char* d, size_t size, int tag){
	Data* data= _neighbors_in.pop_front(src, tag);
	memcpy(d, data->GetBuffer(), size);
	data->Free();
    }
    void Task::Push(TaskName dest, char* d, size_t size, int tag){
        Data* data= new Data(_id, dest, size);
	data->SetTag(tag);
	memcpy(data->GetBuffer(), d, size);
	_outputs.push(data);
    }
}
