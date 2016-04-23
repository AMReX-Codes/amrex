
#include <BoxArray.H>

extern "C" {

    void fi_new_boxarray (void*& ba, int lo[3], int hi[3])
    {
	IntVect small(lo), big(hi);
	ba = (void*) new BoxArray(Box(small,big));
    }

    void fi_delete_boxarray (void* ba)
    {
	delete (BoxArray*) ba;
    }

}
