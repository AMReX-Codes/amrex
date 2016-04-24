
#include <BoxArray.H>

extern "C" {

    void fi_new_boxarray (BoxArray*& ba, int lo[3], int hi[3])
    {
	IntVect small(lo), big(hi);
	ba = new BoxArray(Box(small,big));
    }

    void fi_delete_boxarray (BoxArray* ba)
    {
	delete ba;
    }

    void fi_boxarray_maxsize (BoxArray* ba, int sz)
    {
	ba->maxSize(sz);
    }
}
