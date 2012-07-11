
#ifndef BL_PLIST_H
#define BL_PLIST_H

//
// An enum that controls whether or not class PList<T>
// manages the pointed to objects, as well as the pointer itself.
// Valid values are PListNoManage and PListManage.
//
enum PListPolicy
{
    PListNoManage = 0,
    PListManage   = 1
};

/*
  An Array of Pointers to Objects of Type T

  The PList<T> class is a wrapper for a std::list of pointers. It provides a memory 
  management so that clear and pop operations delete the objects pointed at.
*/

template <class T>
class PList
    :public std::list<T*>
{
public:
    //
    // Constructs a PList<T> with no elements.  The default memory
    // management policy is not to delete the underlying objects
    // upon destruction.
    //
    explicit PList (PListPolicy _managed = PListNoManage);
    //
    // Constructs a PList<T> of length len of null pointers.  
    // The default memory management policy is not to delete the
    // underlying objects upon destruction.
    //
    explicit PList (int          len,
                     PListPolicy _managed = PListNoManage);
    //
    // Overloads the default pop to delete the item if managed.
    //
    void pop_back();
    void pop_front();
    
    //
    // The destructor.  If the memory management policy was
    // set to Manage upon construction, this destructor will
    // also delete all objects referenced by the pointers in
    // this list.
    //
    ~PList ();
    //
    // This function sets each pointer in this PList<T> to null.
    // If the memory management policy was set to Manage
    // upon construction, this function will also delete all objects
    // referenced by the pointers in this list.
    //
    void clear ();

private:
    //
    // The memory management policy.
    //
    PListPolicy managed;

public:
    //
    // Disallow these.
    //
    PList (const PList<T>&);
    // PList<T> operator= (const PList<T>&);
};

template <class T>
PList<T>::PList (const PList<T>& c)
{
    BL_ASSERT(c.size()==0);

    managed = c.managed;
}

template <class T>
PList<T>::PList (PListPolicy _managed)
    :
    managed(_managed)
{}

template <class T>
PList<T>::PList (int          len,
                   PListPolicy _managed)
    :
    std::list<T*>(len,0),
    managed(_managed)
{
}

template <class T>
void
PList<T>::pop_back()
{
    if (managed)
    {
        delete ((T*)(this->back()));
        this->std::list<T*>::pop_back();
    }
    else
    {
        this->std::list<T*>::pop_back();
    }    
}

template <class T>
void
PList<T>::pop_front()
{
    if (managed)
    {
        delete ((T*)(this->front()));
        this->pop_front();
    }
    else
    {
        this->pop_front();
    }    
}

template <class T>
void
PList<T>::clear ()
{
    while (!this->empty())
    {
        this->pop_back();
    }
}

template <class T>
PList<T>::~PList ()
{
    clear();
}


#endif /*BL_PARRAY_H*/
