//BL_COPYRIGHT_NOTICE

#ifndef _LIST_H_
#define _LIST_H_

//
// $Id: List.H,v 1.1 1997-09-12 18:00:10 lijewski Exp $
//

#include <Boolean.H>
#include <Assert.H>
#include <BoxLib.H>

template <class T> class ListLink;
template <class T> class ListIterator;
template <class T> class List;

template <class T>
class ListLink
{
private:
    friend class List<T>;
    friend class ListIterator<T>;

    ListLink (const T&     _val,
              ListLink<T>* _pre,
              ListLink<T>* _suc)
        : val(_val),
          pre(_pre),
          suc(_suc)
    {}

private:
    T            val;
    ListLink<T>* pre;
    ListLink<T>* suc;
};

//
//@Man:
//@Memo: Iterate over a List
/*@Doc:

  The class ListIterator<T> is an iterator over class List<T>.
  
  This class does NOT provide a default constructor or an assignment operator.
*/

template <class T>
class ListIterator
{
public:
    //
    //@ManDoc: Construct a ListIterator<T> to first element of aList.
    //
    inline ListIterator (const List<T>& aList);
    //
    //@ManDoc: The copy constructor.
    //
    inline ListIterator (const ListIterator<T>& rhs);

    /*@ManDoc: Reset this ListIterator<T> to point to the first element
               in the List<T>.
    */
    inline void rewind ();

    /*@ManDoc: Return a constant reference to the object in the List<T>
               currently pointed to by this ListIterator<T>.
    */
    inline const T& operator() () const;

    /*@ManDoc: Return a constant reference to the object in the List<T>
               currently pointed to by this ListIterator<T>.
    */
    inline const T& operator* () const;

    /*@ManDoc: This is a conversion operator that makes the iterator look
               like a pointer.  This operator makes it easy to check if the
               iterator is pointing to an element on the List<T>.  If the
               iterator has been moved off the List<T> or if the List<T> is
               empty, this conversion returns the NULL pointer.
    */
    inline operator void* ();

    /*@ManDoc: Returns true if ListIterator<T> doesn't point to any
               element on the List<T>.
    */
    inline bool operator! () const;

    /*@ManDoc: Return a constant reference to the object in the List<T>
               currently pointed to by the iterator.
    */
    inline const T& value () const;

    /*@ManDoc: This is the prefix auto-increment operator.  It advances the
               ListIterator<T> to point to the next element on the List<T>.
               It then returns a reference to itself to allow for chaining
               with other operators.
    */
    inline ListIterator<T>& operator++ ();

    /*@ManDoc: This is the prefix auto-decrement operator.  It moves the
               ListIterator<T> to point to the previous element on the
               List<T>.  It then returns a reference to itself to allow for
               chaining with other operators.
    */
    inline ListIterator<T>& operator-- ();

    /*@ManDoc: This is the postfix auto-decrement operator.  It moves the
               ListIterator<T> to point to the previous element on the
               List<T>.  It then returns a ListIterator<T> that points to
               the old element to allow for chaining with other operators.

    */
    inline ListIterator<T> operator-- (int);

    /*@ManDoc: This is the postfix auto-increment operator.  It advances the
               ListIterator<T> to point to the next element on the List<T>.
               It then returns a ListIterator<T> that points to the old
               element to allow for chaining with other operators.
    */
    inline ListIterator<T> operator++ (int);

    /*@ManDoc: Do the two ListIterator<T>s point to the same List<T> and
               the same element within the List<T>?
    */
    inline bool operator== (const ListIterator<T>&) const;
    //
    //@ManDoc: Are the ListIterator<T>s not equal?
    //
    inline bool operator!= (const ListIterator<T>&) const;

protected:
    //
    // Construct a ListIterator<T> to a List<T> and object in that List<T>.
    //
    inline ListIterator (const List<T>& _list,
                         ListLink<T>*   _p);
    //
    // A reference to the List<T> to which we point.
    //
    const List<T>& list;
    //
    // A pointer to the element in the List<T> to which we point.
    //
    ListLink<T>* p;

private:
    friend class List<T>;
    //
    // These are disallowed.
    //
    ListIterator ();
    ListIterator<T>& operator= (const ListIterator<T>&);
};

//
//@Man:
//@Memo: A Doubly-Linked List
/*@Doc:

  The List<T> class is a template class that implements a doubly-linked list
  of objects.  A List<T> is a useful container class when the number of
  objects in the collection is not known ahead of time.   A List<T> can
  contain an arbitrary number of elements; operations such as insertion,
  deletion, and catenation are easily implemented and inexpensive.

  The only difficulty when defining a list class is devising a mechanism to
  access the elements.  In an array, an element is accessed using an
  integer index.  Since the elements in a List<T> are ordered by position,
  we could define an integer indexing operation that walks along the
  List<T> links from the beginning until the numbered element is found.
  Unfortunately, this would be very inefficient when accessing elements near
  the  end of a long list.  Another solution is to allow user access to the
  individual link objects that contain the element as well as the forward and
  backward pointers.  This is not a satisfactory solution since it allows
  user access to the internal representation of the class.  The solution
  chosen is to define a ListIterator<T> template class.
  
  Think of a ListIterator<T> as a pointer to an object in the List<T>.  You
  can access the element currently pointed to by the iterator, move the
  iterator forward and backward through the List<T>, and use it as a
  mechanism to define where elements should be inserted and deleted.  If the
  iterator is moved off the end of the list it behaves as a null pointer.

  This is a concrete class, not a polymorphic one.
*/

template <class T>
class List
{
public:
    //
    //@ManDoc: Construct an empty List<T>.
    //
    inline List ();
    //
    //@ManDoc: The copy constructor.
    //
    List (const List<T>& rhs);
    //
    //@ManDoc: The assignment operator.
    //
    List<T>& operator= (const List<T>& rhs);
    //
    //@ManDoc: The destructor.
    //
    inline ~List();
    //
    //@ManDoc: Adds a copy of the value to the beginning of the List<T>.
    //
    inline void prepend (const T& value);
    //
    //@ManDoc: Adds a copy of the value to the end of the List<T>.
    //
    inline void append (const T& value);
    //
    //@ManDoc: Adds a copy of the value to the end of the List<T>.
    //
    void add (const T& value);
    //
    //@ManDoc: Appends a copy of all items in List<T> src to this List<T>.
    //
    void join (const List<T>& src);

    /*@ManDoc: Appends a copy of all items in List<T> src to this List<T>.
               This differs from join() in that it unlinks the objects from
               the List<T> src and glues them to the end of this List<T>,
               leaving List<T> src empty.  This is more efficient that
               join() if src is no longer needed.
    */
    void catenate (List<T>& src);
    //
    //@ManDoc: Removes all objects from the List<T>.
    //
    void clear ();

    /*@ManDoc: Returns a copy of this List<T> on the heap.  It is the user's
               responsibility to delete this when no longer needed.
    */
    inline List<T>* copy () const;
    //
    //@ManDoc: Returns a reference to the first element in the List<T>.
    //
    inline T& firstElement () const;
    //
    //@ManDoc: Returns a reference to the last element in the List<T>.
    //
    inline T& lastElement () const;

    /*@ManDoc: Returns true if the List<T> contains an object identical
               to value.  Type T must have an operator==() defined, or
               be an intrinsic type.
    */
    bool includes (const T& value) const;

    /*@ManDoc: Returns true if the this and rhs are memberwise equal;
               i.e. the two lists are the same size and each of the
               elements in the list compare equal. Type T must have an
               operator==() defined, or be an intrinsic type.
    */
    bool operator== (const List<T>& rhs) const;
    //
    //@ManDoc: Returns true if the this and rhs are not equal.
    //
    bool operator!= (const List<T>& rhs) const;
    //
    //@ManDoc: Returns true if the List<T> is empty.
    //
    inline bool isEmpty () const;
    //
    //@ManDoc: Returns true if the List<T> is not empty.
    //
    inline bool isNotEmpty () const;
    //
    //@ManDoc: Returns the number of objects in the List<T>.
    //
    int length () const;
    //
    //@ManDoc: Removes the first element in the List<T>.
    //
    inline void removeFirst ();
    //
    //@ManDoc: Removes the last element in the List<T>.
    //
    inline void removeLast ();
    //
    //@ManDoc: Returns reference to object pointed to by the ListIterator<T>.
    //
    inline const T& operator[] (const ListIterator<T>& li) const;
    //
    //@ManDoc: Returns reference to object pointed to by the ListIterator<T>.
    //
    inline T& operator[] (const ListIterator<T>& li);
    //
    //@ManDoc: Removes all objects in the List<T> equal to value.
    //
    void remove (const T& value);

    /*@ManDoc: Removes all objects in the List<T> equal to any of the
               values in lst.
    */
    void remove (const List<T>& lst);
    //
    //@ManDoc: Removes the object pointed to by the ListIterator<T>.
    //
    void remove (ListIterator<T>& lit);
    //
    //@ManDoc: Replace the value pointed to by the ListIterator<T> by val.
    //
    inline void replace (ListIterator<T>& li,
                         const T&         val);

    /*@ManDoc: Insert val into List<T> after the object pointed to by
               ListIterator<T>.
    */
    inline void addAfter (ListIterator<T>& lit,
                          const T&         val);

    /*@ManDoc: Insert val into List<T> before the object pointed to by
               ListIterator<T>.
    */
    inline void addBefore (ListIterator<T>& lit,
                           const T&         val);
    //
    //@ManDoc: Returns a ListIterator<T> to the first object in this List<T>.
    //
    inline ListIterator<T> first () const;
    //
    //@ManDoc: Returns a ListIterator<T> to the last object in this List<T>.
    //
    inline ListIterator<T> last () const;

protected:
    //
    // A helper function for removing nodes.
    //
    void remove (ListLink<T> *ln);
    //
    // A helper function for adding nodes.
    //
    ListLink<T>* addBefore (ListLink<T>* ln,
                            const T&     val);
    //
    // A helper function for adding nodes.
    //
    ListLink<T>* addAfter (ListLink<T>* ln,
                           const T&     val);
    //
    // The head of the list.
    //
    ListLink<T>* head;
    //
    // The tail of the list.
    //
    ListLink<T>* tail;
    //
    // Our good and trusted friend.
    //
    friend class ListIterator<T>;
};

//
// Inlines.
//

//
// The ListIterator stuff.
//

template <class T>
inline
ListIterator<T>::ListIterator (const List<T>& _list,
                               ListLink<T>*   _p)
    : list(_list),
      p(_p)
{}

template <class T>
inline
ListIterator<T>::ListIterator (const List<T>& aList)
    : list(aList)
{
    p = list.head;
}

template <class T>
inline
ListIterator<T>::ListIterator (const ListIterator<T>& li)
    : list(li.list),
      p(li.p)
{}

template <class T>
inline
void
ListIterator<T>::rewind ()
{
    p = list.head;
}

template <class T>
inline
const T&
ListIterator<T>::operator() () const
{
    assert(p != 0);
    return p->val;
}

template <class T>
inline
const T&
ListIterator<T>::operator* () const
{
    assert(p != 0);
    return p->val;
}

template <class T>
inline
ListIterator<T>::operator void* ()
{
    return p != 0 ? this : 0;
}

template <class T>
inline
bool
ListIterator<T>::operator! () const
{
    return p == 0 ? true : false;
}

template <class T>
inline
const T&
ListIterator<T>::value () const
{
    assert(p != 0);
    return p->val;
}

template <class T>
inline
ListIterator<T>&
ListIterator<T>::operator++ ()
{
    if (p)
        p = p->suc;
    return *this;
}

template <class T>
inline
ListIterator<T>&
ListIterator<T>::operator-- ()
{
    if (p)
        p = p->pre;
    return *this;
}

template <class T>
inline
ListIterator<T>
ListIterator<T>::operator++ (int)
{
    const ListIterator<T> li = *this;
    ++(*this);
    return li;
}

template <class T>
inline
ListIterator<T>
ListIterator<T>::operator-- (int)
{
    const ListIterator<T> li = *this;
    --(*this);
    return li;
}

template <class T>
inline
bool
ListIterator<T>::operator== (const ListIterator<T>& _li) const
{
    return (&list == &_li.list && p == _li.p) ? true : false;
}

template <class T>
inline
bool
ListIterator<T>::operator!= (const ListIterator<T>& _li) const
{
    return ! ListIterator<T>::operator==(_li);
}

//
// List stuff.
//

template <class T>
inline
List<T>::List ()
    : head(0),
      tail(0)
{}

template <class T>
inline
List<T>::~List ()
{
    clear();
}

template <class T>
inline
void
List<T>::prepend (const T& value)
{
    addBefore(head, value);
}

template <class T>
inline
void
List<T>::append (const T& value)
{
    addAfter(tail, value);
}

template <class T>
inline
List<T>*
List<T>::copy () const
{
    List<T>* newlist = new List<T>(*this);
    if (newlist == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    return newlist;
}

template <class T>
inline
T&
List<T>::firstElement () const
{
    assert(head != 0);
    return head->val;
}

template <class T>
inline
T&
List<T>::lastElement () const
{
    assert(tail != 0);
    return tail->val;
}

template <class T>
inline
bool
List<T>::isEmpty () const
{
    return head == 0 && tail == 0;
}

template <class T>
inline
bool
List<T>::isNotEmpty () const
{
    return !isEmpty();
}

template <class T>
inline
void
List<T>::removeFirst ()
{
    remove(head);
}

template <class T>
inline
void
List<T>::removeLast ()
{
    remove(tail);
}

template <class T>
inline
const T&
List<T>::operator[] (const ListIterator<T>& li) const
{
    assert(li.p != 0);
    return li.p->val;
}

template <class T>
inline
T&
List<T>::operator[] (const ListIterator<T>& li)
{
    assert(li.p != 0);
    return li.p->val;
}

template <class T>
inline
void
List<T>::replace (ListIterator<T>& li,
                  const T&         _val)
{
    assert(li.p != 0);
    li.p->val = _val;
}

template <class T>
inline
void
List<T>::addAfter (ListIterator<T>& lit,
                   const T&         val)
{
    addAfter(lit.p, val);
}

template <class T>
inline
void
List<T>::addBefore (ListIterator<T>& lit,
                    const T&         val)
{
    addBefore(lit.p, val);
}

template <class T>
inline
ListIterator<T>
List<T>::first () const
{
    return ListIterator<T>(*this,head);
}

template <class T>
inline
ListIterator<T>
List<T>::last () const
{
    return ListIterator<T>(*this,tail);
}


#include <Assert.H>

//
// List members
//

template <class T>
List<T>::List (const List<T>& source)
    : head(0),
      tail(0)
{
    if (source.isEmpty())
        tail = head = 0;
    else
        for (ListIterator<T> li(source); li; ++li)
            append(li());
}

//
// This isn't inlined as it's declared virtual.
//

template <class T>
void
List<T>::add (const T& value)
{
    append(value);
}

template <class T>
int
List<T>::length () const
{
    int len = 0;
    for (ListIterator<T> li(*this); li; ++li)
        len++;
    return len;
}

template <class T>
List<T>&
List<T>::operator= (const List<T>& source)
{
    if (!(this == &source))
    {
        clear();
        for (ListIterator<T> li(source); li; ++li)
            append(li());
    }
    return *this;
}

template <class T>
ListLink<T> *
List<T>::addBefore (ListLink<T>* ln,
                    const T&     val)
{
    assert(ln != 0 || head == 0);

    ListLink<T>* newlink;

    if (ln == head)
    {
        head = newlink = new ListLink<T>(val, 0, head);

        if (tail == 0)
            tail = head;
        else
            head->suc->pre = newlink;
    }
    else
    {
        newlink = new ListLink<T>(val, ln->pre, ln);

        ln->pre->suc = newlink;
        ln->pre = newlink;
    }

    if (newlink == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    return newlink;
}

template <class T>
ListLink<T>*
List<T>::addAfter (ListLink<T>* ln,
                   const T&     val)
{
    assert(ln != 0 || tail == 0);

    ListLink<T>* newlink;

    if (ln == tail)
    {
        tail = newlink = new ListLink<T>(val,tail,0);

        if (head == 0)
            head = tail;
        else
            tail->pre->suc = newlink;
    }
    else
    {
        newlink = new ListLink<T>(val, ln, ln->suc);

        ln->suc->pre = newlink;
        ln->suc = newlink;
    }

    if (newlink == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    return newlink;
}

template <class T>
void
List<T>::join (const List<T>& list2)
{
    for (ListIterator<T> li2(list2); li2; ++li2)
        append(li2());
}

template <class T>
void
List<T>::catenate (List<T>& list2)
{
    if (list2.isEmpty())
        //
        // Do nothing.
        //
        ;
    else if (isEmpty())
    {
       head = list2.head;
       tail = list2.tail;
       list2.head = 0;
       list2.tail = 0;
    }
    else
    {
        tail->suc = list2.head;
        list2.head->pre = tail;
        tail = list2.tail;
        list2.head = 0;
        list2.tail = 0;
    }
}

template <class T>
void
List<T>::clear ()
{
    ListLink<T>* next = 0;

    for (ListLink<T>* p = head; p != 0; p = next)
    {
        next = p->suc;
        p->suc = 0;
        delete p;
    }
    tail = head = 0;
}

template <class T>
bool
List<T>::includes (const T& v) const
{
    bool rc = false;
    for (ListIterator<T> li(*this); li && !rc; ++li)
        if (v == li())
            rc = true;
    return rc;
}

template<class T>
bool
List<T>::operator== (const List<T>& rhs) const
{
    if (length() == rhs.length())
    {
        for (ListIterator<T> li(*this), ri(rhs); li; ++li, ++ri)
            if (li() != ri())
                return false;
        return true;
    }

    return false;
}

template<class T>
bool
List<T>::operator!= (const List<T>& rhs) const
{
    return !operator==(rhs);
}

template <class T>
void
List<T>::remove (ListIterator<T>& li)
{
    ListLink<T> *np = li.p->suc;
    remove(li.p);
    li.p = np;
}

template <class T>
void
List<T>::remove (const T& _v)
{
    for (ListIterator<T> litr(*this); litr; ++litr)
        if (litr() == _v)
            remove(litr);
}

template <class T>
void
List<T>::remove (const List<T>& _lv)
{
    for (ListIterator<T> litr(_lv); litr; ++litr)
        remove(litr());
}

template <class T>
void
List<T>::remove (ListLink<T>* ln)
{
    assert(head !=0 && tail != 0);

    if (head == ln && tail == ln)
        head = tail = 0;
    else if (head == ln)
    {
        assert(ln->pre == 0);
        head = ln->suc;
        head->pre = 0;
    }
    else if (tail == ln)
    {
        assert(ln->suc == 0);
        tail = ln->pre;
        tail->suc  = 0;
    }
    else
    {
        assert(ln->suc != 0 && ln->pre != 0);
        ln->suc->pre = ln->pre;
        ln->pre->suc = ln->suc;
    }
    delete ln;
    ln = 0;
}

#endif /*_LIST_H_*/
