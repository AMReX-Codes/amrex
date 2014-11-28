
#ifndef BLMAP_H
#define BLMAP_H

#include <vector>

template <class T> class BLMap;

/*
  This class implements a special purpose and lightweight map.  It is
  somewhat similar to std::map with int as keys.  But there are some
  important differences.  The most important difference is that
  elements in this class must be inserted following the order of keys.
  Another important difference is that unlike std::map thread safety
  is guaranteed in case of simultaneously reading.  Modifying the
  values of different ***existing*** elements is also thread safe.

  The key in this class is usually the global box index.  Note that
  MFIter::index() and FabSetIter.index() return the global box index
  in increasing order.  Thus, it is safe to insert new elements with
  global box index as key inside a loop of MFIter or FabSetIter.
 */

template <class T>
class BLMap
{
public:
    //
    // constructs a BLMap<T> with no elements
    //
    BLMap ();

    int size () const { 
	BL_ASSERT(keys.size() == values.size()); return keys.size(); }
    //
    // Returns a reference to the value at local index (not the key)
    //
    T& at_local (int i) { return values[i]; }
    const T& at_local (int i) const { return values[i]; }
    //
    // Returns local index associated with the key
    //
    int local_index (int k) const;
    //
    // Returns a refence to the value associated with the key
    //
    T& operator[] (int k);
    const T& operator[] (int k) const;
    //
    void reserve (int n) { keys.reserve(n); values.reserve(n); }
private:
    std::vector<int> keys;
    std::vector<T>   values;
};

template <class T>
BLMap<T>::BLMap ()
    : keys(), values()
{}

template <class T>
int BLMap<T>::local_index (int k) const { 
    std::vector<int>::const_iterator low
	= std::lower_bound(keys.begin(), keys.end(), k);
    if (low != keys.end() && *low == k) {
	return low - keys.begin();
    }
    else {
	return -1;
    }
}

template <class T>
T& BLMap<T>::operator[] (int k)
{
    BL_ASSERT(keys.size() == values.size());

    if (keys.empty() || k > keys.back()) {
	keys.push_back(k);
	values.resize(keys.size());
	return values.back();	
    }
    else {
	int i = local_index(k);
	BL_ASSERT(i >= 0);
	return values[i];
    }
}

template <class T>
const T& BLMap<T>::operator[] (int k) const
{
    BL_ASSERT(keys.size() == values.size());
    int i = local_index(k);
    BL_ASSERT(i >= 0);
    return values[i];
}

#endif /*BLMAP_H*/
