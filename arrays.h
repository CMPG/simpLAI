#ifndef _MY_ARRAYS_HPP_
#define _MY_ARRAYS_HPP_

#include "config.h"
#include <string>

//#include <limits.h>  //used to define MAX_INT
#include <climits>  //used to define MAX_INT
///#include <_null.h>   //definition of NULL //ANSI C++
#ifndef NULL
#define NULL    0
#endif

#ifdef _GCC_
#include <iostream>

#include "UniqueRandom.h"
//#include "public.h"
#endif


using namespace std;

//******************************************************************************
//class TArrayAsVector

template<typename T>
class MY_TArrayAsVector {
private:
   T* _array;
   int _size, _ini_size, _used;
   int _begin, _step;
   void quicksort(int l, int r, T* list); //Loro_09_03_00 Sorts the current array
   void heapsort(int r, T* list);
   int Resize(const int& s); //Helper function to resize the array
   int ResizeCarefully(const int& s); //Helper function to resize the array
public:
   MY_TArrayAsVector();
   MY_TArrayAsVector(const int& s, const int& st);

   MY_TArrayAsVector(const MY_TArrayAsVector<T>& x) {
      if (this != &x) {
         _array = NULL;
         *this = x;
      } 
   }

   ~MY_TArrayAsVector() {
      if (_array) delete[] _array; 
   }

   T & operator[](const int& i) {
      return _array[i];
   }

   const T & operator[](const int& i) const {
      return _array[i];
   }
   MY_TArrayAsVector<T> & operator=(const MY_TArrayAsVector<T>& x);
   MY_TArrayAsVector<T> & operator+=(const MY_TArrayAsVector<T>& x);

   int Add(const T& x) {
      if (_used == _size) {
         Resize(_size + _step);
      }
      _array[_used] = x;
      return ++_used;
   }
   int AddCarefully(const T& x) {
      if (_used == _size) {
         ResizeCarefully(_size + _step);
      }
      _array[_used] = x;
      return ++_used;
   }
   int Find(const T& x) const;

   const int& GetItemsInContainer() const {
      return _used;
   }

   const int& ArraySize() const {
      return _size;
   }
   void Flush();
   void FlushFast();
   int HasMember(const T& t) const;
   void Reallocate(const int& len);
   int size() {return _size;}
   void remove_last(); 
   void quickSortArray(); /*{quicksort(0, _used-1, _array);}*/
   void shuffleArray();

   T* begin() const {
      return &_array[0];
   }

   T* elem(int i) const {
      return & _array[i];
   }

   T* end() const {
      return &_array[_size - 1];
   }
   void setStep(const int & newStep) {this->_step=newStep;}
   
   int elemSize() const {return sizeof(T);}
   void directMemCopy(const MY_TArrayAsVector<T>& x);

};

//******************************************************************************
//class TIArrayAsVector

template<typename T>
class MY_TIArrayAsVector {
private:
   T** _array;
   int _size, _used, _ini_size;
   int _begin, _step;
   T* _out_of_b;
   int _array_owns_elements;
   T** _iterator;
public:
   MY_TIArrayAsVector();
   MY_TIArrayAsVector(const int& s, const int& st);

   MY_TIArrayAsVector(const MY_TIArrayAsVector<T>& x) {
      if (this != &x) * this = x;
      return;
   }
   ~MY_TIArrayAsVector();

   T* & operator[](const int& i) {
      if (i < _used) return _array[i];
      else { 
         if (i < _size) {
            _used = i + 1;
            return _array[i];
         } else return _out_of_b;
      }
   }

   T * operator[](const int& i) const {
      if (i < _used) return _array[i];
      else return _out_of_b;
   }

   T * operator++() {
      return *++_iterator;
   }
   MY_TIArrayAsVector<T> & operator=(const MY_TIArrayAsVector<T>& x);
   int Add(T* x);
   int Find(const T* x);
   void quicksort(int l, int r, T** list);
   void heapsort(int r, T** list);

   const int& GetItemsInContainer() const {
      return _used;
   }

   int IsEmpty() const {
      int found = (_used == 0);
      return found;
   }
   int size() {return _size;}
   void Flush();
   void FlushFast();

   void owns(const int& x) {
      _array_owns_elements = x;
   }
   int owns() const {
      return _array_owns_elements;
   }
   int Detach(T* t);
   int Destroy(T* t);
   int Detach(const int& index);
   int Destroy(const int& index);
   void quickSortArray(); 
   void Reallocate(const int& len);
   void deleteArray();

   T * getIterator() {
      if (_array) {
         _iterator = _array;
         return _array[0];
      }
      return NULL;
   }
};

//******************************************************************************
//class TSArrayAsVector

template<typename T>
class MY_TSArrayAsVector {
private:
   T* _array;
   int _size, _used, _ini_size;
   int _begin, _step;
   T _out_of_b;
public:
   MY_TSArrayAsVector();
   MY_TSArrayAsVector(const int& s, const int& st);

   MY_TSArrayAsVector(const MY_TSArrayAsVector<T>& x) {
      if (this != &x) * this = x;
      return;
   }

   ~MY_TSArrayAsVector() {
      if (_array) delete[] _array;
   }

   T & operator[](const int& i) {
      if (i < _used) return _array[i];
      else { 
         if (i < _size) {
            _used = i + 1;
            return _array[i];
         } else return _out_of_b;
      }
   }

   const T & operator[](const int& i) const {
      if (i < _used) return _array[i];
      return _out_of_b;
   }
   MY_TSArrayAsVector<T> & operator=(const MY_TSArrayAsVector<T>& x);
   int Add(const T& x);
   int Find(const T& x);
   void Reallocate(const int& len);

   const int& GetItemsInContainer() const {
      return _used;
   }
   void Flush();
   void FlushFast();
};

//******************************************************************************
//class TISArrayAsVector

template<typename T>
class MY_TISArrayAsVector {
private:
   T** _array;
   int _size, _used, _ini_size;
   int _begin, _step;
   T* _out_of_b;
   int _array_owns_elements;
public:
   MY_TISArrayAsVector();
   MY_TISArrayAsVector(const int& s, const int& st);

   MY_TISArrayAsVector(const MY_TISArrayAsVector<T>& x) {
      if (this != &x) * this = x;
      return;
   }
   ~MY_TISArrayAsVector();

   T* & operator[](const int& i) {
      if (i < _used) return _array[i];
      else { 
         if (i < _size) {
            _used = i + 1;
            return _array[i];
         } else return _out_of_b;
      }
   }

   T * operator[](const int& i) const {
      if (i < _used) return _array[i];
      return _out_of_b;
   }
   MY_TISArrayAsVector<T> & operator=(const MY_TISArrayAsVector<T>& x);
   int Add(T* x);
   int Find(const T* x);

   const int& GetItemsInContainer() const {
      return _used;
   }
   void Flush();

   void owns(const int& x) {
      _array_owns_elements = x;
   }

   int IsEmpty() const {
      return (_used == 0);
   }
};
//------------------------------------------------------------------------------
//member functions
//------------------------------------------------------------------------------

template<typename T>
MY_TArrayAsVector<T>::MY_TArrayAsVector() {
   _size = 11;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = 10;
   try {
      _array = NULL;
      _array = new T[_size];
   } catch (...) {
      if (_array) delete[] _array;
      return;
   }
}

template<typename T>
MY_TArrayAsVector<T>::MY_TArrayAsVector(const int& s, const int& st) {
   _size = s + 1;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = st;
   try {
      _array = new T[_size];
   } catch (...) {
      if (_array) delete[] _array;
      return;
   }
}

template<typename T>
   MY_TArrayAsVector<T>& MY_TArrayAsVector<T>::operator=(const MY_TArrayAsVector<T>& x) {
   if (_array) delete[] _array;
   _size = x._size;
   _used = x._used;
   _begin = x._begin;
   _step = x._step;
   _array = new T[_size];
   //Safer than memcpy
   for (int i=0; i<_used; ++i) {
      _array[i]=x._array[i];
   }
   return *this;
}

//Loro_21_04_15: Should only be used for basic T types, int, long, double, without any pointers.
template<typename T>
void MY_TArrayAsVector<T>::directMemCopy(const MY_TArrayAsVector<T>& x) {
   if (_array && _size!=x._size) {
      delete[] _array;
      _size = x._size;
      _array = new T[_size];
   }  
   else if (!_array) {
      _size = x._size;
      _array = new T[_size];
   }
   _used = x._used;
   _begin = x._begin;
   _step = x._step;
   memcpy(_array, x._array, _size*sizeof(T));
   return;
}

template<typename T>
   MY_TArrayAsVector<T>& MY_TArrayAsVector<T>::operator+=(const MY_TArrayAsVector<T>& x) {
   int totUsed=_used+x._used;
//   int oldSize=_size;
   if (totUsed>=_size) this->ResizeCarefully(totUsed+_step);
   for (int i=0; i<x._used; ++i) {
      _array[_used++]=x._array[i];
   }
   return *this;
}

template<typename T>
int MY_TArrayAsVector<T>::Resize(const int& newSize) {
   T* temp_copy;
//   try {
      temp_copy = NULL;
      temp_copy = new T[newSize];
//   } catch (...) {
//      if (temp_copy) delete[] temp_copy;
//      return 0;
//   }
   memcpy(temp_copy,_array, sizeof (T) * _size);
   if (_array) delete[] _array;
   _array = temp_copy;
   _size = newSize;
   return _size;
}

//Loro_19_09_10 Procedure like ReSize but without memcpy
template<typename T>
int MY_TArrayAsVector<T>::ResizeCarefully(const int& newSize) {
   T* temp_copy;
   try {
      temp_copy = NULL;
      temp_copy = new T[newSize];
   } catch (...) {
      if (temp_copy) delete[] temp_copy;
      return 0;
   }
#pragma omp simd
   for (int i = 0; i < _size; i++) temp_copy[i] = _array[i];
   if (_array) delete[] _array;
   _array = temp_copy;
   _size = newSize;
   return _size;
}

template<typename T>
int MY_TArrayAsVector<T>::Find(const T& x) const {
   for (int i = 0; i < _used; i++) {
      if (_array[i] == x) return i;
   }
   //case nothing found
   return  INT_MAX;
}

template<typename T>
int MY_TArrayAsVector<T>::HasMember(const T& t) const {
   int index = Find(t);
   if (index < INT_MAX) return 1;
   else return 0;
}

template<typename T>
void MY_TArrayAsVector<T>::Flush() {
   if (_array) delete [] _array;
   try {
      _array = new T[_ini_size];
   } catch (...) {
      delete[] _array;
      _size = 0;
      _array = NULL;
      return;
   }
   _size = _ini_size;
   _used = 0;
   _begin = 0;
}

template<typename T>
void MY_TArrayAsVector<T>::FlushFast() {
   _used = 0;
   _begin = 0;
}

template<typename T>
void MY_TArrayAsVector<T>::Reallocate(const int& len) {
   if (_array) delete [] _array;
   try {
      _array = NULL;
      _array = new T[len];
   } catch (...) {
      if (_array) delete[] _array;
      _size = 0;
      _array = NULL;
      return;
   }
   _size = len;
   _used = 0;
   _begin = 0;
}

template<typename T>
void MY_TArrayAsVector<T>::remove_last() {
   if (_used > 0) _used--;
   return;
}

template<typename T> //BUGGED
void
MY_TArrayAsVector<T>::quicksort(int first, int last, T* data) {
   int i, j;
   T temp, Value;

   if (last > first) {
      /* Divide the list in two */
      Value = data[last];
      i = first - 1;
      j = last;
      do {
         while (data[++i] < Value);
         while (data[--j] > Value);
         temp = data[i];
         data[i] = data[j];
         data[j] = temp;
      } while (j > i);
      data[j] = data [i];
      data[i] = data[last];
      data[last] = temp;
      /* Sort the first part */
      quicksort(first, i - 1, data);
      /* Sort the last part */
      quicksort(i + 1, last, data);
   }

}

template<typename T>
void
MY_TArrayAsVector<T>::quickSortArray() {
   //search the smallest element:
   T temp = _array[0];
   int s_index = 0;

   for (int i = 1; i < _used; i++) {
      if (_array[i] < temp) {
         s_index = i;
         temp = _array[i];
      }
   }
   //put smallest at pos 0
   temp = _array[0];
   _array[0] = _array[s_index];
   _array[s_index] = temp;

   heapsort(_used - 1, _array); //sort 1 to _used-1
}
        
template<typename T>
void
MY_TArrayAsVector<T>::heapsort(int n, T* ra) {
   int i, ir, j, l;
   T rra;

   if (n < 2) return;
   l = (n >> 1) + 1;
   ir = n;
   for (;;) {
      if (l > 1) {
         rra = ra[--l];
      } else {
         rra = ra[ir];
         ra[ir] = ra[1];
         if (--ir == 1) {
            ra[1] = rra;
            break;
         }
      }
      i = l;
      j = l + l;
      while (j <= ir) {
         if (j < ir && ra[j] < ra[j + 1]) j++;
         if (rra < ra[j]) {
            ra[i] = ra[j];
            i = j;
            j <<= 1;
         } else j = ir + 1;
      }
      ra[i] = rra;
   }
}

/* (C) Copr. 1986-92 Numerical Recipes Software %8#1&&)L. */

//******************************************************************************
//class TIArrayAsVector

template<typename T>
MY_TIArrayAsVector<T>::MY_TIArrayAsVector() {
   _size = 11;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = 10;
   try {
      _array = new T*[_size];
   } catch (...) {
      deleteArray();
      return;
   }
   _array_owns_elements = 1;
}

template<typename T>
MY_TIArrayAsVector<T>::MY_TIArrayAsVector(const int& s, const int& st) {
   _size = s + 1;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = st;
   try {
      _array = new T*[_size];
   } catch (...) {
      _array = NULL;
      return;
   }
   _array_owns_elements = 1;
}

template<typename T>
MY_TIArrayAsVector<T>::~MY_TIArrayAsVector() {
   deleteArray();
}

template<typename T>
   MY_TIArrayAsVector<T>& MY_TIArrayAsVector<T>::operator=(const MY_TIArrayAsVector<T>& x) {
   //Loro_31_03_14
   deleteArray();
   _size = x._size;
   _used = x._used;
   _begin = x._begin;
   _step = x._step;
   _out_of_b = x._out_of_b;
   try {
      _array = new T*[_size];
      for (int i = 0; i < _used; i++) _array[i] = new T(*(x._array[i]));
   } catch (...) {
      deleteArray();
      return *this;
   }
   return *this;
}

template<typename T>
int MY_TIArrayAsVector<T>::Add(T* x) {
   if (_used == _size) {
      int new_size = _size + _step;
      T** cop;
      try {
         cop = new T*[new_size];
      } catch (...) {
         deleteArray();
         return 0;
      }
      //Loro_01_04_14 potential bug, reclaiming too much memory...
      //memcpy(cop, _array, sizeof (T) * _size);
      memcpy(cop, _array, sizeof (T*) * _size);
      delete [] _array;
      _array = cop;
      _size = new_size;
   }
   _array[_used] = x;
   _used++;
   return _used;
}

template<typename T>
void MY_TIArrayAsVector<T>::Reallocate(const int& len) {
   //Loro_31_03_14
   deleteArray();
   try {
      _array = NULL;
      _array = new T*[len];
      //Loro_10_01_14 safer
      memset(_array, 0, sizeof(T*)*len);
   } catch (...) {
      //Loro_31_03_14
      deleteArray();
      return;
   }
   _size = len;
   _used = 0;
   _begin = 0;
}

template<typename T>
int MY_TIArrayAsVector<T>::Find(const T* x) {
   for (int i = 0; i < _used; i++) {
      if (*(_array[i]) == *x) return i;
   }
   //case nothing found
   return INT_MAX;
}

template<typename T>
void MY_TIArrayAsVector<T>::Flush() {
   //Loro_31_03_14
   deleteArray();
   try {
      _array = new T*[_ini_size];
   } catch (...) {
      deleteArray();
      return;
   }
   _size = _ini_size;
   _used = 0;
   _begin = 0;
}

//Loro_31_03_14
template<typename T>
void MY_TIArrayAsVector<T>::deleteArray() {
   if (_array_owns_elements) {
      for (int i = 0; i < _used; i++) {
         if (_array[i]) delete _array[i];
         _array[i] = NULL;
      }
   }
   delete [] _array;
   _array = NULL;
   _size = 0;
   _used = 0;
   _begin = 0;
}

template<typename T>
void MY_TIArrayAsVector<T>::FlushFast() {
   if (_array_owns_elements) {
      for (int i = 0; i < _used; i++) {
         if (_array[i]) delete _array[i];
         _array[i] = NULL;
      }
   }
   _used = 0;
   _begin = 0;
}

template<typename T> //BUGGED
void
MY_TIArrayAsVector<T>::quicksort(int first, int last, T** data) {
   int i, j;
   T* temp, *Value;
   if (last > first) {
      /* Divide the list in two */
      Value = data[last];
      i = first - 1;
      j = last;
      do {
         while ((*data[++i]) < (*Value));
         while ((*data[--j]) > (*Value));
         temp = data[i];
         data[i] = data[j];
         data[j] = temp;
      } while (j > i);
      data[j] = data [i];
      data[i] = data[last];
      data[last] = temp;
      /* Sort the first part */
      quicksort(first, i - 1, data);
      /* Sort the last part */
      quicksort(i + 1, last, data);
   }
}

template<typename T>
void
MY_TIArrayAsVector<T>::heapsort(int n, T** ra) {
   int i, ir, j, l;
   T* rra;

   if (n < 2) return;
   l = (n >> 1) + 1;
   ir = n;
   for (;;) {
      if (l > 1) {
         rra = ra[--l];
      } else {
         rra = ra[ir];
         ra[ir] = ra[1];
         if (--ir == 1) {
            ra[1] = rra;
            break;
         }
      }
      i = l;
      j = l + l;
      while (j <= ir) {
         if (j < ir && (*ra[j]) < (*ra[j + 1])) j++;
         if ((*rra) < (*ra[j])) {
            ra[i] = ra[j];
            i = j;
            j <<= 1;
         } else j = ir + 1;
      }
      ra[i] = rra;
   }
}

/* (C) Copr. 1986-92 Numerical Recipes Software %8#1&&)L. */


template<typename T>
void
MY_TIArrayAsVector<T>::quickSortArray() {
   //search the smallest element:
   T* temp = _array[0];
   int s_index = 0;

   for (int i = 1; i < _used; i++) {
      if ((*_array[i]) < (*temp)) {
         s_index = i;
         temp = _array[i];
      }
   }
   //put smallest at pos 0
   temp = _array[0];
   _array[0] = _array[s_index];
   _array[s_index] = temp;

   heapsort(_used - 1, _array); //sort 1 to _used-1


}

template<typename T>
int MY_TIArrayAsVector<T>::Detach(T* t) { //not tested
   int index = Find(t);
   if (index < INT_MAX) {
      if (_array_owns_elements) if (_array[index]) delete _array[index];
      for (int i = index; i < _used - 1; ++i) { //we shift the upper objects
         _array[i] = _array[i + 1];
      }
      _used--;
   } else {
      return 0;
   }
   return 1;
}

template<typename T>
int MY_TIArrayAsVector<T>::Destroy(T* t) {

   int index = Find(t);
   if (index < INT_MAX) {
      if (_array_owns_elements) if (_array[index]) delete _array[index]; //we delete the object
      for (int i = index; i < _used - 1; i++) { //we shift the upper objects
         _array[i] = _array[i + 1];
      }
      _used--;
   } else {
      return 0;
   }
   return 1;
}

template<typename T>
int MY_TIArrayAsVector<T>::Detach(const int& index) { //not tested
   if (index < INT_MAX) {
      //Loro_13_07_00
      if (_array_owns_elements) delete _array[index];
      for (int i = index; i < _used - 1; i++) { //we shift the upper objects
         _array[i] = _array[i + 1];
      }
      _used--;
   } else {
      return 0;
   }
   return 1;
}

template<typename T>
int MY_TIArrayAsVector<T>::Destroy(const int& index) {
   if (index < INT_MAX) {
      if (_array[index]) delete _array[index]; //we delete the object
      for (int i = index; i < _used - 1; i++) { //we shift the upper objects
         _array[i] = _array[i + 1];
      }
      _used--;
   } else {
      return 0;
   }
   return 1;
}




//******************************************************************************
//class TSArrayAsVector

template<typename T>
MY_TSArrayAsVector<T>::MY_TSArrayAsVector() {
   _size = 11;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = 10;
   try {
      _array = new T[_size];
   } catch (...) {
      _array = NULL;
      return;
   }
}

template<typename T>
MY_TSArrayAsVector<T>::MY_TSArrayAsVector(const int& s, const int& st) {
   _size = s + 1;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = st;
   try {
      _array = new T[_size];
   } catch (...) {
      _array = NULL;
      return;
   }
}

//template<typename T>
//   MY_TSArrayAsVector<T>& MY_TSArrayAsVector<T>::operator=(const MY_TSArrayAsVector<T>& x) {
//   Flush();
//   if (_array) delete[] _array;
//   _size = x._size;
//   _used = x._used;
//   _begin = x._begin;
//   _step = x._step;
//   _out_of_b = x._out_of_b;
//   try {
//      _array = new T[_size];
//   } catch (...) {
//      delete[] _array;
//      _size = 0;
//      _array = NULL;
//      return *this;
//   }
//   //Loro_20_07_10
//   memcpy(_array, x._array, sizeof (T) * _size);
//   //for (int i = 0; i < _used; i++) _array[i] = x._array[i];
//   return *this;
//}

template<typename T>
int MY_TSArrayAsVector<T>::Add(const T& x) {
   if (_used == _size) {
      int new_size = _size + _step;
      T* cop;
      try {
         cop = new T[new_size];
      } catch (...) {
         delete[] _array;
         _size = 0;
         _array = NULL;
         return 0;
      }

      for (int i = 0; i < _size; i++) cop[i] = _array[i];
      delete[] _array;
      _array = cop;
      _size = new_size;
   }
   _array[_used] = x;
   _used++;
   //now search the right place to insert the element: for the moment it is
   //inserted at the end of the array. We do now some swaps to put it a the
   //right place:
   T temp;
   int i;
   for (i = _used - 2; i >= 0; i--) {
      if (_array[i + 1] < _array[i]) {
         temp = _array[i];
         _array[i] = _array[i + 1];
         _array[i + 1] = temp;
      } else {
         break;
      }
   }
   return (i + 1); //stef_28_3 //the element was inserted at position (i+1)
}

template<typename T>
int MY_TSArrayAsVector<T>::Find(const T& x) {
   for (int i = 0; i < _used; i++) {
      if (_array[i] == x) return i;
   }
   //case nothing found
   return  INT_MAX;

}

template<typename T>
void MY_TSArrayAsVector<T>::Flush() {
   delete [] _array;
   try {
      _array = new T[_ini_size];
   } catch (...) {
      delete[] _array;
      _size = 0;
      _array = NULL;
      return;
   }
   _size = _ini_size;
   _used = 0;
   _begin = 0;
}
template<typename T>
void MY_TSArrayAsVector<T>::FlushFast() {
   _used = 0;
   _begin = 0;
}

template<typename T>
void MY_TSArrayAsVector<T>::Reallocate(const int& len) {
   if (_array) delete [] _array;
   try {
      _array = NULL;
      _array = new T[len];
   } catch (...) {
      if (_array) delete[] _array;
      _size = 0;
      _array = NULL;
      return;
   }
   _size = len;
   _used = 0;
   _begin = 0;
}

//******************************************************************************
//class TIArrayAsVector

template<typename T>
MY_TISArrayAsVector<T>::MY_TISArrayAsVector() {
   _size = 11;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = 10;
   try {
      _array = new T*[_size];
   } catch (...) {
      _array = NULL;
      return;
   }
   _array_owns_elements = 1;
}

template<typename T>
MY_TISArrayAsVector<T>::MY_TISArrayAsVector(const int& s, /*const int& b,*/ const int& st) {
   _size = s + 1;
   _ini_size = _size;
   _used = 0;
   _begin = 0;
   _step = st;
   try {
      _array = new T*[_size];
   } catch (...) {
      _array = NULL;
      return;
   }
   _array_owns_elements = 1;
}

template<typename T>
MY_TISArrayAsVector<T>::~MY_TISArrayAsVector() {
   if (_array_owns_elements) {
      for (int i = 0; i < _used; i++) {
         if (_array[i]) delete _array[i];
         _array[i] = NULL;
      }
   }
   if (_array) delete[] _array;
}

/*
template<typename T>
T*& MY_TISArrayAsVector<T>::operator[](const int& i)
  {
        if(_size<i)
                return _out_of_b;
        else
                return _array[i];

}

template<typename T>
T* MY_TISArrayAsVector<T>::operator[](const int& i) const
  {
        if(_size<i)
                return _out_of_b;
        else
                return _array[i];

}
 */
template<typename T>
   MY_TISArrayAsVector<T>& MY_TISArrayAsVector<T>::operator=(const MY_TISArrayAsVector<T>& x) {
   Flush();
   if (_array) delete[] _array;
   _size = x._size;
   _used = x._used;
   _begin = x._begin;
   _step = x._step;
   _out_of_b = x._out_of_b;
   try {
      _array = new T*[_size];
      for (int i = 0; i < _used; i++) _array[i] = new T(*(x._array[i]));
   } catch (...) {
      delete[] _array;
      _size = 0;
      _array = NULL;
      return *this;
   }

   return *this;
}

template<typename T>
int MY_TISArrayAsVector<T>::Add(T* x) {
   if (_used == _size) {
      int new_size = _size + _step;
      T** cop;
      try {
         cop = new T*[new_size];
      } catch (...) {
         delete[] _array;
         _size = 0;
         _array = NULL;
         return 0;
      }
      for (int i = 0; i < _size; i++) cop[i] = _array[i];
      delete[] _array;
      _array = cop;
      _size = new_size;
   }
   _array[_used] = x;
   _used++;
   //now search the right place to insert the element: for the moment it is
   //inserted at the end of the array. We do now some swaps to put it a the
   //right place:
   T* temp;
   int i;
   for (i = _used - 2; i >= 0; i--) {
      if (*(_array[i + 1]) < *(_array[i])) {
         temp = _array[i];
         _array[i] = _array[i + 1];
         _array[i + 1] = temp;
      } else {
         break;
      }
   }
   return (i + 1); //stef_28_3_98 //return the position Where inserted
}

template<typename T>
int MY_TISArrayAsVector<T>::Find(const T* x) {
   for (int i = 0; i < _used; i++) {
      if ((*(_array[i])) == (*x)) return i;
   }
   //case nothing found
   return INT_MAX;

}

template<typename T>
void MY_TISArrayAsVector<T>::Flush() {
   if (_array_owns_elements) {
      for (int i = 0; i < _used; i++) {
         if (_array[i]) delete _array[i];
         _array[i] = NULL;
      }
   }
   delete [] _array;
   try {
      _array = new T*[_ini_size];
   } catch (...) {
      delete[] _array;
      _size = 0;
      _array = NULL;
      return;
   }

   _size = _ini_size;
   _used = 0;
   _begin = 0;
   //Loro_20_07_00
   //_step=10;
}

template<typename T>
   MY_TSArrayAsVector<T>& MY_TSArrayAsVector<T>::operator=(const MY_TSArrayAsVector<T>& x) {
   if (_array) delete[] _array;
   _size = x._size;
   _used = x._used;
   _begin = x._begin;
   _step = x._step;
   _out_of_b = x._out_of_b;
   _array = new T[_size];
   //Loro_19_12_13 Caused memory leak ... 
//   memcpy(_array, x._array, sizeof (T) * _size);
   for (int i=0; i<x.GetItemsInContainer(); ++i) {
      _array[i]=x._array[i];
   }
   return *this;
}

#endif




