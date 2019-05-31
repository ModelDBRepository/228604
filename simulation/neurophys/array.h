#ifndef ARRAY_H
#define	ARRAY_H

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string.h>
#include <vector>

namespace neurophys {

class ArrayElementOutOfBoundsException: public std::exception {};
class ArraySizeMismatchException: public std::exception {};

// the difference between class member operators and external friends might be
// unnoetig here!
template <class T>
class Array {
    size_t size_;
    T* data_;    
public:
    explicit Array(const size_t size=0);
    Array(const Array<T>& other);
    Array(const std::vector<T>& vec);
#ifdef MOVE_SEMANTICS
    Array(Array<T>&& other);
#endif
    virtual ~Array();
    
    Array<T>& operator=(Array<T> other);
    const size_t size() const { return size_; }
    const T* constdata() const { return data_; } // das geht bestimmt eleganter, oder?
    T* data() { return data_; }
    T& operator[](const int idx);
    const T& operator[](const int idx) const;
    Array<T>& operator*=(const T val);
    Array<T>& operator*=(const Array<T>& other);
    Array<T>& operator/=(const T val);
    Array<T>& operator/=(const Array<T>& other);
    Array<T>& operator+=(const T val);
    Array<T>& operator+=(const Array<T>& other);
    Array<T>& operator-=(const T val);
    Array<T>& operator-=(const Array<T>& other);
    Array<T> operator-() const;
    
    template <class X>
    friend void swap(Array<X>& first, Array<X>& second);
    template <class X>
    friend Array<X> operator*(const Array<X>& first, 
                                 const Array<X>& second);
    template <class X, class Y>
    friend Array<X> operator*(const Array<X>& array, const Y value);
    template <class X, class Y>
    friend Array<X> operator*(const Y value, const Array<X>& array);
    template <class X, class Y>
    friend Array<X> operator/(const Array<X>& array, const Y value);
    template <class X, class Y>
    friend Array<X> operator/(const Y value, const Array<X>& array);
    template <class X, class Y>
    friend Array<X> operator+(const Array<X>& array, const Y value);
    template <class X, class Y>
    friend Array<X> operator+(const Y value, const Array<X>& array);
    template <class X, class Y>
    friend Array<X> operator-(const Array<X>& array, const Y value);
    template <class X, class Y>
    friend Array<X> operator-(const Y value, const Array<X>& array);
    
    template <class X>
    friend Array<X> operator/(const Array<X>& first, 
                                 const Array<X>& second);
    template <class X>
    friend Array<X> operator+(const Array<X>& first, 
                                 const Array<X>& second);
    template <class X>
    friend Array<X> operator-(const Array<X>& first, 
                                 const Array<X>& second);
    template <class X>
    friend const bool operator==(const Array<X>& first, 
                                 const Array<X>& second);
};

template <class T>
Array<T>::Array(const size_t size):
    size_(size), data_(size_ ? new T[size_] : 0) 
{
}

// is it a good idea to have this and allow for implicit construction?
template <class T>
Array<T>::Array(const std::vector<T>& vec):
    size_(vec.size()), data_(vec.size() ? new T[size_] : 0)
{
    std::copy(vec.data(), vec.data() + vec.size(), data_);
}
        
        

template <class T>
Array<T>::Array(const Array<T>& other): 
    size_(other.size_), data_(size_ ? new T[size_] : 0)                
{
    //std::copy(other.data_, other.data_ + size_, data_);
    for (int i = 0; i < size_; i++)
    {
        data_[i] = other.data_[i];
    }
}

#ifdef MOVE_SEMANTICS
template <class T>
Array<T>::Array(Array<T>&& other):
    size_(other.size_), data_(other.data_)       
{
    other.data_ = 0;
}
#endif

template <class T>
Array<T>::~Array()
{
    delete[] data_;
}

template <class T>
Array<T>& Array<T>::operator=(Array other)
{
    swap(*this, other);
    return *this;
}


template <class T>
void swap(Array<T>& first, Array<T>& second)
{
    using std::swap;
    swap(first.size_, second.size_);
    swap(first.data_, second.data_);
}

template <class T>
T& Array<T>::operator[](const int idx)
{
#ifndef NDEBUG
    if ((idx < 0) || (idx >= size_)) 
        throw ArrayElementOutOfBoundsException();
#endif // NDEBUG
    return data_[idx];
}

template <class T>
const T& Array<T>::operator[](const int idx) const
{
    return data_[idx];
}

template <class X, class Y>
Array<X> operator*(const Array<X>& array, const Y val)
{
    Array<X> result(array.size_);
    for (int i = 0; i < array.size_; i++)
    {
        result.data_[i] = array.data_[i] * val;
    }
    return result;
}
template <class X, class Y>
Array<X> operator*(const Y val, const Array<X>& array)
{
    return array * val;
}

template <class X, class Y>
Array<X> operator/(const Array<X>& array, const Y val)
{
    return array * (1./val);
}
template <class X, class Y>
Array<X> operator/(const Y val, const Array<X>& array)
{
    Array<X> result(array.size_);
    for (int i = 0; i < array.size_; i++)
    {
        result.data_[i] = val / array.data_[i];
    }
    return result;
}

template <class X, class Y>
Array<X> operator+(const Array<X>& array, const Y val)
{
    Array<X> result(array.size_);
    for (int i = 0; i < array.size_; i++)
    {
        result.data_[i] = array.data_[i] + val;
    }
    return result;
}
template <class X, class Y>
Array<X> operator+(const Y val, const Array<X>& array)
{
    return array + val;
}

template <class X, class Y>
Array<X> operator-(const Array<X>& array, const Y val)
{
    return array + (-val);
}
template <class X, class Y>
Array<X> operator-(const Y val, const Array<X>& array)
{
    Array<X> result(array.size_);
    for (int i = 0; i < array.size_; i++)
    {
        result.data_[i] = val - array.data_[i];
    }
    return result;
}

template <class T>
Array<T> Array<T>::operator-() const
{
    return(*this * (-1));
}

template <class T>
Array<T>& Array<T>::operator*=(const T val)
{
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] * val;
    }
    return *this;
}

template <class T>
Array<T>& Array<T>::operator*=(const Array<T>& other)
{
    if (size_ != other.size_)
       throw ArraySizeMismatchException(); 
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] * other.data_[i];
    }
    return *this;
}

template <class T>
Array<T>& Array<T>::operator+=(const T val)
{
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] + val;
    }
    return *this;
}

template <class T>
Array<T>& Array<T>::operator+=(const Array<T>& other)
{
    if (size_ != other.size_)
       throw ArraySizeMismatchException(); 
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] + other.data_[i];
    }
    return *this;
}

template <class T>
Array<T>& Array<T>::operator/=(const T val)
{
    return operator*=(1./val);
}

template <class T>
Array<T>& Array<T>::operator/=(const Array<T>& other)
{
    if (size_ != other.size_)
       throw ArraySizeMismatchException(); 
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] / other.data_[i];
    }
    return *this;
}

template <class T>
Array<T>& Array<T>::operator-=(const T val)
{
    return operator+=(-val);
}

template <class T>
Array<T>& Array<T>::operator-=(const Array<T>& other)
{
    if (size_ != other.size_)
       throw ArraySizeMismatchException(); 
    for (int i = 0; i < size_; i++)
    {
        data_[i] = data_[i] - other.data_[i];
    }
    return *this;
}

template <class T>
Array<T> operator*(const Array<T>& first, const Array<T>& second)
{
    if (first.size_ != second.size_)
       throw ArraySizeMismatchException(); 
    Array<T> result(first.size_);
    for (int i = 0; i < first.size_; i++)
    {
        result.data_[i] = first.data_[i] * second.data_[i];
    }
    return result;
}

template <class T>
Array<T> operator/(const Array<T>& first, const Array<T>& second)
{
    if (first.size_ != second.size_)
       throw ArraySizeMismatchException(); 
    Array<T> result(first.size_);
    for (int i = 0; i < first.size_; i++)
    {
        result.data_[i] = first.data_[i] / second.data_[i];
    }
    return result;
}

template <class T>
Array<T> operator+(const Array<T>& first, const Array<T>& second)
{
    if (first.size_ != second.size_)
       throw ArraySizeMismatchException(); 
    Array<T> result(first.size_);
    for (int i = 0; i < first.size_; i++)
    {
        result.data_[i] = first.data_[i] + second.data_[i];
    }
    return result;
}

template <class T>
Array<T> operator-(const Array<T>& first, const Array<T>& second)
{
    if (first.size_ != second.size_)
       throw ArraySizeMismatchException(); 
    Array<T> result(first.size_);
    for (int i = 0; i < first.size_; i++)
    {
        result.data_[i] = first.data_[i] - second.data_[i];
    }
    return result;
}

template <class T>
const bool operator==(const Array<T>& first, const Array<T>& second)
{
    if (first.size_ != second.size_)
       throw ArraySizeMismatchException(); 
    for (int i = 0; i < first.size_; i++)
    {
        if (first.data_[i] != second.data_[i])
            return false;
    }
    return true;
}

}
#endif	/* ARRAY_H */

