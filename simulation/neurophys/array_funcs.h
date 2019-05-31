/* 
 * File:   array_funcs.h
 * Author: fedro
 *
 * Created on December 10, 2014, 3:25 PM
 */

#ifndef ARRAY_FUNCS_H
#define	ARRAY_FUNCS_H

#include "array.h"
#include <assert.h>
#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <complex>
#include <fstream>

namespace neurophys {
namespace arr {

    
    
template<class T>
Array<T> zeros(const size_t size)
{
    // this assumes that whatever the type T, setting the memory byte-wise to
    // zero means the elements become 0
    
    Array<T> result(size);
    memset(result.data(), 0, size*sizeof(T));
    return result;    
}

template<class T>
Array<T> ones(const size_t size)
{
    Array<T> result(size);
    for (int i = 0; i < size; i++)
        result[i] = 1.;
    return result;
}

template<class T>
const T sum(const Array<T>& ar)
{
    T res = 0.;
    for (int i = 0; i < ar.size(); i++)
        res += ar[i];
    return res;
}

template<class T>
const T max(const Array<T>& ar)
{
    T max = ar[0];
    for (int i = 0; i < ar.size(); i++)
        if (ar[i] > max)
            max = ar[i];
    return max;
}

template <class T>
void columns_to_stream(std::ostream& out, const Array<T>& x, const Array<T>& y)
{
    /*if (x.size() != y.size())
        throw ArraySizeMismatchException();*/
    for (int i = 0; i < std::min(x.size(),y.size()); i++)
        out << x[i] << "\t" << y[i] << "\n";
}

template <class T>
void columns_to_stream(std::ostream& out, const Array<T>& x, const Array<T>& y,
                        const Array<T>& z)
{
    /*if (x.size() != y.size())
        throw ArraySizeMismatchException();*/
    for (int i = 0; i < std::min(x.size(),y.size()); i++)
        out << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
}

template <class T>
void columns_to_file(const char filename[], const Array<T>& x, const Array<T>& y)
{
    std::ofstream of; 
    of.open(filename, std::ios::out);
    columns_to_stream(of, x, y);
    of.close();
}

template <class T>
void columns_to_file(const char filename[], const Array<T>& x, const Array<T>& y,
                        const Array<T>& z)
{
    std::ofstream of; 
    of.open(filename, std::ios::out);
    columns_to_stream(of, x, y, z);
    of.close();
}

template <class T> 
Array<T> downsample(const Array<T>& src, const size_t destsize)
{
    assert(destsize <= src.size());
    Array<T> result = zeros<double>(destsize);
    const double f = static_cast<double>(destsize)/src.size();
    for (int i = 0; i < src.size(); i++)
    {
        result[i * f] += src[i] * f;
    }
    return result;
}

double mean(const Array<double>& ar);
double var(const Array<double>& ar);
Array<double> range(const double a, const double b, const double inc=1);
Array<double> cos(const Array<double>& x);
Array<double> sin(const Array<double>& x);
Array<std::complex<double> > conj(const Array<std::complex<double> >& x);
Array<double> real(const Array<std::complex<double> >& x);
Array<double> imag(const Array<std::complex<double> >& x);
Array<double> abs(const Array<std::complex<double> >& x);

}}

#endif	/* ARRAY_FUNCS_H */

