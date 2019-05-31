#ifndef DYNAMIC_HISTOGRAM_H_
#define DYNAMIC_HISTOGRAM_H_

#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <assert.h>

class HistogramException: public std::exception
{
public:
    HistogramException(const char* msg): _msg(msg) {}
    const char* what() const throw() { return _msg; }
private:
    const char* _msg;
};

class HistogramView
{
public:
    HistogramView(double ll, double d):
        l(ll), dx(d) {}
    HistogramView(): l(0), dx(0) {}
    int n(double x) const { return (dx == 0 ? 0 : (int)floor((x - l) / dx)); }
    double x_f(unsigned int n, double f) const { return (n + f) * dx + l; }
    double x_l(unsigned int n) const { return x_f(n, 0.); }
    double x_m(unsigned int n) const { return x_f(n, 0.5); }
    double x_r(unsigned int n) const { return x_f(n, 1.); }
    /*double r() { return l + w(); }
    double w() { return dx * N; }*/
    double l;
    double dx;

};

template <class T>
class Histogram
{
public:
    static const int NORMALIZED = 1;

    Histogram(unsigned int N, double l, double r):
        _N(N), _v(l, (r-l)/N)
    {
        _data = new T[N + 2];
        memset(_data, 0, (N + 2) * sizeof(T));
    }
    virtual ~Histogram()
    {
        delete[] _data;
    }
    void dump(std::ostream& out, const int flags = 0) const
    {   
        if (flags & NORMALIZED)
        {
            const int c = count();
            for (unsigned int i = 0; i < _N; i++)
            out << _v.x_m(i) << "\t" << (double)_data[i] / c / _v.dx << "\n";
        }
        else
        {
            for (unsigned int i = 0; i < _N; i++)
                out << _v.x_m(i) << "\t" << _data[i] << "\n";
        }
    }
    virtual void feed(double x)
    {
        add_to_bin(_v.n(x), 1);
    }
    unsigned int count() const
    {
        unsigned int c = 0;
        for (unsigned int i = 0; i < _N + 2; i++)
        {
            c += _data[i];
        }
        return c;
    }

    T get_underflow() const { return _data[_N]; }
    T get_overflow() const { return _data[_N + 1]; }

    void add_to_bin(int n, unsigned int howmuch)
    {
        _data[(n < 0)? _N : (n >= _N ? _N + 1 : n)] += howmuch;
    }

    void write_ascii(const char filename[], const int flags = 0) const
    {
        std::ofstream of;
	    of.open(filename, std::ios::out);
        dump(of, flags);
        of.close();
    }

    void write(std::ostream& out) const
    {
        out.write((char*)&_N, sizeof(_N));
        const size_t s = sizeof(T);
        const unsigned int c = count();
        out.write((char*)&s, sizeof(s));
        out.write((char*)&c, sizeof(c));
        out.write((char*)&(_v.l), sizeof(_v.l));
        out.write((char*)&(_v.dx), sizeof(_v.dx));
        out.write((char*)_data, (_N + 2) * sizeof(T));
    }

    void read(std::istream& in)
    {
        unsigned int N;
        size_t s;
        unsigned int c;
        in.read((char*)&N, sizeof(N));
        if (N != _N)
            throw HistogramException("reading hist failed: N does"
                   " not match recorded _N"); 
        in.read((char*)&s, sizeof(s));
        if (sizeof(T) != s)
            throw HistogramException("reading hist failed: sizeof(T) does"
                   " not match recorded size"); 
        in.read((char*)&c, sizeof(c));
        in.read((char*)&(_v.l), sizeof(_v.l));
        in.read((char*)&(_v.dx), sizeof(_v.dx));
        in.read((char*)_data, s * (N + 2));
        if (c != count())
            throw HistogramException("reading hist failed: count does"
                   " not match recorded count");
    }
protected:
    unsigned int _N;
    HistogramView _v;
    T* _data;
};


template <class T>
class DynamicHistogram: public Histogram<T> 
{ 
public: 
    DynamicHistogram(unsigned int N, unsigned int k):
        Histogram<T>(N, 0, 0), _k(k) {}

    virtual void feed(double x);
private:
    void resize(const HistogramView& nv);

    unsigned int _k;
};

template <class T>
void DynamicHistogram<T>::feed(double x)
{
    double m = x;
    double var = x * x;
    unsigned int cnt = 1;
    for (unsigned int i = 0; i < this->_N; i++)
    {
        cnt += this->_data[i];
        m += this->_data[i] * this->_v.x_m(i);
        var += this->_data[i] * this->_v.x_m(i) * this->_v.x_m(i);
    }
    m /= cnt;
    var = var / cnt - m * m; // this is a biased estimator, but we don't care
    double s = sqrt(var);
    
    if (cnt == 1)
        this->_v.l = x;

    if (this->_v.dx * this->_N < this->_k * s)
    {
        const double nw = this->_k * s;
        const double nl = fmin(m - nw /2, this->_v.l);
        const double nr = fmax(m + nw /2, this->_v.l + this->_N * this->_v.dx);
        HistogramView nv(nl, (nr - nl) / this->_N);
        resize(nv);
    }
    
    this->add_to_bin(this->_v.n(x), 1);
}

template <class T>
void DynamicHistogram<T>::resize(const HistogramView& nv) 
{
    T* old_data = new T[this->_N];
    memcpy(old_data, this->_data, sizeof(T) * (this->_N));
    memset(this->_data, 0, sizeof(T) * this->_N); 
                // lass over- u. underflow intakt
    for (unsigned int i = 0; i < this->_N; i++)
    {
        const unsigned int n1 = nv.n(this->_v.x_l(i));
        const unsigned int n2 = nv.n(this->_v.x_r(i));
        
        if (n1 == n2)
            this->add_to_bin(n1, old_data[i]);
        else
        {
            this->add_to_bin(n1, (int)ceil((nv.x_r(n1) - this->_v.x_l(i)) / 
                        this->_v.dx * old_data[i]));
            this->add_to_bin(n2, (int)floor((this->_v.x_r(i) - nv.x_l(n2)) / 
                        this->_v.dx * old_data[i]));
        }
    }
    this->_v.l = nv.l;
    this->_v.dx = nv.dx;
    delete[] old_data;
}
    

#endif // DYNAMIC_HISTOGRAM_H_
