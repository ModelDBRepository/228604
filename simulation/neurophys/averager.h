#ifndef AVERAGER_H_
#define AVERAGER_H_

#include <gsl/gsl_statistics_double.h>
#include <string.h>
#include <iostream>

class Averager
{
public:
    virtual ~Averager() {};
    virtual void reset() = 0;
    virtual void feed(const double x) = 0;
    virtual double get_mean() const = 0;
    virtual double get_var() const = 0;
    virtual unsigned int get_n() const = 0;
};

class SimpleAverager: public Averager
{
public:
    SimpleAverager():
        _n(0), _sum(0), _sum_of_squares(0) {}
    virtual void reset() { _n = 0; _sum = 0; _sum_of_squares = 0; }
    virtual void feed(const double x) 
    { 
        _n++; _sum += x; 
        _sum_of_squares += x * x; 
    }
    virtual double get_mean() const { return _sum / _n; }
    virtual double get_var() const 
    { 
        return 1./(_n - 1) * (_sum_of_squares - _n * get_mean()*get_mean()); 
    }
    virtual unsigned int get_n() const { return _n; }
private:
    unsigned int _n;
    double _sum;
    double _sum_of_squares;
};

class BufferAverager: public Averager
{
public:
    BufferAverager(const size_t bufsize):
        bufsize_(bufsize), n_(0), i_(0)
    {
        buf_ = new double[bufsize];
    }
    virtual ~BufferAverager()
    {
        delete[] buf_;
    }
    virtual void feed(const double x)
    {
        if (i_ == bufsize_) 
        {
            i_ = 0;
        }
        if (i_ < bufsize_)
        {
            buf_[i_] = x;
            i_++;
            if (i_ > n_)
                n_++;
        }
    }
    virtual double get_mean() const 
    {
        return gsl_stats_mean(buf_, 1, n_);
    } 
    virtual double get_var() const 
    {
        return gsl_stats_variance(buf_, 1, n_);
    }
    virtual unsigned int get_n() const 
    {
        return n_;
    } 
    virtual void reset()
    {
        memset(buf_, 0, bufsize_ * sizeof(double));
        n_ = 0;
        i_ = 0;
    }
    void dump(std::ostream& out)
    {
        for (int j = 0; j < n_; j++)
            out << buf_[j] << "\n";
    }

private:
    const size_t bufsize_;
    double* buf_;
    unsigned int n_;
    unsigned int i_;
};

class FancyAverager: public Averager
{
public:
    FancyAverager(): 
        _mean_avg(), _sq_avg(), _cub_avg() {}
    virtual void reset()
    {
        _mean_avg.reset();
        _sq_avg.reset();
    }
    virtual void feed(const double x)
    {
        _mean_avg.feed(x);
        _sq_avg.feed(x*x);
        _cub_avg.feed(x*x*x);
    }
    virtual unsigned int get_n() const { return _mean_avg.get_n(); }
    virtual double get_mean() const { return _mean_avg.get_mean(); }
    virtual double get_var() const { return _mean_avg.get_var(); }
    double get_mean_of_square() const { return _sq_avg.get_mean(); }
    double get_var_of_square() const { return _sq_avg.get_var(); }
    double get_mean_of_cube() const { return _cub_avg.get_mean(); }
    double get_var_of_cube() const { return _cub_avg.get_var(); }
private:
    SimpleAverager _mean_avg;
    SimpleAverager _sq_avg;
    SimpleAverager _cub_avg;
};

#endif // AVERAGER_H_
