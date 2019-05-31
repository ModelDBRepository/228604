#ifndef FANO_H_
#define FANO_H_

#include <math.h>
#include "averager.h"
#include "array.h"

class FanoFactorEstimator
{
public:
    FanoFactorEstimator(const size_t N, const double minexp, 
            const double maxexp):
        N_(N), ts_(0), avgs_(0), n_trials_(1), icur_(0), last_t_(0.)
    {
        ts_ = new double[N_];
        avgs_ = new SimpleAverager[N_];
        for (unsigned int i = 0; i < N_; i++)
        {
            ts_[i] = pow(10, minexp+(maxexp-minexp)*i/(N_-1));
        }
    }
    virtual ~FanoFactorEstimator()
    {
        delete[] ts_;
        delete[] avgs_;
    }
    void feed(const double t, const unsigned int count)
    {
        if (t < last_t_) // seems a new trial has begun
        {
            n_trials_++;
            icur_ = 0;
            last_t_ = 0;
        }
        if (icur_ < N_ && t > ts_[icur_] && avgs_[icur_].get_n() < n_trials_)
        {
            avgs_[icur_].feed(count);
            icur_++;
            last_t_ = t;
        }
    }
    void dump(std::ostream& out) const
    {   
        for (unsigned int i = 0; i < N_; i++)
        {
            out << ts_[i] << "\t" << 
                avgs_[i].get_var()/avgs_[i].get_mean() << "\n";
        }    
    }
private:
    const size_t N_;
    double* ts_;
    SimpleAverager* avgs_;
    unsigned int n_trials_;
    unsigned int icur_;
    double last_t_;
};

#endif // FANO_H
