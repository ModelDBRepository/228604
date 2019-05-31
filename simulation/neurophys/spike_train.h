#ifndef SPIKE_TRAIN_H_
#define SPIKE_TRAIN_H_

#include <vector>
#include <iostream>
#include "rng.h"
#include "array.h"

namespace neurophys
{
           
struct Spike
{
    Spike(const double t, const double a):
        time(t), amplitude(a) {}
    Spike(const double t):
        time(t), amplitude(1.) {}
    
    double time;
    double amplitude;
};

class SpikeTrain
{
public:
    SpikeTrain(): T_(-1) {}
    SpikeTrain(const double T): T_(T) {}
    SpikeTrain(const std::vector<double>& spike_times): T_(-1)
    {
        for (std::vector<double>::const_iterator it = spike_times.begin();
                it != spike_times.end(); ++it)
        {
            push_back(Spike(*it));
        }
    }
    
    typedef std::vector<Spike>::iterator iterator;
    typedef std::vector<Spike>::const_iterator const_iterator;
    
    void clear() { spikes_.clear(); }       
    void push_back(Spike sp) 
    {
        if (spikes_.size() == 0 || sp.time >= spikes_.back().time)
            spikes_.push_back(sp);
    }
    size_t size() const { return spikes_.size(); }
    const Spike& front() const { return spikes_.front(); }
    Spike& front() { return spikes_.front(); }
    const Spike& back() const { return spikes_.back(); }
    Spike& back() { return spikes_.back(); }
    iterator begin() { return spikes_.begin(); }
    iterator end() { return spikes_.end(); }
    iterator first_spike_after(const double t);
    iterator last_spike_before(const double t);
    const_iterator begin() const { return spikes_.begin(); }
    const_iterator end() const { return spikes_.end(); }
    const_iterator first_spike_after(const double t) const;
    const_iterator last_spike_before(const double t) const;
    
    double getT() const { return T_ > 0 ? T_ : (size() > 0 ? back().time : 0.); }
private:
    double T_;
    std::vector<Spike> spikes_;
};

namespace spike_train
{
    void to_stream(const SpikeTrain& st, std::ostream& out);
    void from_stream(std::istream& in, SpikeTrain& st);
    void times_from_stream(std::istream& in, SpikeTrain& st);
    
    void from_file(const char* filename, SpikeTrain& st);
    
    Array<double> to_array(const SpikeTrain& st, const size_t N, const double T);
    Array<double> to_array(const SpikeTrain& st, const double tfrom, const double tto, 
            const size_t N);
    SpikeTrain generate_poisson(RNG& rng, const double T, const double r0, const double t0=0.);
    SpikeTrain generate_inhomogeneous_poisson(RNG& rng, const double T, 
            const Array<double>& r, const double rmax, const double t0=0.);
};

};
#endif // SPIKE_TRAIN_H_
