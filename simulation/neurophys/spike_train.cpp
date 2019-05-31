#include "spike_train.h"
#include <string.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include "array_funcs.h"

namespace neurophys {

bool spiketimecmp(const Spike s1, const Spike s2)
{
    return s1.time < s2.time;
}
    
SpikeTrain::iterator SpikeTrain::first_spike_after(const double t) 
{
    return std::upper_bound(spikes_.begin(), spikes_.end(), Spike(t), spiketimecmp);
}

SpikeTrain::iterator SpikeTrain::last_spike_before(const double t) 
{
    SpikeTrain::iterator it = first_spike_after(t);
    if (it == begin()) return it;
    else return --it;
}

// const correctnes schoen und gut, aber das hier sollte doch irgendwie eleganter gehen?

SpikeTrain::const_iterator SpikeTrain::first_spike_after(const double t) const
{
    return std::upper_bound(spikes_.begin(), spikes_.end(), Spike(t), spiketimecmp);
}

SpikeTrain::const_iterator SpikeTrain::last_spike_before(const double t) const
{
    SpikeTrain::const_iterator it = first_spike_after(t);
    if (it == begin()) return it;
    else return --it;
}
    
namespace spike_train {
    void to_stream(const SpikeTrain& st, std::ostream& out)
    {
        for (SpikeTrain::const_iterator it = st.begin(); it != st.end(); it++)
        {
            out << it->time << "\t" << it->amplitude << "\n";
        }
    }
    
    void skip_comments(std::istream& in)
    {
        while (in.peek() == '#')
        {
            in.ignore(1024, '\n');
        }
    }
    
    void from_stream(std::istream& in, SpikeTrain& st)
    {
        st.clear();
        double t = 0;
        double a = 0;
        skip_comments(in);                
        while (in >> t >> a) // this c++ stream syntax really is ugly!
        {
            st.push_back(Spike(t, a));
        }
    }
    
    void times_from_stream(std::istream& in, SpikeTrain& st)
    {
        st.clear();
        double t = 0;
        skip_comments(in);
        while (in >> t) 
        {
            st.push_back(Spike(t));
        }
    }

    Array<double> to_array(const SpikeTrain& st, const size_t N, const double T)
    {
        return to_array(st, 0, T, N);
    }
    
    Array<double> to_array(const SpikeTrain& st, const double tfrom, const double tto, 
            const size_t N)
    {
        Array<double> ar = arr::zeros<double>(N);
        
        const double dt = (tto - tfrom)/N;
        int i = 0;
        for (SpikeTrain::const_iterator it = st.first_spike_after(tfrom); 
                it != st.first_spike_after(tto); it++)
        {
            if ((i = (int)((it->time - tfrom) / dt)) < N)
            {
                ar[i] += it->amplitude/dt;
            }
        }
        
        return ar;
    }

    SpikeTrain generate_poisson(RNG& rng, const double T, const double r0,
            const double t0)
    {
        SpikeTrain st(T);
        double t = t0;
        while (t < T)
        {
            t += gsl_ran_exponential(rng.get(), 1./r0);
            if (t < T)
                st.push_back(Spike(t));
        }
        return st;
    }
    

    SpikeTrain generate_inhomogeneous_poisson(RNG& rng, const double T, 
            const Array<double>& r, const double rmax, const double t0)
    {
        SpikeTrain st(T);
        double t = t0;
        const double dt = T/r.size();
        while (t < T)
        {
            t += gsl_ran_exponential(rng.get(), 1./rmax);
            if ((t < T) && (gsl_rng_uniform(rng.get()) < r[t/dt]/rmax))
                st.push_back(Spike(t));
        }
        return st;
    }
};
};
