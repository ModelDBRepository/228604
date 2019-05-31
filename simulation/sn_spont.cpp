#include <iostream>
#include "autoparams.h"
#include <complex>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include "neurophys/array.h"
#include "neurophys/array_funcs.h"
#include "neurophys/averager.h"
#include "neurophys/signal.h"
#include "neurophys/fourier_transform.h"
#include "neurophys/synapse.h"
#include "neurophys/spike_train.h"
#include "neurophys/histogram_funcs.h"

using namespace std;
using namespace neurophys;
using namespace params;

int main(int argc, char** argv)
{
    if (params::acquire(argc, argv) != 0) return 1;
    RNG rng(params::seed);

    const double T = 1./df;
    const size_t N = round(2 * f_max * T);
	const size_t N_spec = N / 2 + 1;
    const double dt_coarse = T / N;

    Histogram<unsigned int> vhist(vhist_N, vhist_l, vhist_r);
    
    bool pif = (model == string("pif"));
    bool lif = (model == string("lif"));
    bool qif = (model == string("qif"));
    
    double v = vr;
    double dt_thr_hit = -1;
    double t = -T;
    double refr_until = -T;
    bool direct_hit = false;
    int n_spikes = 0;
    int n_direct_hits = 0;
    
    R2CFourierTransform ft(N, T);

    Array<double> vs = arr::zeros<double>(N);
    Array<complex<double> > sxx = arr::zeros<complex<double> >(N_spec);

    double t_next_spike = t + gsl_ran_exponential(rng.get(), 1./rin_e);
    double t_next_sample = t + gsl_ran_exponential(rng.get(), 1./r_sample);
    double t_last_spike = -T;

    SimpleAverager isi_avg;

    for (int ntr = -1; ntr < N_trials; ntr++)
    {
        vs = arr::zeros<double>(N);
        Array<double> x = arr::zeros<double>(N);
        
        while (t < (ntr + 1) * T)
        {
            double t_next_event = min(t_next_spike, t_next_sample);
            // time to next deterministic threshold hit?
            double t_thr_hit = max(t, refr_until);
            if (pif)
            {
                if (mu > 0)
                    t_thr_hit += (vt - v) / mu;
                else t_thr_hit = t_next_event + 1;
            }
            if (lif)
            {
                if (vt < mu)
                    t_thr_hit += log((v - mu)/(vt - mu));
                else t_thr_hit = t_next_event + 1;
            }
            if (qif)
            {
                if (mu > 0)
                    t_thr_hit += 1. / sqrt(mu) *
                        (atan(vt / sqrt(mu)) - 
                         atan(v / sqrt(mu)));
                else if (v > sqrt(-mu))
                    t_thr_hit += -1. / sqrt(-mu) * 
                        (atanh(sqrt(-mu) / vt) -
                         atanh(sqrt(-mu) / v));
                else t_thr_hit = t_next_event + 1;
            }

            if (t_next_event < t_thr_hit)
            {
                // subtract whatever remains of the refractory period
                double dt = max(0., t_next_event - max(refr_until, t)); 
                t = t_next_event;
                // evolve deterministically
                if (pif)
                {
                    v += dt * mu;
                }
                else if (lif)
                {
                    const double expfac = exp(-dt);
                    v = v * expfac + mu * (1. - expfac);
                }
                else if (qif)
                {
                    if (mu > 0)
                        v = sqrt(mu) * tan(sqrt(mu) * dt + atan(v / sqrt(mu)));
                    else if ((v < sqrt(-mu)) && (v > -sqrt(-mu)))
                        v =  sqrt(-mu) * tanh(- sqrt(-mu) * dt 
                                              + atanh(v/sqrt(-mu)));
                    else
                        v = -sqrt(-mu) * 1./tanh(sqrt(-mu) * dt 
                                                 - atanh(sqrt(-mu)/ v));
                }
                // bla
                // next event is an incoming spike?
                if (t_next_spike < t_next_sample) 
                {
                    if (t > refr_until)
                    {
                        // add spike
                        v += (expweights ? 
                                std::min(gsl_ran_exponential(rng.get(), a_e), a_e_cutoff) : a_e);
                        // did this spike kick us accross threshold?
                        if (v > vt)
                        {
                            direct_hit = true;
                            if (ntr > -1)
                            {
                                n_direct_hits++;
                            }
                        }
                    }
                    t_next_spike = t + gsl_ran_exponential(rng.get(), 1./rin_e);
                }
                else // next event is sampling the voltage
                {
                    if (ntr > -1)
                    {
                        if (t > refr_until)
                        {
                            vhist.feed(v);
                        }
                        else
                        {
                            vhist.feed(vhist_r + 1); // let it go to the overlow bin for correct normalization
                        }
                        if (static_cast<int>((t-ntr*T)/dt_coarse) < N)
                        {
                            vs[static_cast<int>((t-ntr*T)/dt_coarse)] = v;
                        }
                    } 
                    t_next_sample = 
                        t + gsl_ran_exponential(rng.get(), 1./r_sample);
                }
            }
            if ((t_thr_hit < t_next_event) || direct_hit)
            {
                if (direct_hit)
                {
                    direct_hit = false;
                }
                else
                {
                    t = t_thr_hit;
                }
                if (ntr > -1)
                {
                    n_spikes++;
                    isi_avg.feed(t-t_last_spike);
                    const int ti = static_cast<int>((t-ntr*T)/dt_coarse);
                    if (ti < N)
                    {
                        x[ti] += 1./dt_coarse;
                    }
                }
                v = vr; // reset
                t_last_spike = t;
                refr_until = t + tr; // do nothing while refractory 
            }
        }

        if (ntr > -1)
        {
            Array<complex<double> > xtilde = ft.transform(x);
            sxx += 1./N_trials * 1./T * xtilde * arr::conj(xtilde);
        }
    }

    vhist.write_ascii("vhist", Histogram<unsigned int>::NORMALIZED);
    arr::columns_to_file("vs", arr::range(0, T, dt_coarse), vs);
    arr::columns_to_file("sxx", arr::range(0, f_max, 1./T), arr::real(sxx));

    cout << "n_spikes = " << n_spikes << endl;
    cout << "al = " << static_cast<double>(n_direct_hits)/n_spikes << endl;
    cout << "T = " << isi_avg.get_mean() << endl;
    cout << "T2 = " << isi_avg.get_var() << endl;
    cout << "r0 = " << 1./isi_avg.get_mean() << endl;
    cout << "cv = " << sqrt(isi_avg.get_var())/isi_avg.get_mean() << endl;

    return 0;
}
