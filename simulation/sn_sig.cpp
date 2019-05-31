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
#include "integrator.h"

using namespace std;
using namespace neurophys;
using namespace params;

int main(int argc, char** argv)
{
    if (params::acquire(argc, argv) != 0) return 0;
    RNG rng(params::seed);

    const double T = 1./df;
    const size_t N = round(2 * f_max * T);
	const size_t N_spec = N / 2 + 1;
    const double dt_coarse = T / N;
    const size_t N_sim = round(T / dt);

    R2CFourierTransform ft(N, T);
    C2RFourierTransform ift(N_sim, T);
    
    Histogram<unsigned int> vhist(vhist_N, vhist_l, vhist_r);
    
    Array<complex<double> > S_ii = arr::zeros<complex<double> >(N_spec);
    Array<complex<double> > S_xx = arr::zeros<complex<double> >(N_spec);
    Array<complex<double> > S_sx = arr::zeros<complex<double> >(N_spec);
    Array<complex<double> > S_ss = arr::zeros<complex<double> >(N_spec);
    
    bool lif = (model == string("lif"));
    bool eif = (model == string("eif"));

    double v = vr;
    double fv = 0;
    Array<double> vs(N_sim);
    Array<double> xin;
    SpikeTrain st;
    Array<double> r;
    int i = N_sim;
    int n_spikes = 0;
    Array<double> s;
    if (f_sig > 0)
    {
        f_sig = f_sig - fmod(f_sig, df);
        std::cout << "f_sig_c = " << f_sig << std::endl;
        // this needs to be at least one period too long, so let's take two, just to be sure:
        s = arr::cos(2 * M_PI * f_sig * (arr::range(0, T+2./f_sig, dt))); 
    }
    int i_phase = 0;

    for (int ntr = -1; ntr < N_trials; ntr++)
    {
        Array<double> x = arr::zeros<double>(N_sim);
        vs = arr::ones<double>(N_sim) * vr;
        if (f_sig < 0)
        {
            s = signal::generate_gaussian(rng, ift, 
                BandLimitedFlatSpectrum(1./(2*f_c), 0., f_c), N_sim, T);
        }
        else
        {
            i_phase = gsl_rng_uniform(rng.get()) * 1. / f_sig / dt;
        }
        r = rin_e * (1 + eps_r * s);
        const double rmax = arr::max(r);
        
        st = spike_train::generate_inhomogeneous_poisson(rng, T, r, rmax);
        for (SpikeTrain::iterator it = st.begin(); it != st.end(); ++it)
        {
            // assign each input spike an exp. distributed weight
            it->amplitude = gsl_ran_exponential(rng.get(), a_e);
        }
        xin = spike_train::to_array(st, N_sim, T);
        
        for (i = i - N_sim; i < N_sim; i++)
        {
            if (lif)
                fv = mu - v;
            else if (eif)
                fv = mu - v + d * exp((v-vtb)/d);
            v += dt * ((fv + eps_v * s[i_phase + i]) + xin[i]);
            vs[i] = v;
            if (v > vt)
            {
                v = vr;
                x[i] = 1./dt;
                i += static_cast<int>(tr/dt);
                n_spikes++;
            }
        }
        
        if (ntr > -1)
        {
            Array<double> ssh(N_sim);
            memcpy(ssh.data(), s.data() + i_phase, sizeof(double) * N_sim);
            const Array<complex<double> > xintilde = 
                ft.transform(arr::downsample(xin, N));
            const Array<complex<double> > xtilde = 
                ft.transform(arr::downsample(x, N));
            const Array<complex<double> > stilde = 
                ft.transform(arr::downsample(ssh, N));
            S_ii += 1./N_trials * 1./T * (xintilde * arr::conj(xintilde));
            S_xx += 1./N_trials * 1./T * (xtilde * arr::conj(xtilde));
            S_sx += 1./N_trials * 1./T * (stilde * arr::conj(xtilde));
            S_ss += 1./N_trials * 1./T * (stilde * arr::conj(stilde));
            hist::feed_array_to_hist(vs, &vhist);
        }
    }
    
    // write spike times as a " "-separated list
    ofstream of;
    of.open(".spikes", ios::out);
    for (SpikeTrain::const_iterator it = st.begin(); it != st.end(); it++)
    {
        of << it->time << " ";
    }
    of << "\n";
    of.close();
    
    cout << "r0 = " << n_spikes/(T * N_trials) << endl;

    const Array<double> fs = arr::range(0, N_spec*df, df);
    const Array<double> ts = arr::range(0, T, dt);
    arr::columns_to_file("ii.spec", fs, arr::real(S_ii));      
    arr::columns_to_file("xx.spec", fs, arr::real(S_xx));      
    arr::columns_to_file("ss.spec", fs, arr::real(S_ss));      
    arr::columns_to_file("sx.spec", fs, arr::real(S_sx), arr::imag(S_sx));
    arr::columns_to_file("vs", ts, vs);  
    arr::columns_to_file("xin", ts, xin);  
    arr::columns_to_file("poprate", arr::range(0, T, dt_coarse), r);    
    
    if (f_sig > 0)
    {
        const double e = (eps_v > 0 ? eps_v : eps_r);
        cout << "resus = " << 
            std::real(S_sx[static_cast<int>(f_sig/df+0.5)] / (e/4/df)) << endl;
        cout << "imsus = " << 
            std::imag(S_sx[static_cast<int>(f_sig/df+0.5)] / (e/4/df)) << endl;
    }
    
    vhist.write_ascii("vhist", Histogram<unsigned int>::NORMALIZED);
    
    return 0;
}
