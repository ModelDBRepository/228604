#ifndef SPIKE_TRAIN_SPECTRUM_H__
#define SPIKE_TRAIN_SPECTRUM_H__

#include <complex.h>
#include <fftw3.h>
#include "fourier_transform.h"
#include "spectra.h"

// a class that is continously fed spike times and whenever a time window T has
// passed uses the spikes collected in that window to calculate the power 
// spectrum
// N is the size of the timeseries, the size of the spectrum will be N/2+1!

class SpikeTrainSpectrumCalculator
{
public:
    SpikeTrainSpectrumCalculator(const R2CFourierTransform& ft, const size_t N,
        const double T, const double t0);
    virtual ~SpikeTrainSpectrumCalculator();
    void feed(const double t);
    size_t get_num_trials() { return n_trials_; };
    void get_spec(fftw_complex* out_spec); 
private:
    CrossSpecCalculator csc_;
    const size_t N_;
    const double T_;
    double t0_;
    double* timeseries_;
    fftw_complex* summed_spec_;
    size_t n_trials_;
};

#endif
