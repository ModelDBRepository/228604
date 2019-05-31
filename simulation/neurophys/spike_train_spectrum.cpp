#include <math.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include "spike_train_spectrum.h"

SpikeTrainSpectrumCalculator::SpikeTrainSpectrumCalculator(
        const R2CFourierTransform& ft, const size_t N,
        const double T, const double t0):
    csc_(ft, N, T), N_(N), T_(T), t0_(t0), n_trials_(0)
{
    timeseries_ = (double*)fftw_malloc(sizeof(double) * N);
    memset(timeseries_, 0, sizeof(double) * N);
    summed_spec_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N/2+1));
    memset(summed_spec_, 0, sizeof(fftw_complex) * (N/2+1));
}

SpikeTrainSpectrumCalculator::~SpikeTrainSpectrumCalculator()
{
    fftw_free(timeseries_);
    fftw_free(summed_spec_);
}

void SpikeTrainSpectrumCalculator::feed(const double t)
{
    if (t > t0_)
    {
        while (t > t0_ + T_)
        {
            csc_.calc(timeseries_, timeseries_, summed_spec_, 
                    CrossSpecCalculator::ADDITIVE);
            n_trials_++;
            memset(timeseries_, 0, sizeof(double) * N_);
            t0_ += T_;//ceil(t/T_)*T_;
        }
        const double dt = T_/N_;
        timeseries_[(int)((t-t0_)/T_ * N_)] += 1./dt;
    }
}

void SpikeTrainSpectrumCalculator::get_spec(fftw_complex* out_spec)
{
    const size_t N_tr = get_num_trials();
    for (unsigned int i = 0; i < N_/2+1; i++)
        out_spec[i] = summed_spec_[i]/N_tr;
}
