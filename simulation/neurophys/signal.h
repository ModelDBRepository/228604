#ifndef SIGNAL_H_
#define SIGNAL_H_

#include <gsl/gsl_statistics_double.h>
#include <complex>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include "fourier_transform.h"
#include "array.h"
#include "rng.h"
#include <assert.h>

// this is slowly turning into sth that should maybe rather be called
// "stochastic processes" rather than signal

namespace neurophys {

class SpectrumDescriptor
{
public:
    virtual double get(const double f) const = 0;
};

class WhiteSpectrum: public SpectrumDescriptor
{
public:
    WhiteSpectrum(const double s0): s0_(s0) {}
    virtual double get(const double f) const { return s0_; }
private:
    const double s0_;
};

class BandLimitedFlatSpectrum: public SpectrumDescriptor
{
public:
    BandLimitedFlatSpectrum(const double s0, const double f0, const double f1):
        s0_(s0), f0_(f0), f1_(f1) 
        {}
    virtual double get(const double f) const 
    { 
        return (f >= f0_ && f < f1_) ? s0_ : 0.; 
    } 
private:
    const double s0_;
    const double f0_;
    const double f1_;
};

class PowerLawSpectrum: public SpectrumDescriptor
{
public:
    PowerLawSpectrum(const double expo, const double f0, const double f1): 
        expo_(expo), f0_(f0), f1_(f1), A_(0)
    {
        // normalization, so that variance = 1
        A_ = 1./2 * (1.-expo) / (pow(f1,1.-expo) - expo*pow(f0,1.-expo));
    }
    virtual double get(const double f) const 
    {
        if (f < f0_) return A_ * pow(f0_, -expo_);
        else if (f < f1_) return A_ * pow(f, -expo_);
        else return 0.;
    }
private:
    const double expo_;
    const double f0_;
    const double f1_;
    double A_;
};

namespace signal {

Array<double> generate_gaussian(RNG& rng, 
        const C2RFourierTransform& c2rft, const SpectrumDescriptor& spec, 
        const size_t N, const double T);

};

}

#endif /* SIGNAL_H */
