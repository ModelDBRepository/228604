#ifndef FOURIER_TRANSFORM_H_
#define FOURIER_TRANSFORM_H_

#include <complex>
#include <fftw3.h>
#include "array.h"

namespace neurophys {

class R2CFourierTransform {
public:
    R2CFourierTransform(const size_t N, const double T);
    virtual ~R2CFourierTransform();
    Array<std::complex<double> > transform(const Array<double>& src) const;
private:
    const size_t N_;
	const double T_;
	fftw_plan fp_;
};

class C2RFourierTransform {
public:
    C2RFourierTransform(const size_t N, const double T);
    virtual ~C2RFourierTransform();
    Array<double> transform(const Array<std::complex<double> >& src) const;
private:
    const size_t N_;
	const double T_;
	fftw_plan fp_;
};

}
#endif // FOURIER_TRANSFORM_H_

