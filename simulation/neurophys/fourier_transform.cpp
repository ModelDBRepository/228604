#include "fourier_transform.h"

using namespace std;

namespace neurophys {

R2CFourierTransform::R2CFourierTransform(const size_t N, const double T):
    N_(N), T_(T)
{
    Array<double> tmpsrc(N);
    Array<complex<double> > tmpdst(N / 2 + 1);
    fp_ = fftw_plan_dft_r2c_1d(N, tmpsrc.data(), 
            reinterpret_cast<fftw_complex*>(tmpdst.data()), FFTW_MEASURE);
}

Array<complex<double> > R2CFourierTransform::transform(const Array<double>& src) const
{
    if (N_ != src.size())
        throw ArraySizeMismatchException();
    Array<complex<double> > dest(N_ / 2 + 1);
    fftw_execute_dft_r2c(fp_, const_cast<double*>(src.constdata()), 
            reinterpret_cast<fftw_complex*>(dest.data()));
    dest *=  T_/N_;
    return dest;
}

R2CFourierTransform::~R2CFourierTransform()
{
    fftw_destroy_plan(fp_);
}

C2RFourierTransform::C2RFourierTransform(const size_t N, const double T):
    N_(N), T_(T)
{
    Array<complex<double> > tmpsrc(N / 2 + 1);
    Array<double> tmpdst(N);
    fp_ = fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(tmpsrc.data()),
            tmpdst.data(), FFTW_MEASURE);
}

Array<double> C2RFourierTransform::transform(const Array<complex<double> >& src) const
{
    if (N_ / 2 + 1 != src.size())
        throw ArraySizeMismatchException();
    Array<double> dest(N_);
    fftw_execute_dft_c2r(fp_, 
            reinterpret_cast<fftw_complex*>(const_cast<complex<double>*>(src.constdata())), 
            dest.data());
    dest *=  1./T_;
    return dest;
}

C2RFourierTransform::~C2RFourierTransform()
{
    fftw_destroy_plan(fp_);
}

}
