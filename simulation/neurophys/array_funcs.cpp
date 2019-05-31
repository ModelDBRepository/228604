#include "array_funcs.h"
#include <cmath>

namespace neurophys {
namespace arr {


double mean(const Array<double>& ar)
{
    return gsl_stats_mean(const_cast<double*>(ar.constdata()), 1, ar.size());
}
    
double var(const Array<double>& ar)
{
    return gsl_stats_variance(const_cast<double*>(ar.constdata()), 1, ar.size());
}

Array<double> range(const double a, const double b, const double inc)
{
    Array<double> result(std::max(0, static_cast<int>(std::ceil((b-a)/inc))));
    
    double d = a;
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = d;
        d += inc;
    }
    return result;
}

// if i add more of these, this could be done with less code bloat
Array<double> cos(const Array<double>& x)
{
    Array<double> result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::cos(x[i]);
    }
    return result;
}

Array<double> sin(const Array<double>& x)
{
    Array<double> result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::sin(x[i]);
    }
    return result;
}

Array<std::complex<double> > conj(const Array<std::complex<double> >& x)
{
    Array<std::complex<double> > result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::conj(x[i]);
    }
    return result;
}

Array<double> real(const Array<std::complex<double> >& x)
{
    Array<double> result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::real(x[i]);
    }
    return result;
}
Array<double> imag(const Array<std::complex<double> >& x)
{
    Array<double> result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::imag(x[i]);
    }
    return result;
}
Array<double> abs(const Array<std::complex<double> >& x)
{
    Array<double> result(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        result[i] = std::abs(x[i]);
    }
    return result;
}


}}
