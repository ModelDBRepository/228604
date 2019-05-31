#include "signal.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <string.h>
#include "array_funcs.h"

using namespace std;

namespace neurophys {

namespace signal {
    
Array<double> generate_gaussian(RNG& rng, 
        const C2RFourierTransform& c2rft, const SpectrumDescriptor& spec, 
        const size_t N, const double T)
{
    Array<complex<double> > stilde = arr::zeros<complex<double> >(N/2 + 1);
    stilde[0] = gsl_ran_gaussian(rng.get(), sqrt(2 * spec.get(0) * T / 2));
    for (int i = 1; i < N/2; i++) 
    {
        const double sigma = sqrt(spec.get(i*1./T) * T / 2);
        stilde[i] = complex<double>(gsl_ran_gaussian(rng.get(), sigma), 
                                    gsl_ran_gaussian(rng.get(), sigma));
    }
    stilde[N/2] = gsl_ran_gaussian(rng.get(), sqrt(2 * spec.get(N/2*1./T) * T / 2));
    return c2rft.transform(stilde);
}

}

}
