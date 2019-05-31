#include <gsl/gsl_errno.h>
#include "autoparams.h"
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <time.h>
#include "integrator.h"

using namespace params;
using namespace integrator;

const double EPS = 1e-2; // relative accuracy of each integration
double D = 0;

bool warned_already = false;
void my_error_handler(const char* reason, const char* file, int line, int gsl_errno)
{
     // ignore roundoff errors
     if (gsl_errno == GSL_EROUND)
     {
        if (!warned_already)
        {
            warned_already = true;
            cerr << "GSL_EROUND ignored!" << endl;
        }
     }
     else 
     {
         gsl_set_error_handler(0);
         gsl_error(reason, file, line, gsl_errno);
     }
}

inline double U(const double x)
{
    return -mu*x-x*x*x/3;
}

double T1_i2(const double y, const void* p)
{
    const double x = *(static_cast<const double*>(p));
    return exp((U(x)-U(y))/D);
}


double T1_i1(const double x, const void* p)
{
    return integrate_infl(T1_i2, x, EPS, &x);
}

double T1()
{
    return 1./D * integrate_inf(T1_i1, EPS); 
}


double T1_BruLat_i1(const double z, const void* p)
{
    return exp(-mu*z*z-D*D/12*z*z*z*z*z*z);
}

double T1_BruLat()
{
    return pow(M_PI,0.5) * integrate_inf(T1_BruLat_i1, EPS);
}

double dT2_i3(const double z, const void* p)
{
    const double x = *(static_cast<const double*>(p));
    return exp((U(x)-U(z))/D);
}

double dT2_i2(const double y, const void* p)
{
    const double x = *(static_cast<const double*>(p));
    return exp((U(y)-U(x))/D);
}

double dT2_i1(const double x, const void* p)
{
    const double i3 = integrate_infl(dT2_i3, x, EPS, &x);
    return integrate_infu(dT2_i2, x, EPS, &x) * i3*i3;
}
double dT2()
{
    return 2./D/D * integrate_inf(dT2_i1, EPS);
}

int main(int argc, char** argv)
{

	if (params::acquire(argc, argv) != 0) return 1;

    const double infty = fabs(std::min(vr, -(params::integ_infty+1)));
    mu = mu + a_e * rin_e;
    D = a_e * a_e * rin_e;
    
    gsl_set_error_handler(my_error_handler);

    const double t1 = T1();
    const double dt2 = dT2();
    const double r = 1./t1;
    const double cv = sqrt(dt2)/t1;
    cout.precision(4);
    
    cout << "# first moment" << endl;
    cout << "T1 = " << t1 << endl;    
    cout << "# firing rate" << endl;
    cout << "r0 = " << r << endl;
    cout << "# CV" << endl;
    cout << "cv = " << cv << endl;
    return 0;
}
