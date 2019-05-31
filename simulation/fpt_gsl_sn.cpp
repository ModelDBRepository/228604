/*
 * C++ implementation for evaluating the integrals involved in the first and 
 * second moment of the FPT density. uses the gnu scientific library for i
 * integration
 *
 * XXX this cannot yet deal with a refractory period!
 *
 * the integrals have been transformed into a more suitable form on paper,
 * as I remember this was kind of messy, and I should definitely tex it
 *
 */

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
const double deltav = 0; // how close to problematic points do we allow 
                            // evaluation?
double vmin = 0;
double al = 1.;
double s = 0;

// we restrict ourselfs to neuron models with a most one stable and one
// unstable fixed point in the "-" dynamics, with the stable on left of the 
// unstable one, both of which makes sense in any case i can think of
class NeuronModel
{
public:
    virtual double f(const double v) = 0; // the nonlinearity
    virtual double dfdv(const double v) = 0; // its derivative
    virtual double phi(const double v) = 0; 
        // phi(v)=\int^v dx kp/(f(x)+s) + rin_e/(f(x)-s)
    virtual bool has_stable_fp() = 0; // stable fixed point in "-" dyn?
    virtual bool has_unstable_fp() = 0; // unstable FP in "-" dyn?
    virtual double get_stable_fp() = 0; // return stable FP
    virtual double get_unstable_fp() = 0; // return unstable FP
};

class PIF: public NeuronModel
{
public:
    virtual double f(const double v) { return mu; }
    virtual double dfdv(const double v) { return 0; }
    virtual double phi(const double v) { return v*(1./a_e+rin_e/mu); } 
    virtual bool has_stable_fp() { return false; }
    virtual bool has_unstable_fp() { return false; }
    virtual double get_stable_fp() { return 0.; }
    virtual double get_unstable_fp() { return 0.; }
};

class LIF: public NeuronModel
{
public:
    virtual double f(const double v) { return mu-v; }
    virtual double dfdv(const double v) { return -1.; }
    virtual double phi(const double v) 
    { 
        return v/a_e-rin_e*log(fabs(mu-v));
    }
    virtual bool has_stable_fp() { return mu<vt; }
    virtual bool has_unstable_fp() { return false; }
    virtual double get_stable_fp() { return mu; }
    virtual double get_unstable_fp() { return 0.; }
};

class QIF: public NeuronModel
{
public:
    virtual double f(const double v) { return mu+v*v; }
    virtual double dfdv(const double v) { return .5*v; }
    virtual double phi(const double v) 
    { 
        return v/a_e+rin_e*
            (mu > 0 ? 
                1./sqrt(mu)*atan(v/sqrt(mu)) :
                -1./sqrt(-mu) * 
                    ((sqrt(-mu) > v && v > -sqrt(-mu)) ? 
                        atanh(v/sqrt(-mu)) : 
                        .5*log((v/sqrt(-mu)+1.)/(v/sqrt(-mu)-1.)) // == acoth(v/sqrt(-mu))
            ));
    }
    virtual bool has_stable_fp() { return mu<0; }
    virtual bool has_unstable_fp() { return (mu<0 && sqrt(-mu) < vt); }
    virtual double get_stable_fp() { return -sqrt(-mu); }
    virtual double get_unstable_fp() { return sqrt(-mu); }
};

NeuronModel* neu = 0;

// heaviside step function
inline double heav(double x)
{
    return (x < 0 ? 0 : 1);
}

inline double phi(const double v)
{
    if (neu != 0) 
        return neu->phi(v);
    return 0.;
}

inline double f(const double v)
{
    if (neu != 0)
        return neu->f(v);
    return 0.;
}


// structure to pass some data along with function pointers that called by the
// integration routines
struct NVStruct {
    unsigned int n;
    double v;
    unsigned int i;
};


struct Interval
{
    double l;
    double r;
    double c;
};
unsigned int n_ints = 0;
Interval intervals[3];

double H0i(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const double v = nv->v;
    return (exp(phi(x)-phi(v)) * heav(x-vr))/a_e;
}


double H0(unsigned int i, const double v)
{
    NVStruct nv = {0,v,i};
    const double ci = intervals[i].c;

    return 1/f(v)* 
            ((heav(vr-ci)*heav(v-vr) - heav(ci-vr)*heav(vr-v)) * exp(phi(vr)-phi(v)) * (-1.)  
            - integrate(H0i, ci, v, EPS, &nv)); 
}

double J1i(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    return -H0(nv->i,x);
}

double J1(const double v)
{
    //Interval intervals[3];
    //const unsigned int n_ints = construct_intervals(intervals);
    // what interval is v in? 
    unsigned int iv = 0;
    for (iv = 0; iv < n_ints; iv++)
    {
        if (intervals[iv].r >= v)
            break;
    }
    double res = 0;
    unsigned int i = 0;
    for (i = 0; i < iv + 1; i++)
    {
        NVStruct nv = {1,v,i};
        res += integrate(J1i, intervals[i].l, std::min(intervals[i].r, v), EPS, &nv);
    }
    //cout << "in J1: " << i << " integrals evaluated\n";
    return res;
}

double J1alt_i3(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    return exp(phi(nv->v)-phi(x))/f(x);
}

double J1alt_i2(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const unsigned int i = nv->i;
    const double cibar = intervals[i].l + intervals[i].r - intervals[i].c;
    const NVStruct nvinner = {1,x,i};
    return heav(x-vr)/a_e*integrate(J1alt_i3, x, cibar, EPS, &nvinner);
}

double J1alt_i1(const double x, const void* p)
{
    return exp(phi(vr)-phi(x))/f(x);
}

double J1alt_vt()
{
    const double c0bar = intervals[0].l + intervals[0].r - intervals[0].c;
    double res = integrate(J1alt_i1, vr, c0bar, EPS);
    for (unsigned int i = 0; i < n_ints; i++)
    {
        NVStruct nv = {1,0,i};
        res +=  integrate(J1alt_i2, intervals[i].l, intervals[i].r, EPS, &nv);
    }
    return res;

}

unsigned int construct_intervals(Interval* intervals)
{
    unsigned int n_ints = 0;

    if (f(vr) > 0) 
    {
        intervals[0].l = vr;
        intervals[0].c = vr;
        if (neu->has_stable_fp()) {
            intervals[1].l = intervals[0].r = neu->get_stable_fp();
            if (neu->has_unstable_fp()) {
                n_ints = 3;
                intervals[2].l = intervals[1].r = intervals[2].c 
                    = intervals[1].c = neu->get_unstable_fp();
                intervals[2].r = vt;
            } else {
                n_ints = 2;
                intervals[1].r = intervals[1].c = vt;
            }
        } else {
            n_ints = 1;
            intervals[0].r = vt;
        }
    }
    else
    {
        intervals[0].l = neu->get_stable_fp();
        if (neu->has_unstable_fp()) {
            n_ints = 2;
            intervals[1].l = intervals[0].r = intervals[1].c 
                = intervals[0].c = neu->get_unstable_fp();
            intervals[1].r = vt;
        } else {
            n_ints = 1;
            intervals[0].r = intervals[0].c = vt;
        }
    }
    for (int i = 0; i < n_ints; i++)
        cout << "# " << intervals[i].l << " " << intervals[i].r << " " << intervals[i].c << endl;
    return n_ints;
}

inline double cbar(const Interval& inter)
{
    return (inter.l == inter.c ? inter.r : inter.l);
}

double J2i_inner(const double x, const void* p)
{

    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const double v = nv->v;
    return exp(phi(v)-phi(x))/f(x);
}

double a_i(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const unsigned int i  = nv->i;
    NVStruct nvinner = {2,x,i}; 
    //const double cbar = (intervals[i].l + intervals[i].r - intervals[i].c);
    return 1./a_e*integrate(J2i_inner, x, cbar(intervals[i]), EPS, &nvinner);
}

double J2_vt_faster_i(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const unsigned int i  = nv->i;
    NVStruct nvi = {2,x,i}; 
    
    double a_intsum = integrate(a_i, x, intervals[i].r, EPS, &nvi);
    for (unsigned int j = i+1; j < n_ints; j++)
    {
        NVStruct nvj = {2,vt,j};
        a_intsum += integrate(a_i, intervals[j].l, intervals[j].r, EPS, &nvj);
    }
    //const double cbar = (intervals[i].l + intervals[i].r - intervals[i].c);
    const double innerint = integrate(J2i_inner, x, cbar(intervals[i]), EPS, &nvi);
    return -H0(i,x)*(a_intsum + innerint);
}


double J2_vt_faster()
{
    double res = 0;
    unsigned int i = 0;
    for (i = 0; i < n_ints; i++)
    {
        NVStruct nv = {2,vt,i};
        res += integrate(J2_vt_faster_i, intervals[i].l, intervals[i].r, EPS, &nv);
    }
    return 2*res;
}


double J2alt_i1(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const double v = nv->v;
    const unsigned int i = nv->i;
    NVStruct nvi = {2,x,i}; 
    double a_intsum = integrate(a_i, x, intervals[i].r, EPS, &nvi);
    for (unsigned int j = i+1; j < n_ints; j++)
    {
        NVStruct nvj = {2,vt,j};
        a_intsum += integrate(a_i, intervals[j].l, intervals[j].r, EPS, &nvj);
    }
    const double cbar = (intervals[i].l + intervals[i].r - intervals[i].c);
    const double innerint = integrate(J2i_inner, x, cbar, EPS, &nvi);
    return exp(phi(v)-phi(x))/f(x) * (a_intsum + innerint);
}

double J2alt_i2(const double x, const void* p)
{
    const NVStruct* nv = static_cast<const NVStruct*>(p);
    const unsigned int i = nv->i;
    const double cbar = (intervals[i].l + intervals[i].r - intervals[i].c);
    NVStruct nvi = {2,x,i}; 
    return heav(x-vr)/a_e*integrate(J2alt_i1,x,cbar,EPS,&nvi);
}


double J2alt_vt()
{
    const double c0bar = intervals[0].l + intervals[0].r - intervals[0].c;
    NVStruct nv1 = {2,vr,0};
    double res = integrate(J2alt_i1, vr, c0bar, EPS, &nv1);
    for (unsigned int i = 0; i < n_ints; i++)
    {
        NVStruct nv = {2,vt,i};
        res += integrate(J2alt_i2, intervals[i].l, intervals[i].r, EPS, &nv);
    }
    return 2*res;
}


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

int main(int argc, char** argv)
{

	if (params::acquire(argc, argv) != 0) return 1;

    const double infty = fabs(std::min(vr, -(params::integ_infty+1)));
    //integrator::infty = params::integ_infty;
    
    if (model == string("pif"))
        neu = new PIF();
    else if (model == string("lif"))
        neu = new LIF();
    else if (model == string("qif"))
        neu = new QIF();
    else throw invalid_argument("unknown neuron model");

    if (tr > 1e-12) throw invalid_argument("refr. period tr not implemented yet!");

    n_ints = construct_intervals(intervals);

    vmin = f(vr)>0 ? vr : 
        (neu->has_stable_fp() ? neu->get_stable_fp() : -infty);
    if (vmin < -infty) vmin = -infty;
    
    gsl_set_error_handler(my_error_handler);

    // calculate alpha if needed
    /*
    if (f(vt)>0)
    {
        if (neu->has_unstable_fp())
            al = 1.+1./a_e*exp(-phi(vt))*integrate(
                [](double x, const void* p){ return exp(phi(x)); }
                , vt, neu->get_unstable_fp(), EPS);
        else
            al = 1.-exp(-phi(vt))*(exp(phi(vr))+1./a_e*integrate(
                [](double x, const void* p){ return exp(phi(x)); }
                , vr, vt, EPS));
    }*/
    const double T1 = J1(vt);
    //const double T1 = J1alt_vt();
    //const double T2 = J2(vt);
    //const double T2=0;
    //const double T2 = J2alt_vt();
    const double T2 = J2_vt_faster();
    const double r = 1./T1;
    const double cv = sqrt(T2-T1*T1)/T1;
/*
    cout << "# mu\trin_e\ta_e\tvr\tvt\tD\t\t\tT1\tT2\tr\tcv\tal\n";
    cout.precision(4);
    cout << mu << "\t" << rin_e << "\t" << a_e << "\t" << vr << "\t"
        << vt << "\t" << rin_e*a_e*a_e << "\t\t\t" 
        << T1 << "\t" << T2  << "\t"<< r << "\t" << cv << "\t" << al << endl;
*/
    cout << "# first moment" << endl;
    cout << "T1 = " << T1 << endl;    
    cout << "# second moment" << endl;
    cout << "T2 = " << T2 << endl;    
    cout << "# firing rate" << endl;
    cout << "r0 = " << r << endl;
    cout << "# CV" << endl;
    cout << "cv = " << cv << endl;


    delete neu;
    return 0;
}
