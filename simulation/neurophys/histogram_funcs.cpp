#include "histogram_funcs.h"

namespace neurophys {

namespace hist {    
    
void feed_array_increments_to_hist(const Array<double>& x, 
        Histogram<unsigned int>* hist, const double dt, const double tinc)
{
    const int inclen = tinc/dt;
    const int nincs = x.size()/inclen;
    for (int i = 0; i < nincs; i++)
    {
        double v = 0;
        for (int j = 0; j < inclen; j++)
        {
            v += x[i * inclen + j];
        }
        hist->feed(v*dt);
    }
}

void feed_array_to_hist(const Array<double>& x, Histogram<unsigned int>* hist)
{
    for (int i = 0; i < x.size(); i++)
    {
        hist->feed(x[i]);
    }
}

}
}
