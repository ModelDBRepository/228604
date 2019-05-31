/* 
 * File:   rng.h
 * Author: fedro
 *
 * Created on November 3, 2014, 5:33 PM
 */

#ifndef RNG_H
#define	RNG_H

#include <gsl/gsl_rng.h>

class RNG
{
public:
    RNG(const gsl_rng_type* rng_type, const int seed=0)
    {
        rng_ = gsl_rng_alloc(rng_type);
        gsl_rng_set(rng_, seed);
    }
    RNG(const int seed=0)
    {
        rng_ = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(rng_, seed);
    }
    virtual ~RNG()
    {
        gsl_rng_free(rng_);
    }
    gsl_rng* get() { return rng_; }
private:
    gsl_rng* rng_;
};

#endif	/* RNG_H */

