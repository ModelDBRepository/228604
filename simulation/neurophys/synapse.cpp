#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "synapse.h"
#include "spike_train.h"
#include <iostream>
#include <math.h>

namespace neurophys {

SpikeTrain DeterministicFDSynapse::filter(const SpikeTrain& in)
{
    SpikeTrain out;
    
    double t = t_;
    double FC = FC0_;
	double D = D0_;
	double F = F0_ + (1. / (1. / (1. - F0_) + 1. / FC));
	double A = F * D;

    for (SpikeTrain::const_iterator it = in.begin(); it != in.end(); ++it)
    {
        double tn = it->time;
        
        if (flags_ & FACILITATING) 
		{
			FC = FC * exp(-(tn - t) / tau_f_);
			F = F0_ + (1. / (1. / (1. - F0_) + 1. / FC));
		} else F = F0_;
		if (flags_ & DEPRESSING)
			D = 1 - (1 - D) * exp(-(tn - t) / tau_d_);

		else D = 1.;
		
		A = F * D;
        
        out.push_back(Spike(it->time, it->amplitude * A));
        
        if (flags_ & FACILITATING)
				FC += delta_ * it->amplitude;
		if (flags_ & DEPRESSING)
				D -= A * it->amplitude;
        
        t = tn;
    }
    
    if (flags_ & KEEP_STATE_VARIABLES)
	{
		D0_ = D;
		FC0_ = FC;
	}	
    
    t_ = t - in.getT();
    return out;
}

SpikeTrain StochasticFDSynapse::filter(const SpikeTrain& in)
{
    SpikeTrain out;
    
    double t = t_;
    double FC = FC0_;
	double F = F0_ + (1. / (1. / (1. - F0_) + 1. / FC));
  
	unsigned int n_ves_present = n_ves_present_;

    for (SpikeTrain::const_iterator it = in.begin(); it != in.end(); ++it)
    {
        double tn = it->time;
    		
		// how many vesicles have we picked up in the time interval tn-t?
		unsigned int n_takeup = 0;
		double p_u = 1;	
		if (flags_ & DEPRESSING)
		{
			p_u = 1.-exp(-(tn-t)/tau_d_);
			n_takeup = gsl_ran_binomial(rng_, p_u, 
                    n_release_sites_ - n_ves_present);
		}
		else n_takeup = n_release_sites_ - n_ves_present;	
		n_ves_present += n_takeup;
		// how has the probability for a docked vesicle to be released upon 
        // spike arrival changed through facilitation?	
		if (flags_ & FACILITATING)
		{
			FC = FC * exp(-(tn - t) / tau_f_);
			F = F0_ + (1. / (1. / (1. - F0_) + 1. / FC));
		} else F = F0_;
		
		// now lets deal with the incoming spike
		
		// deal with facilitation wie gehabt
		if (flags_ & FACILITATING)
			FC += delta_ * it->amplitude;
        // releeease your vesicles
		unsigned int n_release = gsl_ran_binomial(rng_, F, n_ves_present);
		n_ves_present -= n_release;
		
        out.push_back(Spike(tn, (double)n_release/n_release_sites_));
	
		t =  tn;
	}	
	if (flags_ & KEEP_STATE_VARIABLES)
	{
		n_ves_present_ = n_ves_present;
		FC0_ = FC;
	}	        
    t_ = t - in.getT();
	return out;
}

}