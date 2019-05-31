#ifndef SYNAPSE_H_
#define SYNAPSE_H_

#include "spike_train.h"
#include <gsl/gsl_rng.h>

namespace neurophys {

class DynamicSynapse {
public:
	virtual SpikeTrain filter(const SpikeTrain& in) = 0;
};

class FDSynapse: public DynamicSynapse {
public:
	static const int STATIC = 0;
	static const int DEPRESSING = 1;
	static const int FACILITATING = 2;
	static const int GENERAL = DEPRESSING | FACILITATING;
	static const int KEEP_STATE_VARIABLES = 4; // no reset between trials
	
	FDSynapse(double F0, double tau_f, double tau_d, double delta, int flags): 
		F0_(F0), tau_f_(tau_f), tau_d_(tau_d), delta_(delta), flags_(flags), 
        FC0_(0), t_(0) {};

protected:
	const double F0_;
	const double tau_f_;
	const double tau_d_;
	const double delta_;
	const int flags_;
    double FC0_;
    double t_;
};

class DeterministicFDSynapse: public FDSynapse {
public:
	DeterministicFDSynapse(double F0, double tau_f, double tau_d, double delta, 
           int flags):
		FDSynapse(F0, tau_f, tau_d, delta, flags), D0_(0) {};
	
	virtual SpikeTrain filter(const SpikeTrain& in);
private:
    double D0_;
};

// rng in constructor is not so nice -- make explicit as param of filter and 
// screw poplymorphism?
class StochasticFDSynapse: public FDSynapse {
public:
	StochasticFDSynapse(double F0, double tau_f, double tau_d, double delta, 
            int flags, gsl_rng* rng, unsigned int n_release_sites):
		FDSynapse(F0, tau_f, tau_d, delta, flags), rng_(rng), n_release_sites_(n_release_sites), n_ves_present_(n_release_sites) {};

	virtual SpikeTrain filter(const SpikeTrain& in);
private:
	gsl_rng* rng_;
	const unsigned int n_release_sites_;
	unsigned int n_ves_present_;	
};

}

#endif // SYNAPSE_H_
