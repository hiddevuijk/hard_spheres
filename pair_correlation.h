#ifndef GUARD_PAIR_CORRELATION_H
#define GUARD_PAIR_CORRELATION_H


#include <vector>
#include <math.h>
#include <fstream>
#include <string>

#include "xyz.h"
#include "particle.h"


class PairCorrelation {

public:
	// constructor
	PairCorrelation(unsigned int NN, double LL, unsigned int Nbinn, double bss);

	// sample
	void sample( const std::vector<Particle> &particles );

	// normalize
	void normalize();

	// get current number of samples taken
	int get_Nsamp() { return Nsamp;}

	// write to out stream
	void write(std::ostream & out);
	// write to file called outname
	void write(std::string outname);
private:

	unsigned int N;	// number of particles
	double L;		// box size. Use pbc if L>0,
					// if L<=0 don't use pbc.
	unsigned int Nbin;	// number of bins
	double bs;		// bin size

	std::vector<double> bins;	// mid values of the bins 
	std::vector<double> pc;		// pair correlation function
	
	unsigned int Nsamp;	// number of samples
	bool is_normalized;	// if true, pc is normalized;
						// sampling is no longer possible.
};

#endif
