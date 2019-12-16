#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>

#include "xyz.h"
#include "particle.h"
#include "potential.h"

#include <vector>
#include <cmath>
#include <boost/random.hpp>

class System {
public:
    System( int seed, unsigned int N, double L, Potential potential, double d, double rv);

	// random number generator
    // uniform distribution [-1,1]
	const boost::uniform_real<double> udist11;
    // uniform distribution [0,1]
	const boost::uniform_real<double> udist01;
    // random int [0,N]  
    const boost::random::uniform_int_distribution<int> u_int_dist;
    
	int seed;
	boost::mt19937 rng;		
	boost::variate_generator<boost::mt19937&,
		boost::uniform_real<double> > rudist11;
	boost::variate_generator<boost::mt19937&,
		boost::uniform_real<double> > rudist01;

    // initialize on a lattice
    void init_random();
    void init_random_shift();

    // attemt to move a particle
    void mc_move();
    void mc_move_verlet();

     
    void update_verlet_list();

    unsigned int N;      // number of particles
    double L;   // system size
    Potential potential; 
    double d;  // move size 
    double rv;  // radius of Verlet sphere

    std::vector<Particle> particles;
    std::vector<Particle> particles_before_update;;
    std::vector<std::vector<int> > verletList;
    std::vector<unsigned int> neighbour_number_list;


    double max_diff; // maximum distance moved before update
                     // of the Verlet list

    long int Ntry; // number of attempted moves
    long int Nacc;  // number of accepted moves
    long int Nverlet; // number of verlet list updates

};

    


#endif
