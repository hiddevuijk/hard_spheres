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
    System(Potential potential, int N, double L, double d, double rv, int seed);

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

    // initialize on a fcc lattice
    void init_random();
   
    // attemt to move a particle
    void mc_move();

     
    void update_verlet_list();

    Potential potential; 
    std::vector<Particle> particles;
    std::vector<Particle> particles_before_update;;
    std::vector<std::vector<int> > verletList;

    int N;      // number of particles
    double L;   // system size
    double rv;  // radius of Verlet sphere
    double d;  // move size 

    long int Ntry; // number of attempted moves
    long int Nacc;  // number of accepted moves
    long int Nverl; // number of verlet list updates

};

    

void System::mc_move()
{
    // random displacement
    XYZ dr( d*rudist11(), d*rudist11(), d*rudist11() );

    // random particle index 
    int i = u_int_dist(rng); 

    Ntry += 1;

    double deltaU = 0;
    bool overlap = false;
    // ni loops over the Verlet list of particle i
    // and calculate difference in potential energy 
    //for(int ni=0; ni<N; ++ni) {
    //    if(ni == i) continue;
    //    overlap = potential.get_overlap(particles[i].r+dr, particles[ni].r,L);
    //    if( overlap) break;
    //}
    //if(overlap == false) {
    //    particles[i].r += dr;
    //    particles[i].r.pbc(L);
    //    Nacc += 1;
    //}
    for(int ni = 0; ni < verletList[i].size(); ++ni) {
        // combine ??
        overlap = potential.get_overlap(particles[i].r + dr, particles[ verletList[i][ni] ].r, L );
        if( overlap ) break;

        deltaU += potential.U( particles[i].r + dr, particles[ verletList[i][ni] ].r, L)
                  - potential.U(particles[i].r, particles[ verletList[i][ni] ].r , L) ;
     
    }
    if( !overlap ) {
        if( rudist01() < std::exp( -potential.beta*deltaU ) ) {
            particles[i].r += dr;

            // apply periodic boundary conditions
            //particles[i].r.pbc(L);

            Nacc += 1;

            // check if the Verlet list needs to be updated
            //update_verlet_list();
            if( xyz::dist_pbc(particles[i].r, particles_before_update[i].r,L) > 0.1*(rv-potential.rco) ) update_verlet_list();
        } 
    }
}

System::System(Potential potential, int N, double L, double d, double rv, int seed)
: potential(potential), N(N), L(L), d(d), rv(rv),
   udist01(0,1),udist11(-1,1), u_int_dist(0,N-1), rng(seed), rudist01(rng,udist01),
    rudist11(rng,udist11), particles(N), particles_before_update(N)
{
    Ntry = 0;
    Nacc = 0;
    Nverl = 0;

    // initialize the particle on a square
    init_random();

    // initialize Verlet list
    update_verlet_list();
    

}

void System::update_verlet_list()
{
    Nverl += 1;
    verletList = std::vector<std::vector<int> >(N);    
    for(int i=0; i<N; ++i) {
        particles_before_update[i] = particles[i];
        for(int j=i+1; j<N; ++j) {
            if( xyz::dist_sq_pbc(particles[i].r, particles[j].r,L) < rv ){
                verletList[i].push_back(j);
                verletList[j].push_back(i);
            }
        }
    }
}

void System::init_random()
{
	// node_pd: nodes per dim
	int node_pd = ceil(pow(1.*N,1./3));
	int Nnodes = node_pd*node_pd*node_pd;
	double  node_dist = L/node_pd;

	std::vector<XYZ> nodes(Nnodes);
	int i=0;
	for(int xi =0;xi<node_pd;++xi) {
		for(int yi=0;yi<node_pd;++yi) {
			for(int zi=0;zi<node_pd;++zi) {
				nodes[i].x = xi*node_dist;
				nodes[i].y = yi*node_dist;
				nodes[i].z = zi*node_dist;
				++i;
			}
		}
	}
    
    // shuffle nodes
    for(unsigned int i=0; i<N; ++i) {
        int j = i + (int)( rudist01()*(N-i) );
        XYZ temp = nodes[i];
        nodes[i] = nodes[j];
        nodes[j] = temp;
    }
	for(unsigned int i=0;i<N; ++i) {
        particles[i].r = nodes[i];
        particles[i].index = i;
        particles_before_update[i] = particles[i];
    }

}



#endif
