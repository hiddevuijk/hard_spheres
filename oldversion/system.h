#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <vector>
#include <assert.h>
#include <boost/random.hpp>

#include "xyz.h"

class System{
public:

	System (unsigned int  N, double L,double beta, double rhs, 
            double alpha, double A, double deltaa, double rn,
             unsigned int seed);

    const boost::normal_distribution<double> ndist;
    const boost::uniform_real<double> udist;

    unsigned int seed;
    boost::mt19937 rng;
    boost::random::uniform_int_distribution<> ridist;
    boost::variate_generator<boost::mt19937&,
           boost::normal_distribution<double> > rndist;
    boost::variate_generator<boost::mt19937&,
            boost::uniform_real<double> > rudist; 


	void move();
	void move_nl();
	void initialize();

	int get_Naccepted() const {return Naccepted;}
	int get_Nmoves() const {return Nmoves;}
	int get_fraction_accepted_moves() const 
		{ return ( (double) Naccepted)/( (double) Nmoves); }
	const std::vector<XYZ>& get_positions() const	
		{ return positions;}

	void pbc();
	void neighbour_update();
//private:

	// system params
	unsigned int N; // number of disks
	double L;	// system size
    double beta; // 1/temperature
    double rhs; // hardsphere radius
    double alpha; // yukawa decay length
    double A; // yukawa amplitude
     

	// max step size
	double delta;
	double rn;
	unsigned int Nmoves;
	unsigned int Naccepted;	

	std::vector<XYZ> positions;
    std::vector<XYZ> positions_before_update;
	// neighbour data
	std::vector<std::vector<unsigned int> > neighbour_index;
	std::vector<unsigned int> neighbour_number;

	// use for trail moves
	XYZ new_position;
	// index of trial disk
	unsigned int index;	
	bool overlap;
    double dU; 
};

void System::move_nl() 
{
        index = ridist(rng);

        
        //system_func::xyz_random_uniform(new_position,rudist,delta);
        //new_position += positions[index];
        new_position.x = positions[index].x + ( 0.5-rudist() )*delta;
        new_position.y = positions[index].y + ( 0.5-rudist() )*delta;
        new_position.z = positions[index].z + ( 0.5-rudist() )*delta;

		XYZ temp;
		overlap = false;
        dU = 0;
		for(unsigned int i=0;i<neighbour_number[index];++i) {
			temp = positions[neighbour_index[index][i] ];

			if(xyz::dist_sq_pbc(new_position,temp,L) < 1. ){
				overlap = true;
				break;
			}
		}

		if( !overlap ) {
            if( rudist() < exp(-beta*dU) ) {
                positions[index] = new_position;
                ++Naccepted;
                // if particle moved to far, update neighbour list
                if( xyz::dist_pbc(positions[index],positions_before_update[index],L )
                        > 0.5*(rn-rhs)-delta )  {
                    neighbour_update();
                }
            } 
		}


		++Nmoves;
}

void System::neighbour_update()
{
	double dist_sq;
	std::fill(neighbour_number.begin(),
			neighbour_number.end(),0);
	for(unsigned int i=0;i<N;++i) {
		for(unsigned int j=i+1;j<N;++j) {
			dist_sq = xyz::dist_sq_pbc(positions[i],positions[j],L);
			if(dist_sq < 1){ std::cerr<< "FUCK \n";
                std::cerr << dist_sq << '\n';
                std::cerr << Nmoves << '\n';
            }
            assert(!(dist_sq<1));
			if(dist_sq < rn*rn) {
				neighbour_index[i][ neighbour_number[i] ] = j;
				neighbour_index[j][ neighbour_number[j] ] = i;
				++neighbour_number[i];
				++neighbour_number[j];
			}
		}
	}

    for(unsigned int i=0;i<N;++i)
        positions_before_update[i] = positions[i];
}

void System::pbc()
{
	for(unsigned int i=0;i<N;++i)
		positions[i].pbc(L);
}

System::System(unsigned int NN, double LL, double betaa, double rhss, double alphaa, double AA,
                 double dd, double rn, unsigned int seed)
: ndist(0,1), udist(0,1), seed(seed), rng(seed),ridist(0,NN-1),
    rndist(rng,ndist), rudist(rng, udist),
N(NN), L(LL), beta(betaa), rhs(rhss), alpha(alphaa), A(AA), delta(dd), rn(rn), Nmoves(0), Naccepted(0),
	positions(NN), positions_before_update(NN), neighbour_index(NN,std::vector<unsigned int>(NN,0)),
	neighbour_number(NN,0)
{}

void System::move() 
{

        //system_func::xyz_random_uniform(new_position,rudist,delta);
        //new_position += positions[index];
        new_position.x = positions[index].x + ( 0.5-rudist() )*delta;
        new_position.y = positions[index].y + ( 0.5-rudist() )*delta;
        new_position.z = positions[index].z + ( 0.5-rudist() )*delta;



		overlap = false;	
		for(unsigned int i=0;i<positions.size(); ++i) {
			if(i != index and 
				xyz::dist_sq_pbc(new_position,positions[i],L) < 1.) {
					overlap = true;
			}
			if(overlap) break;
		}

		if( !overlap ) {
			positions[index] = new_position;
			++Naccepted;
		}
		++Nmoves;
}

void System::initialize()
{


	unsigned int p_per_dim = floor(L);
	int k = 0;
	std::vector<XYZ> vertices(p_per_dim*p_per_dim*p_per_dim);
	std::vector<unsigned int> indices(p_per_dim*p_per_dim*p_per_dim);	
	for(unsigned int i=0;i<p_per_dim; ++i) {
		for(unsigned int j=0;j<p_per_dim; ++j) {
            for(unsigned int l=0;l<p_per_dim; ++l) {
                indices[ i*p_per_dim*p_per_dim + j*p_per_dim + l] = i*p_per_dim*p_per_dim + j*p_per_dim + l;	
                vertices[k].x = 0.5 + i;
                vertices[k].y = 0.5 + j;
                vertices[k].z = 0.5 + l;
                ++k;
            }
		}
	}

    

	// randomize indices
	unsigned int temp;
	unsigned int j;
	for(unsigned int i=0;i<indices.size()-1;++i ) {
		//j = i + ( rand() % (indices.size() - i) );
        j = i + (int) ( rudist()*( indices.size() - i ) );
		temp = indices[i];
		indices[i] = indices[j];
		indices[j] = temp;
	}
	// assign positions
	for(unsigned int i=0; i<N; ++i) {
		positions[i] = vertices[indices[i]];
        positions_before_update[i] = positions[i];
	}


}
#endif
