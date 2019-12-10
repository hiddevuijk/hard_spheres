#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <vector>
#include <assert.h>
#include <boost/random.hpp>

#include "xyz.h"

class System{
public:

	System (unsigned int  N, double L,double deltaa, double rn, unsigned int seed);
	System (unsigned int  N, double L,double deltaa, double rn);

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
};

void System::move_nl() 
{
        index = ridist(rng);

        
        //system_func::xyz_random_uniform(new_position,rudist,delta);
        //new_position += positions[index];
        new_position.x = positions[index].x + ( 1-2*rudist() )*delta;
        new_position.y = positions[index].y + ( 1-2*rudist() )*delta;
        new_position.z = positions[index].z + ( 1-2*rudist() )*delta;

		XYZ temp;
		overlap = false;
		for(unsigned int i=0;i<neighbour_number[index];++i) {
			temp = positions[neighbour_index[index][i] ];

			if(xyz::dist_sq_pbc(new_position,temp,L) < 1. ){
				overlap = true;
				break;
			}
		}

		if( !overlap ) {
			positions[index] = new_position;
			++Naccepted;
            // rn - 0.5 ?????????????? 
            if( xyz::dist_pbc(positions[index],positions_before_update[index],L )
                    > 0.5*(rn-.5)-delta )  {
                neighbour_update();
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

System::System(unsigned int NN, double LL, double dd, double rn, unsigned int seed)
: ndist(0,1), udist(0,1), seed(seed), rng(seed),ridist(0,NN-1),
    rndist(rng,ndist), rudist(rng, udist),
N(NN), L(LL), delta(dd), rn(rn), Nmoves(0), Naccepted(0),
	positions(NN), positions_before_update(NN), neighbour_index(NN,std::vector<unsigned int>(NN,0)),
	neighbour_number(NN,0)
{}

void System::move() 
{

        //system_func::xyz_random_uniform(new_position,rudist,delta);
        //new_position += positions[index];
        new_position.x = positions[index].x + ( 1-2*rudist() )*delta;
        new_position.y = positions[index].y + ( 1-2*rudist() )*delta;
        new_position.z = positions[index].z + ( 1-2*rudist() )*delta;



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
		j = i + ( rand() % (indices.size() - i) );
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
