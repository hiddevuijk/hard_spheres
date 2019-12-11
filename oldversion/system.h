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

    double U(const XYZ& r1, const XYZ& r2)
    {    double dist = xyz::dist_pbc(r1,r2,L);
         return -A*exp( -alpha*(dist-dhs)/dist );
    }
    double U(double r) { return -A*exp(-alpha*(r-dhs)/r); }

	int get_Naccepted() const {return Naccepted;}
	int get_Nmoves() const {return Nmoves;}
	int get_fraction_accepted_moves() const 
		{ return ( (double) Naccepted)/( (double) Nmoves); }
	const std::vector<XYZ>& get_positions() const	
		{ return positions;}

	void pbc();
	void neighbour_update();
	void neighbour_update_cell();
    void add_neighbour(int ni, int nj);

    bool check_overlap();
//private:

	// system params
	unsigned int N; // number of disks
	double L;	// system size
    double beta; // 1/temperature
    double rhs; // hardsphere radius
    double dhs; // hardsphere diameter = 2*rhs
    double alpha; // yukawa decay length
    double A; // yukawa amplitude
     

	// max step size
	double delta;
    double max_step;
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
        
        new_position.x = positions[index].x + ( 0.5-rudist() )*delta;
        new_position.y = positions[index].y + ( 0.5-rudist() )*delta;
        new_position.z = positions[index].z + ( 0.5-rudist() )*delta;

		XYZ temp;
		overlap = false;
        dU = 0;
		for(unsigned int i=0;i<neighbour_number[index];++i) {
			temp = positions[neighbour_index[index][i] ];
            double old_dist = xyz::dist_pbc(positions[index], temp,L);     
            double new_dist = xyz::dist_pbc(new_position, temp,L);     
			if( new_dist < dhs ){
				overlap = true;
				break;
			}
            dU += U(new_dist);
            dU -= U(old_dist);
		}

		if( !overlap ) {
            if( rudist() < exp(-beta*dU) ) {
                positions[index] = new_position;
                positions[index].pbc(L);
                ++Naccepted;
                // if particle moved to far, update neighbour list
                double max = 0.5*(rn-rhs-max_step);
                double dist = xyz::dist_pbc(positions[index],positions_before_update[index],L ) ;
                if( dist > max)  {
                    //neighbour_update();
                    neighbour_update_cell();
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
            //add_neighbour(i,j);
			dist_sq = xyz::dist_sq_pbc(positions[i],positions[j],L);
			if(dist_sq < 1){ std::cerr<< "FUCK \n";
                std::cerr << dist_sq << '\n';
                std::cerr << Nmoves << '\n';
            }
            assert( !( dist_sq < dhs*dhs ) );
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

class CellIndex {
public:
    int i,j,k;
};

void System::neighbour_update_cell()
{

	std::fill(neighbour_number.begin(),
			neighbour_number.end(),0);

    int n = 1;
    while( L/n > rn ) n++;
    --n;
    double l = L/n;
    
   
    int nmaxincell = N;
    std::vector<CellIndex> cell_indices(N);
    std::vector<std::vector<std::vector<std::vector<int> > > > particle_index_in_cell(n,
                        std::vector<std::vector<std::vector<int> > >(n,
                            std::vector<std::vector<int> >(n,
                                std::vector<int>(nmaxincell) ) ) );
    std::vector<std::vector<std::vector<int> > > n_particle_in_cell(n,
                        std::vector<std::vector<int> >(n,
                            std::vector<int>(n,0) ) );
    int i,j,k;
    for(unsigned int ni=0; ni<N; ++ni) {
        positions[ni].pbc(L);
        i = (int) (positions[ni].x/l);
        j = (int) (positions[ni].y/l);
        k = (int) (positions[ni].z/l);
        cell_indices[ni].i = i;
        cell_indices[ni].j = j;
        cell_indices[ni].k = k;

        unsigned int ncell = n_particle_in_cell[i][j][k];
        n_particle_in_cell[i][j][k]++;
        particle_index_in_cell[i][j][k][ncell]  = ni;
    }
    
    int i_pbc, j_pbc, k_pbc;
    unsigned int ncell;
    int nk;
    for(unsigned int ni=0;ni<N;++ni) {
        i = (int) (positions[ni].x/l);
        j = (int) (positions[ni].y/l);
        k = (int) (positions[ni].z/l);
       
        // this cell (0,0,0)
        i_pbc = i; 
        j_pbc = j; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
            nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           if( ni==nk ) continue;
           add_neighbour(ni,nk);
        }

        //  cell (0,1,0)
        i_pbc = i; 
        j_pbc = (j+1)%n; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (0,0,1)
        i_pbc = i; 
        j_pbc = j; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (0,1,1)
        i_pbc = i; 
        j_pbc = (j+1)%n; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (0,-1,0)
        i_pbc = i; 
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (0,0,-1)
        i_pbc = i; 
        j_pbc = j;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (0,-1,-1)
        i_pbc = i; 
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        // right cell (1,0,0)
        i_pbc = (i+1)%n;
        j_pbc = j; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,1,0)
        i_pbc = (i+1)%n;
        j_pbc = (j+1)%n; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,0,1)
        i_pbc = (i+1)%n;
        j_pbc = j; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,1,1)
        i_pbc = (i+1)%n;
        j_pbc = (j+1)%n; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,-1,0)
        i_pbc = (i+1)%n;
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,0,-1)
        i_pbc = (i+1)%n;
        j_pbc = j;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (1,-1,-1)
        i_pbc = (i+1)%n;
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        // right cell (-1,0,0)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = j; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,1,0)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = (j+1)%n; 
        k_pbc = k; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,0,1)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = j; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,1,1)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = (j+1)%n; 
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,-1,0)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = (k+1)%n; 

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,0,-1)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = j;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }

        //  cell (-1,-1,-1)
        i_pbc = i>0 ? i-1 : n-1;
        j_pbc = j>0 ? j-1 : n-1;
        k_pbc = k>0 ? k-1 : n-1;

        ncell = n_particle_in_cell[i_pbc][j_pbc][k_pbc];
        for(unsigned int nj=0; nj<ncell; ++nj) {
           nk = particle_index_in_cell[i_pbc][j_pbc][k_pbc][nj];
           add_neighbour(ni,nk); 
        }



	}



    for(unsigned int i=0;i<N;++i)
        positions_before_update[i] = positions[i];

}
void System::add_neighbour(int ni, int nj)
{ 
    double dist_sq = xyz::dist_sq_pbc(positions[ni],positions[nj],L);
    if(dist_sq < 1){ std::cerr<< "FUCK \n";
        std::cerr << dist_sq << '\n';
        std::cerr << Nmoves << '\n';
    }
    assert(!(dist_sq<1));
    if(dist_sq < rn*rn) {
        neighbour_index[ni][ neighbour_number[ni] ] = nj;
        neighbour_index[nj][ neighbour_number[nj] ] = ni;
        ++neighbour_number[ni];
        ++neighbour_number[nj];
    }
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
N(NN), L(LL), beta(betaa), rhs(rhss), dhs(2*rhss), alpha(alphaa), A(AA), delta(dd), max_step(sqrt(3*dd*dd) ),
         rn(rn), Nmoves(0), Naccepted(0),
	positions(NN), positions_before_update(NN), neighbour_index(NN,std::vector<unsigned int>(NN,0)),
	neighbour_number(NN,0)
{}

void System::move() 
{

        index = ridist(rng);
        //system_func::xyz_random_uniform(new_position,rudist,delta);
        //new_position += positions[index];
        new_position.x = positions[index].x + ( 0.5-rudist() )*delta;
        new_position.y = positions[index].y + ( 0.5-rudist() )*delta;
        new_position.z = positions[index].z + ( 0.5-rudist() )*delta;



		overlap = false;	
		for(unsigned int i=0;i<positions.size(); ++i) {
			if(i != index and 
				xyz::dist_sq_pbc(new_position,positions[i],L) < 2*rhs) {
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




bool System::check_overlap()
{
    double dist;
    for(int ni=0;ni<N;++ni) {
        for(int nj=ni+1;nj<N; ++nj) {
            dist =  xyz::dist_pbc(positions[ni], positions[nj], L);
            if(dist < dhs )
                return true;
        }
    }
    return false;



}
#endif
