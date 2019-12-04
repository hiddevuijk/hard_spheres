#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H


#include "xyz.h"
#include "particle.h"
#include "cell_index.h"

#include <vector>


class System {
public:

    
    // attemt to move particle i in the dr 
    bool mc_move(int i, XYZ dr);


    // set cell indices of particle p
    void get_cell_index(Particle& p);
    // give the particles the right cell index
    void assign_cell_index();

    void set_cell_neighbors();
    
    std::vector<int> particles;
    //std::vector<list> cells;
    std::vector<Cindex> neighbor_cells;

    int N;      // number of particles
    double L;   // system size
    double Lcell;   // cell size
    double Ncell;   // number of cells
    double ncell; // ncell^3 = Ncell
    double co;  // cut off length
    double d;  // move size 

    long int Nmove; // number of attempted moves
    long int Nacc;  // number of accepted moves
};


bool System::mc_move(int i, XYZ dr)
{

    double delta_E = 0;

    // loop over neighbor cells
    for(int nci = 0; nci < Ncell; ++nci) {
        // loop over particles in cell nci

        // particle_ptr = 
        // while( particle_ptr != end ){
        //      delta_E += U(particles[i], particles[particle_ptr->index])
        //              - U(particles[i]+dr, particles[particle_ptr->index]);
        // }

        
    } 
}


void System::set_cell_neighbors()
{
    
    neighbor_cells = std::vector<Cindex>(Ncell);
    int i=0;
    int j=0;
    int k=0;
    for(int ni = 0; ni<Ncell; ++ni) {
        neighbor_cells[ni] = Cindex(i,j,k);
        i += 1;
        if( i == ncell ) {
            i = 0;
            j += 1;
        } 
        if( j == ncell ) {
            j = 0;
            k += 1;
        }
    }

}


#endif
