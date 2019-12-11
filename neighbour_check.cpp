#include "neighbour_check.h"

void NeighbourCheck::first(const std::vector<Particle>& particles)
{
    int N = particles.size(); 

    for(int i=0;i<N; ++i) {
        for(int j=0; j<N; ++j) {
            if( i == j) continue;
            if( xyz::dist_pbc(particles[i].r,
                            particles[j].r, L) ) {

                matrix1[i].push_back(j);
            }
        }
    }


}


void NeighbourCheck::second(const std::vector<Particle>& particles)
{

    int N = particles.size(); 

    for(int i=0;i<N; ++i) {
        for(int j=0; j<N; ++j) {
            if( i == j) continue;
            if( xyz::dist_pbc(particles[i].r,
                            particles[j].r, L) ) {

                matrix2[i].push_back(j);
            }
        }
    }


}
