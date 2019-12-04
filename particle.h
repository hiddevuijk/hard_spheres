#ifndef GUARD_PARTICLE_H
#define GUARD_PARTICLE_H

#include "xyz.h"
#include "cell_index.h"

class Particle {
public:
    // position
    XYZ r;

    // particle index
    int index;
    

    // triplet of indices corresponding to the cell the particle is in
    Cindex cell_index;
};


#endif
