#ifndef GUARD_NEIGHBOUR_CHECK_H
#define GUARD_NEIGHBOUR_CHECK_H

#include "xyz.h"
#include "particle.h"

#include <vector>

class NeighbourCheck {
public:
    NeighbourCheck( double l, double L)
        :l(l), L(L) {}

    void first(const std::vector<Particle>& particles);
    void second(const std::vector<Particle>& particles);

    double l;
    double L;
    std::vector<std::vector<int> > matrix1;
    std::vector<std::vector<int> > matrix2;
};





#endif
