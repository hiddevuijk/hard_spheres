#ifndef GUARD_POTENTIAL
#define GUARD_POTENTIAL

#include "xyz.h"

class Potential {
public:
    Potential( double rhs, double rco, double beta)
        : rhs(rhs), dhs(2*rhs), dhs_sq(4*rhs*rhs),
          rhs_sq(rhs*rhs), rco(rco), beta(beta)
        {}

    
    double U( const XYZ& r1, const XYZ& r2, double L) { return 0;}
    bool get_overlap(const XYZ& r1, const XYZ& r2, double L) ;

    double rhs,dhs, dhs_sq, rhs_sq, rco, beta;
};


bool Potential::get_overlap(const XYZ& r1, const XYZ& r2, double L)
{
    if( xyz::dist_sq_pbc(r1,r2,L) < dhs_sq ) return true;
    else return false;
}

#endif
