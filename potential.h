#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

#include "xyz.h"

class Potential {
public:
    Potential(double A, double alpha,  double rhs, double rco)
        : A(A), alpha(alpha), rhs(rhs), dhs(2*rhs), dhs_sq(4*rhs*rhs),
          rhs_sq(rhs*rhs), rco(rco) {}


    // Yukawa potiential    
    double U( double  r) { 
            if( r < rco ) {
                return -A*exp( -alpha*(r-dhs) )/r; 
            } else {
                return 0;
            }
    }

    double U( const XYZ& r1, const XYZ& r2, double L) { 
            double r = xyz::dist_pbc(r1,r2,L);
            return U(r);
    }

    bool get_overlap(const XYZ& r1, const XYZ& r2, double L)
        { return xyz::dist_sq_pbc(r1, r2, L) < dhs_sq; }

    bool get_overlap(double r) { return r < dhs; }

    double A, alpha;
    double rhs,dhs, dhs_sq, rhs_sq, rco;
};




#endif
