//////////////////////////////
//
// 3 dimensional vector class
// The +-*/ are defined as elementwise
// operations.
// 
//////////////////////////////

#ifndef GUARD_XYZ_H
#define GUARD_XYZ_H

#include <math.h>
#include <iostream>

class XYZ {
public:
	XYZ(double xx, double yy, double zz): x(xx), y(yy), z(zz) {}
	XYZ(): x(0.), y(0.), z(0.) {}

	// only two variable
	double x,y,z;

	// member functions	
	double length() const { return sqrt(x*x + y*y + z*z);}
	double length_sq() const {return x*x+y*y+z*z;}
	void normalize(double d=1.);
	void pbc(double L);
	XYZ operator+=(const XYZ& r);
	XYZ operator+=(const double& add);
	XYZ operator-=(const XYZ& r);
	XYZ operator-=(const double &sub);
	XYZ operator*=(const XYZ& r);
	XYZ operator*=(const double &mult);
	XYZ operator/=(const XYZ& r);
	XYZ operator/=(const double &div);
	
};

/////////////////////////
//
// NON MEMBER FUNCTIONS 
//
/////////////////////////


// out stream
std::ostream& operator<<(std::ostream &out, XYZ const& r);
// overload addition
XYZ operator+(const XYZ& l, const XYZ& r);
XYZ operator+(const XYZ& l, const double& add);
XYZ operator+(const double& add, const XYZ& r);
// overload subtraction
XYZ operator-(const XYZ& l, const XYZ& r);
XYZ operator-(const XYZ& l, const double& sub);
XYZ operator-(const double& sub, const XYZ& r);
// overload multiplication
XYZ operator*(const XYZ& l, const XYZ& r);
XYZ operator*(const XYZ& l, const double& mult);
XYZ operator*(const double& mult, const XYZ& r);
// overload division
XYZ operator/(const XYZ& l, const XYZ& r);
XYZ operator/(const XYZ& l, const double& div);


/////////////////////////
//
// NON MEMBER FUNCTIONS 
// IN NAME SPACE xyz
//
/////////////////////////

namespace xyz {
// distance between c1 and c2
double dist(const XYZ &c1, const XYZ &c2);

// same as dist, but with periodic boundary conditions
double dist_pbc(const XYZ &c1, const XYZ &c2, double L);


// returns the distance squared
double dist_sq(const XYZ &c1, const XYZ &c2);

// same as dist_sq, but with periodic boundary conditions
double dist_sq_pbc(const XYZ &c1, const XYZ &c2, double L);





// dot produc between c1 and c2
double dot(const XYZ &c1, const	XYZ &c2);

// the cross product, returns a double
XYZ cross(const XYZ &c1, const XYZ &c2);


};

#endif
