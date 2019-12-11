#include "xyz.h"


/////////////////////////
//
// MEMBER FUNCTIONS
//
/////////////////////////

void XYZ::normalize(double d)
{	 double len = this->length();
	x *= d/len; y *= d/len; z*= d/len;
}

void XYZ::pbc(double L)
{
    x -= L*floor(x/L);
    y -= L*floor(y/L);
    z -= L*floor(z/L);
}

XYZ XYZ::operator+=(const XYZ& r)
	{ x += r.x; y += r.y; z += r.z; return *this; }

XYZ XYZ::operator+=(const double& add)
	{ x += add; y += add; z+= z; return *this; }

XYZ XYZ::operator-=(const XYZ& r)
	{ x -= r.x; y -= r.y; z -= r.z; return *this; }

XYZ XYZ::operator-=(const double& sub)
	{ x -= sub; y -= sub; z -= sub; return *this; }

XYZ XYZ::operator*=(const XYZ& r)
	{ x *= r.x; y *= r.y; z *= r.z; return *this; }

XYZ XYZ::operator*=(const double& mult)
	{ x *= mult; y *= mult; z *= mult; return *this; }

XYZ XYZ::operator/=(const XYZ& r)
	{ x /= r.x; y /= r.y; z /= r.z; return *this; }

XYZ XYZ::operator/=(const double& div)
	{ x /= div; y /= div; z /= div; return *this; }


	
/////////////////////////
//
// NON MEMBER FUNCTIONS 
//
/////////////////////////


// out stream
std::ostream& operator<<(std::ostream &out, XYZ const& r)
{ out << r.x << '\t' << r.y << '\t' << r.z; return out; }

// overload addition
XYZ operator+(const XYZ& l, const XYZ& r)
{ return XYZ(l.x+r.x, l.y+r.y, l.z+r.z); }

XYZ operator+(const XYZ& l, const double& add)
{ return XYZ(l.x+add, l.y+add, l.z+add); }

XYZ operator+(const double& add, const XYZ& r)
{ return XYZ(r.x+add, r.y+add, r.y+add); }

// overload subtraction
XYZ operator-(const XYZ& l, const XYZ& r)
{ return XYZ(l.x-r.x, l.y-r.y, l.z-r.z); }

XYZ operator-(const XYZ& l, const double& sub)
{ return XYZ(l.x-sub, l.y-sub, l.z-sub); }

XYZ operator-(const double& sub, const XYZ& r)
{ return XYZ(sub-r.x,sub-r.y, sub-r.z); }

// overload multiplication
XYZ operator*(const XYZ& l, const XYZ& r)
{ return XYZ(l.x*r.x, l.y*r.y, l.z*r.z); }

XYZ operator*(const XYZ& l, const double& mult)
{ return XYZ(l.x*mult, l.y*mult, l.z*mult); }

XYZ operator*(const double& mult, const XYZ& r)
{ return XYZ(mult*r.x,mult*r.y,mult*r.z); }


// overload division
XYZ operator/(const XYZ& l, const XYZ& r)
{ return XYZ(l.x/r.x, l.y/r.y, l.z/r.z); }

XYZ operator/(const XYZ& l, const double& div)
{ return XYZ(l.x/div, l.y/div, l.z/div); }




/////////////////////////
//
// NON MEMBER FUNCTIONS 
// IN NAME SPACE xyz
//
/////////////////////////

// distance between c1 and c2
double xyz::dist(const XYZ &c1, const XYZ &c2)
{ return sqrt( (c1.x - c2.x)*(c1.x - c2.x) +
			(c1.y - c2.y)*(c1.y - c2.y) +
            (c1.z - c2.z)*(c1.z - c2.z) ); }

// same as dist, but with periodic boundary conditions
double xyz::dist_pbc(const XYZ &c1, const XYZ &c2, double L)
{
	XYZ d( c1.x-c2.x, c1.y-c2.y, c1.z-c2.z );

	d.x -= L*round(d.x/L);
	d.y -= L*round(d.y/L);
    d.z -= L*round(d.z/L);

	return sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}


// returns the distance squared
double xyz::dist_sq(const XYZ &c1, const XYZ &c2)
{ return  (c1.x - c2.x)*(c1.x - c2.x) +
		  (c1.y - c2.y)*(c1.y - c2.y) + 
		  (c1.z - c2.z)*(c1.z - c2.z);
}

// same as dist_sq, but with periodic boundary conditions
double xyz::dist_sq_pbc(const XYZ &c1, const XYZ &c2, double L)
{

	XYZ d( c1.x-c2.x, c1.y-c2.y, c1.z-c2.z );

	d.x -= L*round(d.x/L);
	d.y -= L*round(d.y/L);
    d.z -= L*round(d.z/L);

	return d.x*d.x + d.y*d.y + d.z*d.z;
}


// dot produc between c1 and c2
double xyz::dot(const XYZ &c1, const	XYZ &c2) 
{ return c1.x*c2.x + c1.y*c2.y + c1.z*c2.z; }

// the cross product, returns a double
XYZ xyz::cross(const XYZ &c1, const XYZ &c2) 
{
    XYZ temp( c1.y*c2.z - c1.z*c2.y,
              c1.z*c2.x - c1.x*c2.z,
              c1.x*c2.y - c1.y*c2.x );

    return temp;
}


