#include "xyz.h"

#include <iostream>

using namespace std;


int main()
{

    double L = 10.;
    XYZ r1(9,1,3);
    XYZ r2(-12,1,2);

    double d = xyz::dist_pbc(r1,r2,L);
    r2.pbc(L);
    cout << r2 << endl;

    cout << d << endl;
    return 0;
}


