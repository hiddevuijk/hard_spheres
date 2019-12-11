
#include "xyz.h"
#include "pair_correlation.h"
#include "potential.h"
#include "system.h"

#include <vector>
#include <iostream>

using namespace std;

int main()
{

    int T = 50000;
    int T_init = 20000;
    int T_sample = 1;
    int Tmc = 1000;
    int print_every = 50;

    int seed = 123456789;
    int N = 500;
    double L = 20.;
    double d = 0.1;

    double rhs = .5;
    double rc = 3.;
    double A = 0.;
    double alpha = 0;
    double rv = 5;

    int Nbin = 500;
    double bs = (1.*L)/(1.*Nbin);

	PairCorrelation pc(N,L,Nbin,bs);
    Potential potential(A,alpha, rhs, rc);
    System system(seed, N, L,potential, d, rv);


    for(int t=0; t<T_init; ++t) {

        if(t%print_every == 0) cout << t << endl;

        for(int tmc=0; tmc<Tmc; ++tmc) {
            //system.mc_move();
            system.mc_move_verlet();
        }
    }
    
    
    for(int t=0; t<T; ++t) {

        if(t%print_every == 0) cout << t << endl;

        for(int tmc=0; tmc<Tmc; ++tmc) {
            //system.mc_move();
            system.mc_move_verlet();
        }

        if( t%T_sample == 0)
            pc.sample(system.particles);
    }
    cout <<( (double) system.Nacc)/( (double) system.Ntry) << endl;

    cout << system.Nverlet << endl;

    pc.normalize();
    pc.write("gr.dat");

    return 0;
}
