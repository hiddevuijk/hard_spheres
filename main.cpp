
#include "xyz.h"
#include "system.h"
#include "potential.h"
#include "pair_correlation.h"

#include <vector>
#include <iostream>

using namespace std;

int check_overlap(const vector<Particle>& particles, double rhs, Potential p, double L)
{
    int N = particles.size();
    int n = 0;
    for( int i=0;i<N;++i) {
        for(int j=i+1; j<N; ++j) {
            if( p.get_overlap(particles[i].r, particles[j].r,L) ) {
                ++n;
                cout << particles[i].r << endl;
                cout << particles[j].r << endl;
                cout << xyz::dist_sq_pbc(particles[i].r, particles[j].r,L) - 4*rhs*rhs << endl;
                cout << endl;
            }
        
        }
    }
    return n;
}

int main()
{

    int T = 50;
    int T_sample = 10;
    int seed = 123456789;
    int N = 1000;
    double L = 16.666667;
    double d = 0.5;
    cout << "rho " << 1.*N/(1.*L*L*L) << endl;

    double rhs = .5;
    double rc = 1.2;
    double beta = 1.;
    double rv = 5.;

    int Nbin = 100;
    double bs = (1.*L)/(1.*Nbin);
	PairCorrelation pc(N,L,Nbin,bs);
    Potential potential(rhs, rc, beta);
    System system(potential, N, L, d, rv, seed);

    for(int t=0; t<T; ++t) {
        cout << t << endl;
        for(int tmc=0; tmc<N; ++tmc) {
            system.mc_move();
            //check_overlap(system.particles, rhs, potential, L);
        }
        pc.sample(system.particles);
    }

    cout << check_overlap(system.particles, rhs, potential, L) << endl;
    cout <<  1.*system.Nacc/system.Ntry << endl;
    cout << 1.*system.Nverl/system.Ntry << endl;
    pc.normalize();
    pc.write("gr.dat");

    return 0;
}
