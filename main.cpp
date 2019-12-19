
#include "xyz.h"
#include "pair_correlation.h"
#include "potential.h"
#include "system.h"
#include "configFile.h"


#include <vector>
#include <iostream>

using namespace std;

int main()
{

	Config config("input.txt");

    int T = config.get_parameter<int>("T");
    int T_init = config.get_parameter<int>("T_init");
    int T_sample = config.get_parameter<int>("T_sample");
    int Tmc = config.get_parameter<int>("Tmc");
    int print_every = config.get_parameter<int>("print_every");

    int seed = config.get_parameter<int>("seed");
    int N = config.get_parameter<int>("N");
    double L = config.get_parameter<double>("L");
    double d = config.get_parameter<double>("d");

    double rhs = config.get_parameter<double>("rhs");
    double rc = config.get_parameter<double>("rc");
    double A = config.get_parameter<double>("A");
    double alpha = config.get_parameter<double>("alpha");
    double rv = config.get_parameter<double>("rv");

    int Nbin = config.get_parameter<int>("Nbin");

    double bs = (1.*L)/(1.*Nbin);

	PairCorrelation pc(N,L,Nbin,bs);
    Potential potential(A,alpha, rhs, rc);
    System system(seed, N, L,potential, d, rv);


    for(int t=0; t<T_init; ++t) {

        if(t%print_every == 0) {
            cout <<T+T_init << '\t' <<  t << endl;
            cout << system.Nacc/( (double) system.Ntry ) << endl << endl;
        }

        for(int tmc=0; tmc<Tmc; ++tmc) {
            //system.mc_move();
            system.mc_move_verlet();
        }
    }
    
    
    for(int t=0; t<T; ++t) {

        if(t%print_every == 0) {
            cout <<T+T_init << '\t' <<  t+T_init << endl;
            cout << system.Nacc/( (double) system.Ntry ) << endl << endl;
        }

        for(int tmc=0; tmc<Tmc; ++tmc) {
            //system.mc_move();
            system.mc_move_verlet();
        }

        if( t%T_sample == 0)
            pc.sample(system.particles);
    }
    cout <<( (double) system.Nacc)/( (double) system.Ntry) << endl;

    cout << system.Nverlet << endl;
    cout << system.Nacc/( (double) system.Ntry ) << endl;

    pc.normalize();
    pc.write("gr.dat");

    return 0;
}
