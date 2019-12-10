
#include "xyz.h"
#include "system.h"
#include "pair_correlation.h"
#include "configFile.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

int main()
{ 

	Config config("input.txt");
	double rn = config.get_parameter<double>("rn");
	unsigned int N = config.get_parameter<unsigned int>("N");
	double L = config.get_parameter<double>("L");
    double rhs = config.get_parameter<double>("rhs");
    double beta = config.get_parameter<double>("beta");
    double alpha = config.get_parameter<double>("alpha");
    double A = config.get_parameter<double>("A");
    cout << N/(L*L*L) << endl;

	// mc params
	double d = config.get_parameter<double>("d");
	unsigned int t_unit = config.get_parameter<unsigned int>("multiplier");
	unsigned int T_init = config.get_parameter<unsigned int>("T_init");
	unsigned int T = config.get_parameter<unsigned int>("T");
	unsigned int T_sample = config.get_parameter<unsigned int>("T_sample");
	unsigned int T_print = config.get_parameter<unsigned int>("T_print");
    unsigned int seed = config.get_parameter<unsigned int>("seed");
    
    string name = config.get_parameter<string>("name");

	// pc params
	unsigned int Nbin = config.get_parameter<unsigned int>("Nbin");
	double bs = 1.*L/Nbin;
	
	PairCorrelation pc(N,L,Nbin,bs);

	System system(N,L,beta,rhs,alpha, A, d,rn,seed);
	system.initialize();
	system.neighbour_update();
	for(unsigned int ti = 0;ti<T_init; ++ti) {
        // print status
		if( (ti%T_print) ==0){
			cout << (T_init+T) << '\t' << ti << endl;
			cout << ( (double) system.get_Naccepted() ) /
                    ( (double) system.get_Nmoves() ) << endl;
			cout << endl;
		}

        // make moves using neighbour list
		for( unsigned int tti=0;tti<t_unit;++tti)  
			system.move_nl();

	}	
    
	for(unsigned int ti = 0;ti<T; ++ti) {
        // print status
		if( (ti%T_print) ==0) {
			cout << (T_init+T) << '\t' << ti+T_init << endl;
			cout << ( (double) system.get_Naccepted() ) /
                    ( (double) system.get_Nmoves() ) << endl;
			cout << endl;
		}

        // make moves using neighbour list
		for(unsigned int tti=0;tti<t_unit;++tti)
			system.move_nl();

		//sample pair correlation
		if(  (ti%T_sample) == 0 )
			pc.sample(system.get_positions());

	}

    // print fraction of accepted moves
	cout << ( (double) system.get_Naccepted() ) /
			( (double) system.get_Nmoves() ) << endl;

	// normalize pair correlation function
	pc.normalize();

	// write pair correlation function
	pc.write(name+"_pc.dat");	


	return 0;
}


