
#include "pair_correlation.h"



PairCorrelation::PairCorrelation(unsigned int NN, double LL, 
		unsigned int Nbinn, double bss) :
			N(NN), L(LL),  Nbin(Nbinn),
			bs(bss), bins(Nbinn), pc(Nbinn,0.), Nsamp(0),is_normalized(false)
{

	// set mid points of the bins
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (0.5+i)*bs;
}			

void PairCorrelation::sample( const std::vector<Particle> &particles )
{

	// only sample before normalization 
	if(is_normalized){
		// quit
	}

	double d;

	for(unsigned int i=0;i<N;++i) {
		for(unsigned int j=i+1;j<N;++j) {
			//d = xyz::dist_pbc(particles[i].r,particles[j].r,L);
			d = xyz::dist_pbc(particles[i].r,particles[j].r,L) - 1;
			if(d< 1.*L) pc[ (int)(d/bs) ] += 2.;
		}
	}

	++Nsamp;

}
void PairCorrelation::normalize()
{
	const double pi = 4*atan(1.);
	// number of ideal gas particles in bin.
	double Nid;
	// bulk density
	double rhob = N/(L*L*L);
  
	for(unsigned int i=0;i<Nbin;++i ) {
        //double vb = ( (i+1)*(i+1)*(i+1) - i*i*i) * bs*bs*bs;
        double vb = (1+(1+i)*bs)*(1+(1+i)*bs)*(1+(1+i)*bs) - (1+i*bs)*(1+i*bs)*(1+i*bs) ; 
        Nid = (4./3.)*pi*vb*rhob;
		pc[i] /= N*Nid*Nsamp; 
	}	

	is_normalized = true;
}


void PairCorrelation::write(std::ostream & out)
{
	for(unsigned int i=0;i<Nbin;++i) {
		out << bins[i] << ' ';
		out << pc[i] << '\n';
	}

}


void PairCorrelation::write(std::string outname)
{

	std::ofstream out;
	out.open(outname);
	for(unsigned int i=0;i<Nbin;++i) {
		out << bins[i] << ' ' << pc[i];
		if( i < Nbin-1) out <<  '\n';
	}
	out.close();
	
}
