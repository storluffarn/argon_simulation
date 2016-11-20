
// Molecular dynamics simulation of argon liquid n = 864 @ T = 120 K, for simulation methods
// Originally performed by A. Rahman, Phys letters, vol 136:2A
//
// SI-units
//

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//
// Preliminaries

// includes

#include <iostream>                                             // basic io stuff
#include <iterator>
#include <fstream>                                              // writing to file
#include <string>
#include <cmath>                                                // for math functions
#include <random>                                               // for random generation
#include <numeric>												// for accumulate
#include <algorithm>                                            // for count_if
#include <vector>                                               // for using the vector type members
#include <chrono>                                               // for measuring execution times
#include <ctime>                                                // for random seed time

using namespace std;                                            // tidy up

// constants

static const double umass = 1.660539040e-27;                    // atomic mass
static const double argonmass = umass*39.948;                   // argon mass
static const double resargonmass = 1/argonmass;                 // resiprocal argon mass
static const double pi = atan(1)*4;                             // for using pi
static const double respi = atan(1)*4;                          // resiprocal pi
static const double kB = 1.38064852e-23;                        // Boltzmann's constant
static const double sigma = 3.4e-10;                            // interaction term
static const double epsilon = 120*kB;                           // energy term
static const double dcutoff = 2.25*sigma;                       // interaction cutoff 
static const double dcutoffsquare = pow(dcutoff,2);             // interaction cutoff 

// objects

class dataclass													// data to be analysed
{
	// private section

	// data members
	
	string filename;											// name of file to read
	vector <double> data;										// data read
	int size;													// numbers written/read

	// function members
	
	vector<double> readdata(string);							// function for reading data
	void writedata(vector<double>*,string);						// function for writing data
	
	// public section
	
	public:
	
	// constructor
	dataclass(string s)
		: filename(s)											// setting initial values
	{
		data = readdata(filename);								// read data from file
		size = data.size();										// save size of data
	}
	
	//function members
	
	// inline functions
	double calcmean()											// calculates mean of data
	{
		return accumulate(data.begin(),data.end(),0.0)/data.size();	
	}
	
	// out of line functions	
	void acorrelation();										// auto correlation function
	void speeddist(double,int);									// speed distribution
	void radialdistfunc(double,int,int);						// radial distribution function
	
};

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
//
// Functions on files

vector <double>  dataclass::readdata(string filename)			// reads data from file
{
	ifstream readdata(filename);								// read from start to end to vector
	istream_iterator<double> start(readdata), end;
	vector<double> data(start, end);
	cout << "Read " << filename << " " << data.size() << " numbers" << endl;	// inform user of read date

	return data;
}

void dataclass::writedata(vector <double>* writedata, string filename)	// wrotes data to file
{
	ofstream filestream;
	filestream.open(filename);

	//int it = 0;
	for(unsigned int k = 0; k < writedata->size(); k++)			// write all data to file
	{
		filestream << (*writedata)[k] << endl;						
		//it++;
	}
	
	filestream.close();

	cout << "Wrote " << writedata->size() << " numbers to " << filename << endl;	// inform user of what happens
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
//
// Statstics

void dataclass::speeddist(double maxspeed, int bins)			// calculates the velocity distrobution
{
	vector <double> speeddist(bins);
	double resbins = 1/(static_cast<double> (bins));			// reciprocal velbins

	double shell1 = 0, shell2 = 0;
	for(double k = 0; k < bins; k++)							// count speeds into bins
	{ 
		shell2 = (k + 1)*resbins*maxspeed;
		speeddist[k] = count_if(data.begin(), data.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;});
		shell1 = shell2;
	}

	writedata(&speeddist,"speeddist.txt");
}

void dataclass::radialdistfunc(double boxlength, int natoms, int bins)       // calculates the radial distribution
{
	vector <double> rdf (bins);
	double resbins = 1/(static_cast<double> (bins));			// reciprocal densbins
	double halfbox = 0.5*boxlength;								// half box length

	vector<int>hist(bins);
	double shell1 = 0, shell2 = 0;

	for(double k = 0; k < bins; k++)							// bin distances up to half box length
	{ 
		shell2 = (k + 1)*resbins*halfbox;
		hist[k] = count_if(data.begin(), data.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;}); 
		shell1 = shell2;
	}

	double dens, norm;											
	double natomssqr = pow(natoms,2);							// number of atoms squared
	double resvol = pow(boxlength,-3);							// resiprocal volume

	for (double k = 0; k < bins; k++)							// calculates radial distribution function
	{
		dens = natomssqr*resvol;
		norm = 4*pi*pow((k+1)*resbins*halfbox,2)*dens*(halfbox*resbins);
		rdf[k] = hist[k]/norm;
	}	

	writedata(&rdf,"radialdistfunc.txt");
}

void dataclass::acorrelation()
{
	unsigned int tau = data.size();
	vector <double> corrlist(tau);

	double norm, corr;
	double mean = calcmean();									// mean of data
	
	for (unsigned int k = 0; k < tau; k++)						// calculate all autocorrelations
	{

		norm = 0;
		corr = 0;
				
		for (unsigned int l = 0; l < data.size(); l++)
			norm += pow(data[l]-mean,2);
			
		for ( unsigned int l = 0; l < data.size() - k; l++)
		{
			corr += (data[l] - mean)*(data[l+k] - mean);
		}
		
		corrlist[k] = corr/norm;
	
	}

	writedata(&corrlist, "correlation.txt");
}

// main

int main()
{
	// correlations
	dataclass particle("particle.txt");
	particle.acorrelation();

	dataclass kinetic("kinetic.txt");

	// speed distribution
	dataclass speeds("speeds.txt");
	
	int bins = 100;
	double maxspeed = 800;
	speeds.speeddist(maxspeed, bins);

	// radial distribution function
	dataclass rdf("distances.txt");
	
	bins = 500;
	double boxlength = 10.229*(3.4e-10);
	int npart = 864;
	rdf.radialdistfunc(boxlength, npart, bins);

}














































