
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
	double average;

	// function members
	
	vector<double> readdata(string);							// function for reading data
	void writedata(vector<double>*,string);						// function for writing data
	
	// public section
	
	friend void radialdistfunc(double,int,int,dataclass*);										// radial distribution function	
	friend void heatcapacity(dataclass*,dataclass*,dataclass*);
	friend void speeddist(double,int,dataclass*);
	friend void blockaverage(double, dataclass*);

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

	void calcavg()											// calculates mean of data
	{
		average = accumulate(data.begin(),data.end(),0.0)/data.size();	
	}

	vector <double> getdata()
	{
		return data;
	}	
	
	// out of line functions	
	void acorrelation();										// auto correlation function
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

// Out of line functions

void dataclass::acorrelation()
{
	unsigned int tau = data.size();
	vector <double> corrlist(tau);

	double norm, corr;
	
	for (unsigned int k = 0; k < tau; k++)						// calculate all autocorrelations
	{

		norm = 0;
		corr = 0;
				
		for (unsigned int l = 0; l < data.size(); l++)
			norm += pow(data[l]-average,2);
			
		for ( unsigned int l = 0; l < data.size() - k; l++)
		{
			corr += (data[l] - average)*(data[l+k] - average);
		}
		
		corrlist[k] = corr/norm;
	
	}

	writedata(&corrlist, "correlation.txt");
}

// Non class member functions
	
vector <double> calcerror(vector <double>* data)
{
	double average = accumulate(data->begin(),data->end(),0.0)/data->size();	
	
	vector <double> shifted (data->size());

	for (auto& el: shifted)
	{
		auto k = &el - &shifted[0];

		el = pow((*data)[k] - average,2);
	}

	double variance = accumulate(shifted.begin(),shifted.end(),0.0)/shifted.size();
	double stdev = sqrt(variance);
	double error = stdev/sqrt(data->size()-1);

	vector <double> retvals = {average,variance,stdev,error};
	
	return retvals;
}

void speeddist(double maxspeed, int bins, dataclass* speeds)			// calculates the velocity distrobution
{
	vector <double> speeddist(bins);
	vector <double> speederror(bins);
	vector <double> binerror(bins);
	vector <double> binsizes(bins);	

	int binsize = maxspeed/static_cast<double> (bins);

	double resbins = 1/(static_cast<double> (bins));			// reciprocal velbins

	double shell1 = 0, shell2 = 0;
	for(double k = 0; k < bins; k++)							// count speeds into bins
	{ 
		shell2 = (k + 1)*resbins*maxspeed;
		//speeddist[k] = count_if(speeds->data.begin(), speeds->data.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;});

		int count = 0;
		vector <double> speedbin = {};

		for (auto& el : speeds->data)
		{
			if (shell1 < el && el < shell2)
			{
				speedbin.push_back(el);
				count++;
			}
		}

		vector <double> stats = calcerror(&speedbin);

		speeddist[k] = count;
		speederror[k] = stats.back()*binsize;
		if (count == 0 || count == 1)
			speederror[k] = 0;
		
		if (count != 0)
			binerror[k] = sqrt(count)/count;

		binsizes[k] = k*binsize;

		shell1 = shell2;
	}
	
	speeds->writedata(&speeddist,"speeddist.txt");
	speeds->writedata(&speederror,"speederror.txt");
	speeds->writedata(&binerror,"spdbinerror.txt");
	speeds->writedata(&binsizes,"spdbinsize.txt");
}

void radialdistfunc(double boxlength, int natoms, int bins, dataclass* distances)       // calculates the radial distribution
{
	vector <double> rdf (bins);
	vector <double> rdferror(bins);
	vector <double> binerror(bins);
	vector <double> binsizes(bins);	

	double halfbox = 0.5*boxlength;								// half box lengt
	double binsize = halfbox/static_cast<double> (bins);
	double resbins = 1/(static_cast<double> (bins));			// reciprocal densbins
	
	double shell1 = 0, shell2 = 0;

	for(double k = 0; k < bins; k++)							// bin distances up to half box length
	{ 
		shell2 = (k + 1)*resbins*halfbox;
		
		int count = 0;
		vector <double> rdfbin = {};

		for (auto& el : distances->data)
		{
			if (shell1 < el && el < shell2)
			{
				rdfbin.push_back(el);
				count++;
			}
		}

		vector <double> stats = calcerror(&rdfbin);

		rdf[k] = count;
		rdferror[k] = stats.back(); //*binsize;
		if (count == 0 || count == 1)
			rdferror[k] = 0;
		if (count != 0)
			binerror[k] = sqrt(count)/count;

		binsizes[k] = k*binsize;
		
		
		shell1 = shell2;
	}

	double dens, norm;											
	double natomssq = pow(natoms,2)*0.5;							// number of atoms squared
	double resvol = pow(boxlength,-3);							// resiprocal volume

	for (double k = 0; k < bins; k++)							// calculates radial distribution function
	{
		dens = natomssq*resvol;
		norm = 4*pi*pow((k+1)*resbins*halfbox,2)*dens*(halfbox*resbins);
		rdf[k] /= norm;
		binerror[k] /= norm;
	}	

	distances->writedata(&rdf,"radialdistfunc.txt");
	distances->writedata(&rdferror,"rdferror.txt");
	distances->writedata(&binerror,"rdfbinerror.txt");
	distances->writedata(&binsizes,"rdfbinsize.txt");
}

void blockaverage(double blocks, dataclass* datac)
{
	double average = accumulate(datac->data.begin(),datac->data.end(),0.0)/datac->data.size();

	double blocksize = datac->data.size()/blocks;
	double blockavg = 0;
	double accavg = 0;
	int pivot = 0;

	for (unsigned int k = 0; k < blocks; k ++)
	{
		for (int l = pivot; l < pivot + blocksize; l++)
		{
			accavg += pow(datac->data[l],2);
		}
		blockavg += accavg/blocksize-pow(average,2);
		cout << accavg/blocksize << " " << blockavg << endl;
		accavg = 0;
		pivot += blocksize;
	}

	double blockerror = sqrt(blockavg/(blocks-1.0));

	cout << "blockerror: " << blockerror << endl;
}

// main

int main()
{
	// correlations
	dataclass particle("particle.txt");
	particle.acorrelation();

	//kinetic energy
	dataclass kinetic("kinetic.txt");

	//vector <double> tmp = {1000,500,250,100,50,10,5,2};

	//for (unsigned int k = 0; k < tmp.size(); k++)
	//	blockaverage(tmp[k], &kinetic);
	blockaverage(10, &kinetic);

	vector <double> kindata = kinetic.getdata();
	vector <double> kinerror = calcerror(&kindata);
	cout << kinerror[2] << endl;

	// speed distribution
	dataclass speeds("speeds.txt");
	
	int bins = 100;
	double maxspeed = 500;
	speeddist(maxspeed, bins, &speeds);

	// radial distribution function
	dataclass rdf("distances.txt");
	
	bins = 500;
	double boxlength = 10.229*(3.4e-10);
	int npart = 864;
	radialdistfunc(boxlength, npart, bins, &rdf);

}














































