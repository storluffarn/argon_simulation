// 
// Molecular dynamics simulation of argon liquid n = 864 @ T = 120 K, for simulation methods
// Originally performed by A. Rahman, Phys letters, vol 136:2A
//
// SI-units
//
//
// Bad things that should be fixed/thought about:
//
//		x save all particle data every step, do all data analysis in different file
//		x get code nicer plotting rutine, possibly use SVC standard or such
//
//
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//
// Preliminaries

// includes

#include <iostream>												// basic io stuff
#include <fstream>												// writing to file
#include <cmath>												// for math functions
#include <random>												// for random generation
#include <algorithm>											// for count_if
#include <vector>												// for using the vector type members
#include <chrono>												// for measuring execution times
#include <ctime>												// for random seed time
#include <string>

using namespace std;											// tidy up

// constants

static const double umass = 1.660539040e-27;					// atomic mass
static const double argonmass = umass*39.948;					// argon mass
static const double resargonmass = 1/argonmass;					// resiprocal argon mass
static const double pi = atan(1)*4;								// for using pi
static const double respi = atan(1)*4;							// resiprocal pi
static const double kB = 1.38064852e-23;						// Boltzmann's constant
static const double sigma = 3.4e-10;							// interaction term
static const double epsilon = 120*kB;							// energy term
static const double dcutoff = 2.25*sigma;						// interaction cutoff 
static const double dcutoffsquare = pow(dcutoff,2);				// interaction cutoff 
static const double potjump = 4*epsilon*(-6*pow(sigma,6)*pow(dcutoffsquare,-4) + 12*pow(sigma,12)*pow(dcutoffsquare,-7) );


// classes

class cloud														// class to hold the set of all particles and functions on these
{
	// private section
	
class particle													// class to hold particle data and functions
{
	// private members

	// data members
	double pox;													// positions
	double poy;
	double poz;
	double vex;													// velocities
	double vey;
	double vez;
	double fox;													// forces
	double foy;
	double foz;
	double spd;													// speed
	double spdsqr;												// speed squared
	double kin;													// potential energy
	double pot;													// kinetic energy

	// member functions (empty)

	// freinds with bennifits
	friend cloud;												// to let cloud access private members
	friend vector <double> uglyerror(cloud*);								// this should go in the data analysis program, but alas, it does not =(

	// public members
	public:
	
	// public function members
	void posverlet(double,particle*);							// verlet velocity function
	void velverlet(double,particle*);							// vertet position function
	void langevin(double,double,double,particle*);				// for langevin dynamics
	void kinetic(particle*);									// calculating kinetic energy
	void speed(particle*);										// calculating speed
	void speedsquare(particle*);								// calculating the speed square to save CPU time	
};

	// data members
	
	// box data
	int natomsside;												// number of atoms on one side
	int natoms;													// total number of atoms
	int pairs;													// total number of particle pairs
	double boxlength;											// box side length
	double resboxlength;										// resiprocal boxlength
	double resnatoms;
	double gridsize;											// initial separation of particles
//	double volume;												// total volume of cloud
//	double density;												// density of cloud
//	double partdens;											// particle density of cloud
	double langfric;											// langevin friction constant
	
	// object member
	vector<particle> particleArray{static_cast<unsigned int> (natoms)};		// vector of all particles

	// function members
	
	double calcpairs()											// calculate the number of pairs
	{
		int sum = 0.5*(natoms - 1)*(natoms);					// arithemtic sum
		return sum;
	}
//	double calcvolume()											// caclulate volume
//	{
//		return pow(boxlength,3);
//	}
//
//	double calcdensity()										// calculate density
//	{
//		return natoms * argonmass / volume;
//	}
//
//	double calcpartdens()										// calculate number density
//	{
//		return natoms / volume;
//	}
//	
//	void calcavgkin()											// calculate average kinetic energy
//	{
//		avgkin = totkin/natoms;
//	}
//
//	void calcavgpot()											// calculate average potential energy
//	{
//		avgpot = totpot/natoms;
//	}
	
// friends										

	friend vector <double> uglyerror(cloud*);								// this should go in the data analysis program, but alas, it does not =(

// public section
	
	public:
	
	// data members

	// cloud data
	double temp;												// temperature
	double heatcapacity;
	double totkin;												// total kinetic energy
	double totpot;												// total potential energy

	// constructor
	cloud(int t1, unsigned int t2, double t3, double t4, double t5, double t6)			
		: natomsside(t1), natoms(t2), boxlength(t3), gridsize(t4), langfric(t5), temp(t6)		// setting initial conditions
	{															
		// constructing functions goes here
		resboxlength = 1.0/boxlength;
		resnatoms = 1.0/natoms;	
//		pairs = calcpairs();
//		volume = calcvolume();
//		density = calcdensity();
//		partdens = calcpartdens();
	}	
	
	// function members

	// accessors
	void setpos(int,double,double,double);				// for setting positions
	void setvel(int,double,double,double);				// for setting velcities
	void setforce(int,double,double,double);				// for setting forces
	void setpot(int,double);								// for setting potential energy
	vector<double> getpos(int);							// for getting positions	
	vector<double> getvel(int);							// for getting velocities
	double getpot(int);									// for getting potential energy
	double getkin(int);									// for getting kinetic energy

	// functions on particles
	double distance(particle*,particle*);						// distance between particles
	double distancesquare(particle*,particle*);					// distance square between particles
	vector<double> direction(particle*,particle*);				// calculates direction
	void force(particle*,particle*);							// force between particles
	void potential(particle*,particle*);						// calculating potential energy

	// function on object
	void forcecloud();									// total force on all particls
	void velverletcloud(double);							// verlet velocity step
	void posverletcloud(double);							// verlet position step
	void langevincloud(double);							// langevin dynamics step
	int metropolis(double);
	void speedcloud();									// calculating speeds
	void speedsquarecloud();								// caclulating speed squared
	void potcloud();										// potential energy
	void kincloud();										// kinetic energy
	void tempcalc();											// calculates temperature
	void heatcalc();
	void mcheatcalc();
	void zeropotforce();

	// data collection
	void initPositions();									// initial positions
	void initVelocities();								// initial velocities
	void writedata(vector <double>*,string);
	void printPositions();								// print positions
	void printDistances();								// print distance between all pairs
	void printSpeeds();									// print speeds
	void display(int);									// print stuff to terminal
};


// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
//
// Accessors

void cloud::setpos(int i, double k, double l, double m)		// for setting positons
{
	particleArray[i].pox = k;
	particleArray[i].poy = l;
	particleArray[i].poz = m;
}

void cloud::setvel(int i, double k, double l, double m)		// for setting velocities
{
	particleArray[i].vex = k;
	particleArray[i].vey = l;
	particleArray[i].vez = m;
}

void cloud::setforce(int i, double k, double l, double m)	// for setting forces
{
	particleArray[i].fox = k;
	particleArray[i].foy = l;
	particleArray[i].foz = m;
}

void cloud::setpot(int i, double k)						// for setting potential energy
{
	particleArray[i].pot = k;
}

vector<double> cloud::getpos(int i)						// for getting positions
{
	vector<double> positions = {particleArray[i].pox, particleArray[i].poy, particleArray[i].poz};

	return positions;
}

vector<double> cloud::getvel(int i)						// for getting velocities
{
	vector<double> velocities = {particleArray[i].vex, particleArray[i].vey, particleArray[i].vez};

	return velocities;
}

double cloud::getkin(int i)								// for getting kinetic energy	
{
	return particleArray[i].kin;
}

double cloud::getpot(int i)								// for getting potenial energy
{
	return particleArray[i].pot;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
// Initialization


void cloud::initPositions()								// setting initial positions with face centered packing
{		
	int k = 0;
	for (int kz = 0; kz < natomsside; kz++)
		for (int ky = 0; ky < natomsside; ky++)
			for (int kx = 0; kx < natomsside; kx++, k += 4)		// placing four particles a time
			{   
				setpos(k,kx*gridsize, ky*gridsize, kz*gridsize);
				setpos(k + 1, (kx + 0.5)*gridsize, (ky + 0.5)*gridsize, kz*gridsize);
				setpos(k + 2, (kx + 0.5)*gridsize, ky*gridsize, (kz + 0.5)*gridsize);
				setpos(k + 3, kx*gridsize, (ky + 0.5)*gridsize, (kz + 0.5)*gridsize);
			}   
}

void cloud::initVelocities()							// setting initial velocities, as gaussian in the components
{
	double mean = 0;											// mean velocity
	double std = sqrt(kB*temp/argonmass);							// standard deviation

	default_random_engine generator;							// random number generator to be used
	normal_distribution<double> veldist(mean,std);				// probability distribution to be used

	for (int k = 0; k < natoms; k++)							// setting velocites
	{
		setvel(k,veldist(generator), veldist(generator), veldist(generator));
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
// Functions on the particles

vector<double> cloud::direction (particle* p1, particle* p2)	// calculates the distance vector between two particles
{
	vector<double> direction;
	
	double x,y,z;
	
	x = p1->pox - p2->pox;	
	y = p1->poy - p2->poy;
	z = p1->poz - p2->poz;
		
	x -= static_cast<int> (x >= 0 ? x*resboxlength + 0.5 : x*resboxlength - 0.5)*boxlength;		// periodic boundary conditions
	y -= static_cast<int> (y >= 0 ? y*resboxlength + 0.5 : y*resboxlength - 0.5)*boxlength;
	z -= static_cast<int> (z >= 0 ? z*resboxlength + 0.5 : z*resboxlength - 0.5)*boxlength;

	//cout << "dist1 " << dx << " " << dy << " " << dz << "\n";

	direction = {x,y,z};

	//cout << "dist2 " << direction[0] << " " << direction[1] << " " << direction[2]	<< "\n";

	return direction;
}	

double cloud::distance (particle* p1, particle* p2)				// distance between particles
{
	double r, dx, dy, dz;
	
	vector<double> dir = direction(p1,p2);

	dx = dir[0];
	dy = dir[1];
	dz = dir[2];
	
	r = sqrt(dx*dx + dy*dy + dz*dz);

	return r;
}

double cloud::distancesquare (particle* p1, particle* p2)		// calculating square of distance to save computer time
{
	double r, dx, dy, dz;
	
	vector<double> dir = direction(p1,p2);

	dx = dir[0];
	dy = dir[1];
	dz = dir[2];
	
	r = dx*dx + dy*dy + dz*dz; 

	return r;
}

void cloud::particle::speed (particle* p1)						// calculating speeds
{
	double s = sqrt(pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2));

	 p1 -> spd = s;
}

void cloud::particle::speedsquare (particle* p1)				// caclulating square of speed to save computer time
{
	double s = pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2);
	
	p1 -> spdsqr = s;
}

void cloud::force (particle* p1, particle* p2)					// force between two particles
{
	vector<double> dir = direction(p1,p2);

	double r = distancesquare(p1, p2);
	double F = 4*epsilon*(-6*pow(sigma,6)*pow(r,-4) + 12*pow(sigma,12)*pow(r,-7) ) + potjump;		// Lennard-Jones Force with squared distances to save compter time
	
	double dx,dy,dz;
	
	dx = dir[0];
	dy = dir[1];
	dz = dir[2];

	p1->fox += dx*F;											// calculating forces on p1 additive
	p1->foy += dy*F;	
	p1->foz += dz*F;

	p2->fox -= dx*F;											// calculating forces on p2 additive
	p2->foy -= dy*F;
	p2->foz -= dz*F;
}

void cloud::particle::posverlet(double timestep,particle* p1)	// taking a verlet position step
{
	p1->pox += timestep*p1->vex;
	p1->poy += timestep*p1->vey;
	p1->poz += timestep*p1->vez;
}

void cloud::particle::velverlet(double timestep, particle* p1)	// taking a verlet velocity step, also calculating accelerations
{	
	p1->vex += timestep*p1->fox*resargonmass*0.5;
	p1->vey += timestep*p1->foy*resargonmass*0.5;
	p1->vez += timestep*p1->foz*resargonmass*0.5;
}

void cloud::particle::langevin(double langfric, double langacc, double timestep, particle* p1)	// one langevin time step
{
	p1->vex += timestep*p1->fox*resargonmass*0.5 - p1->vex*langfric + langacc;
	p1->vey += timestep*p1->foy*resargonmass*0.5 - p1->vey*langfric + langacc;
	p1->vez += timestep*p1->foz*resargonmass*0.5 - p1->vez*langfric + langacc;
}


void cloud::particle::kinetic(particle* p1)						// calculates the kinetic energy of a particle
{
	p1->kin = 0.5*argonmass*p1->spdsqr;
}

void cloud::potential(particle* p1, particle* p2)				// calculates the pontential energy between two particles
{
	double r = distancesquare(p1, p2);
	double pot = 4*epsilon*(pow(sigma,12)*pow(r,-6)-pow(sigma,6)*pow(r,-3));		// LJ-potential with squared distance to save CPU time
	
	p1->pot += pot;
	p2->pot += pot;
}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//
// Functions on the cloud

void cloud::speedcloud()									// calculates all speeds
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speed(&particleArray[k]);
}

void cloud::speedsquarecloud()							// calculates all speed squares
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speedsquare(&particleArray[k]);
}

void cloud::forcecloud()								// calculates all forcews
{
	for(int k = 0; k < natoms; k++)								// reset forces
	{
		setforce(k,0,0,0);
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)						// loop over all pairs
	{
		for(int p2 = p1 + 1; p2 < natoms; p2++)

		{
			if (distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)		// distances not computed for r < rc
			{
				force(&particleArray[p1],&particleArray[p2]);
			}
		}
	}
}

void cloud::kincloud()									// calculating kinetic energies
{
	speedsquarecloud();

	for(int k = 0; k < natoms; k++)
	{
		particleArray[k].kinetic(&particleArray[k]);
	}
	
	totkin = 0;
	for(int k = 0; k < natoms; k++)								// calculates total kinetic energy
	{
		totkin += getkin(k);
	}
	
	//cout << totkin << endl;
}

void cloud::potcloud()									// calculating potenital energies
{
	for(int k = 0; k < natoms; k++)								// resetting
	{
		setpot(k,0);
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)						// loop over all pairs
	{
		for (int p2 = p1 + 1; p2 < natoms; p2++)
		{   
			if (distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)		// calculates E if r < rc
			{
				potential(&particleArray[p1],&particleArray[p2]);
			}
		} 
	}

	totpot = 0;
	for(int k = 0; k < natoms; k++)								// calculate total potential energy
	{
		totpot += getpot(k);
	}

	//cout << c->totpot << endl;
}

void cloud::tempcalc()										// caclulate temperature
{
	// which of the two ways is best?

	temp = 2*totkin/(3*natoms*kB);
	
	//double velsqr = 0;

	//for (int k = 0; k < natoms; k++)
	//{
	//	velsqr += c->particleArray[k].spdsqr;
	//}
	
	//c->T = argonmass/(3*natoms*kB)*velsqr;
}

void cloud::heatcalc()
{
	double accensq = 0;

	for (auto& el : particleArray)
	{
		accensq += pow(el.pot+el.kin,2);
	}
	
	double avgensq = accensq*resnatoms;
	double sqavgen = pow((totkin+totpot)*resnatoms,2);
	
	heatcapacity = (avgensq - sqavgen)/(kB*temp*temp);

	//cout << avgensq << " " << sqavgen << " " << temp << " " << heatcapacity << endl;
}

void cloud::mcheatcalc()
{
	double accensq = 0;
	
	for (auto& el : particleArray)
	{
		accensq += pow(el.pot,2);
	}
	
	double avgensq = accensq*resnatoms;
	double sqavgen = pow(totpot*resnatoms,2);
	
	heatcapacity = (avgensq - sqavgen)/(kB*temp*temp) + 1.5*kB;

	//cout << avgensq << " " << sqavgen << " " << temp << " " << heatcapacity << endl;
}

void cloud::posverletcloud(double timestep)				// new positions using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].posverlet(timestep, &particleArray[k]);
}

void cloud::velverletcloud(double timestep)				// new velocities using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].velverlet(timestep, &particleArray[k]);
}

void cloud::langevincloud(double timestep)				// langevin dynamics
{
	double mean = 0;											// mean velocity
	double std = sqrt(2*langfric*argonmass*kB*temp);				// standard deviation
	
	default_random_engine generator;							// random number generator to be used
	normal_distribution<double> langacc(mean,std);				// gaussian
	
	for(int k = 0; k < natoms; k++)								// apply dynamics
		particleArray[k].langevin(langfric, langacc(generator), timestep, &particleArray[k]);
}

int cloud::metropolis(double maxdisp)
{

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-maxdisp,maxdisp);

	uniform_int_distribution<> randpart(0,natoms-1);

	int k = randpart(gen);

	vector <double> old = getpos(k);
	
	double potold = totpot;
		
	double dx = dis(gen);
	double dy = dis(gen);
	double dz = dis(gen);

	setpos(k,old[0]+dx,old[1]+dy,old[2]+dz);

	potcloud();

	double diff = exp(-(totpot - potold)/(kB*temp*temp));
	bernoulli_distribution accprob(diff);

	if (totpot <= potold)
		return 1;
	else if(accprob(gen))
		return 1;
	else
		{setpos(k,old[0],old[1],old[2]); return 0;}
}

void cloud::zeropotforce()
{
	for (auto& el : particleArray)
	{
		el.pot = 0;
		el.fox = 0;
		el.foy = 0;
		el.foz = 0;
	}
	totkin = 0;
	totpot = 0;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//
// Printing functions

void cloud::writedata(vector <double>* data, string filename)
{
	ofstream writestream;
	writestream.open(filename);

	for (auto& el : *data)
	   writestream << el << endl;

	writestream.close();
}

void cloud::printPositions()									// printing positions
{
	ofstream filestream;
	filestream.open("positions.txt");
	for(int k = 0; k < natoms; k++)
		filestream << particleArray[k].pox << " "  << particleArray[k].poy << " " << particleArray[k].poz << "\n";
	filestream.close();
}

void cloud::printSpeeds()										// printing speeds
{
	ofstream filestream;
	filestream.open("speeds.txt");
	double vel;
	for(int k = 0; k < natoms; k++)
	{
		vel = particleArray[k].spd;
		filestream << vel << "\n";
	}
	filestream.close();
}

void cloud::printDistances()									// print distance of all pairs
{
	ofstream filestream;
	filestream.open("distances.txt");
	
	for(int k = 0; k < natoms - 1; k++)
	{   
		for(int l = k + 1; l < natoms; l++)
		{   
			filestream << distance(&particleArray[k],&particleArray[l]) << endl;		// build a vector of all distances
		}
	}

	filestream.close();
}



void cloud::display(int k)									// for printing to screen
{
	//cout << "pos: " << particleArray[k].pox << " "  << particleArray[k].poy << " " << particleArray[k].poz << "\n"; 
	cout << "vel: " << particleArray[k].vex << " "  << particleArray[k].vey << " " << particleArray[k].vez << "\n"; 
	//cout << "force: " << particleArray[k].fox << " "  << particleArray[k].foy << " " << particleArray[k].foz << "\n"; 
	//cout << "Ekin: " << particleArray[k].kin << "\n"; 
	//cout << "Epot: " << particleArray[k].pot << "\n";
	//cout << "Epot: " << particleArray[k].pot << "\n";
}

// Non object member functions

// Warning: really, really ugly code ahead, this will all look better after christmas when I get all data analysis to be in the same program. Did a bad move in my implementation earlier bevause I didn't plan properly for error analysis and this is the quickfix...

vector <double>  uglyerror(cloud* c)												
{
	double avgkin = c->totkin/c->natoms;

	vector <double> shifted (c->natoms);

	for (auto& el: shifted)
	{
		auto k = &el - &shifted[0];
		el = pow(c->particleArray[k].kin - avgkin,2);

	}

	double kinvariance = accumulate(shifted.begin(),shifted.end(),0.0)/shifted.size();
	double kinstdev = sqrt(kinvariance);
	double kinerror = kinstdev/sqrt(c->natoms-1);

	//double testerror = sqrt(1.0/(c->natoms-1.0)*accumulate(shifted.begin(),shifted.end(),0.0));

	double avgpot = c->totpot/c->natoms;

	for (auto& el: shifted)
	{
		auto k = &el - &shifted[0];
		el = pow(c->particleArray[k].pot - avgpot,2);

	}

	double potvariance = accumulate(shifted.begin(),shifted.end(),0.0)/shifted.size();
	double potstdev = sqrt(potvariance);
	double poterror = potstdev/sqrt(c->natoms-1);
	
	double avgen = (c->totkin+c->totpot)/c->natoms;
	
	for (auto& el: shifted)
	{
		auto k = &el - &shifted[0];
		el = pow(c->particleArray[k].kin + c->particleArray[k].pot - avgen,2);

	}

	double envariance = accumulate(shifted.begin(),shifted.end(),0.0)/shifted.size();
	double enstdev = sqrt(envariance);
	double enerror = enstdev/sqrt(c->natoms-1);

	double temperror = 2*kinerror*c->natoms/(3*c->natoms*kB);

	double accensq = 0;

	for (auto& el : c->particleArray)
	{
		accensq += pow(el.pot + el.kin + enerror,2);
	}
	
	double avgensq = accensq*c->resnatoms;
	double sqavgen = pow((c->totkin + c->totpot + c->natoms*enerror)*c->resnatoms,2);
	
	double shiftheat = (avgensq - sqavgen)/(kB*pow(c->temp+temperror,2));
	
	double accmcensq = 0;
	
	for (auto& el : c->particleArray)
	{
		accmcensq += pow(el.pot + poterror,2);
	}
	
	double avgmcensq = accmcensq*c->resnatoms;
	double sqavgmcen = pow((c->totpot + c->natoms*poterror)*c->resnatoms,2);
	
	double shiftheatmc = (avgmcensq - sqavgmcen)/(kB*pow(c->temp+temperror,2));

	vector <double> retvec = {kinerror, kinvariance, poterror, potvariance, enerror, envariance, temperror, shiftheat, shiftheatmc};

	return retvec;
}

\

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------

// Main program

int main()
{

	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();		//for measuring execution time
	
	int side = 6;													// box side in number of atoms
	int atoms = 4*pow(side,3);										// total number of atoms
	double length = 10.229*sigma;									// box length
	double grid = length/side;										// initial separation of atoms
	double langfric = 10e-6;
	double temperature = 94.4;

	// creating particle cloud
	cloud particleCloud(side, atoms, length, grid, langfric, temperature);
	
	// set initial conditions
	srand(time(NULL));												// using current time as random seed
	particleCloud.initPositions();					// initializing positions
	particleCloud.initVelocities();					// initializing velocities
	particleCloud.forcecloud();						// calculating initial force
	particleCloud.potcloud();							// calculating potential energies

	// simulation specifics
	double timestep = 1e-14;										// time step to be used
	int duration = 2000;											// duration of simulation
	double maxdisp = 2.25e-1*grid;
	double count = 0;
	
	ofstream fskin, fspot, fsetot, fsheat, fsmcheat, fstemp, fsparticle, fserror;
	fskin.open("kinetic.txt"), fspot.open("potential.txt"), fsetot.open("totalenergy.txt"), fstemp.open("temperature.txt"), fsheat.open("heatcapacity.txt"), fsmcheat.open("mcheatcap.txt"), fsparticle.open("particle.txt"), fserror.open("error1.txt");
	int particleid = 0.5*864;

	for(int k = 0; k< duration; k++)								// main simulation loop
	{
		//particleCloud.velverletcloud(timestep);	// calculating new velocities
		//particleCloud.langevincloud(timestep);
		//particleCloud.posverletcloud(timestep);		// caluclating new positions
		//particleCloud.forcecloud();					// calculating new forces
		//particleCloud.velverletcloud(timestep);	// calculating new velocities
		//particleCloud.langevincloud(timestep);
		
		int comp = particleCloud.metropolis(maxdisp);
		
		if (comp == 1)
			count++;

		//particleCloud.kincloud();						// calculating kinetic energies
		//particleCloud.potcloud();						// calculating potential energies
		//particleCloud.tempcalc();
		particleCloud.heatcalc();
		//particleCloud.mcheatcalc();
		vector <double> errors = uglyerror(&particleCloud);
		
		errors[7] = -errors[7] + particleCloud.heatcapacity;
		errors[8] = errors[8] - particleCloud.heatcapacity;

		for (auto& el : errors)
			fserror << el << " " ;
		fserror << endl;
		
		fsparticle << particleCloud.getvel(particleid)[1] << endl;
		fskin << particleCloud.totkin/atoms << endl;
		fspot << particleCloud.totpot/atoms << endl;
		fsetot << particleCloud.totkin + particleCloud.totpot << endl;
		fstemp << particleCloud.temp << endl;
		fsheat << particleCloud.heatcapacity << endl;
	}

	fskin.close(), fspot.close(), fsetot.close(), fstemp.close(), fsheat.close(), fsmcheat.close(), fsparticle.close(), fserror.close();	

	particleCloud.speedcloud();						// calculating speeds
	particleCloud.speedsquarecloud();					// calculating square speeds
	
	particleCloud.printPositions();					// print positions
	particleCloud.printDistances();
	particleCloud.printSpeeds();						// print raw speeds

	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	//ofstream noice;
	//noice.open("particle.txt");
	//double rval = 1;
	//for ( int k = 1; k < 100; k++)
	//{
	//	//rval = rval*0.95 + 2*((double) rand() / (RAND_MAX)) - 1;
	//	rval -= 0.00001;
	//	noice << rval << endl;
	//}
	//noice.close();		

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	double accratio = count / duration;
	
	cout << accratio << endl;
	cout << benchtime << "\n" ;

	//while(true) // interface for simple debugging
	//{
	//	int k;
	//	cout << "give k\n";
	//	cin >> k;
	//	particleCloud.display(k,&particleCloud);
	//}


//	cloud testcloud(side, 2, length, grid, langfric, temperature);
//
//	testcloud.setpos(0,0,0,0*sigma);
//	testcloud.setpos(1,0,0,3*sigma);
//	testcloud.setvel(0,0,0,100);
//	testcloud.setvel(1,0,0,-100);
//
//	duration = 5000;
//	timestep = 1e-15;
//
//	ofstream test; test.open("test.txt");
//
//	for(int k = 0; k < duration; k++)
//	{
//		testcloud.velverletcloud(timestep);	// calculating new velocities
//		testcloud.posverletcloud(timestep);		// caluclating new positions
//		testcloud.forcecloud();					// calculating new forces
//		testcloud.velverletcloud(timestep);	// calculating new velocities
//		
//		testcloud.kincloud();						// calculating kinetic energies
//		testcloud.potcloud();						// calculating potential energies
//
//		test << testcloud.totkin << endl;
//		test << testcloud.totpot << endl;
//		test << testcloud.totkin + testcloud.totpot << endl;
//	}
//
//	test.close();

	return 0;
}









