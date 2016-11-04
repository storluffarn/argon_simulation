// 
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

#include <iostream>												// basic io stuff
#include <fstream>												// writing to file
#include <cmath>												// for math functions
#include <random>												// for random generation
#include <algorithm>											// for count_if
#include <vector>												// for using the vector type members
#include <chrono>												// for measuring execution times
#include <ctime>												// for random seed time

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
	friend cloud;											// a class of particles

	// public members
	public:
	
	// public function members
	void posverlet(double,particle*);							// verlet velocity function
	void velverlet(double,particle*);							// vertet position function
	void langevin(double,double,double,particle*);
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
	double gridsize;											// initial separation of particles
	double volume;
	double density;
	double partdens;
	double langfric;
	
	// object member
	vector<particle> particleArray{static_cast<unsigned int> (natoms)};		// vector of all particles

	// function members
	
	double calcpairs()									// calculate the number of pairs
	{
		int sum = 0.5*(natoms - 1)*(natoms);						// arithemtic sum
		return sum;
	}

	double calcvolume()
	{
		return pow(boxlength,3);
	}

	double calcdensity()
	{
		return natoms * argonmass / volume;
	}

	double calcpartdens()
	{
		return natoms / volume;
	}
	
	void calcavgkin()
	{
		avgkin = totkin/natoms;
	}

	void calcavgpot()
	{
		avgpot = totpot/natoms;
	}
	
// public section
	
	public:
	
	// data members

	// cloud data
	double T = 94.4;											// temperature
	double totkin;												// total kinetic energy
	double totpot;												// total potential energy
	double avgkin;
	double avgpot;

	// constructor
	
	cloud(int t1, unsigned int t2, double t3, double t4, double t5)			// setting up cloud
		: natomsside(t1), natoms(t2), boxlength(t3), gridsize(t4), langfric(t5)
	{
		resboxlength = 1 / boxlength; 
		pairs = calcpairs();
		volume = calcvolume();
		density = calcdensity();
		partdens = calcpartdens();
	}	
	
	// function members

	// accessors
	void setpos(int,double,double,double,cloud*);				// for setting positions
	void setvel(int,double,double,double,cloud*);				// for setting velcities
	void setforce(int,double,double,double,cloud*);				// for setting forces
	void setpot(int,double,cloud*);								// for setting potential energy
	vector<double> getpos(int,cloud*);							// for getting positions	
	vector<double> getvel(int,cloud*);							// for getting velocities
	double getpot(int,cloud*);									// for getting potential energy
	double getkin(int,cloud*);									// for getting kinetic energy

	// functions on particles
	double distance(particle*,particle*);						// distance between particles
	double distancesquare(particle*,particle*);					// distance between particles
	vector<double> direction(particle*,particle*);				// calculates direction
	void force(particle*,particle*);							// force between particles
	void potential(particle*,particle*);						// calculating potential energy

	// function on object
	void forcecloud(cloud*);									// total force on all particls
	void velverletcloud(double,cloud*);							// verlet velocity step
	void posverletcloud(double,cloud*);							// verlet position step
	void langevincloud(double,cloud*);
	void speedcloud(cloud*);									// calculating speeds
	void speedsquarecloud(cloud*);								// caclulating speed squared
	void potcloud(cloud*);										// potential energy
	void kincloud(cloud*);										// kinetic energy
	void temp(cloud*);
	void radialdistfunc(int,vector<double>*,cloud*);			// print density over distanc
	void veldist(double,int,vector<double>*,cloud*);			// print velocity distribution

	// data collection
	void initPositions(cloud*);									// initial positions
	void initVelocities(cloud*);								// initial velocities
	void printPositions(cloud*);								// print positions
	void printSpeeds(cloud*);									// print speeds
	void printDensArr(int,vector<double>*);
	void printSpeedDist(int,vector<double>*);
	void display(int,cloud*);									// print stuff to terminal
};


// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
//
// Accessors

void cloud::setpos(int i, double k, double l, double m, cloud*)			// for setting positons
{
	particleArray[i].pox = k;
	particleArray[i].poy = l;
	particleArray[i].poz = m;
}

void cloud::setvel(int i, double k, double l, double m, cloud*)			// for setting velocities
{
	particleArray[i].vex = k;
	particleArray[i].vey = l;
	particleArray[i].vez = m;
}

void cloud::setforce(int i, double k, double l, double m, cloud*)		// for setting forces
{
	particleArray[i].fox = k;
	particleArray[i].foy = l;
	particleArray[i].foz = m;
}

void cloud::setpot(int i, double k, cloud*)						// for setting potential energy
{
	particleArray[i].pot = k;
}

vector<double> cloud::getpos(int i, cloud*)						// for getting positions
{
	vector<double> positions = {particleArray[i].pox, particleArray[i].poy, particleArray[i].poz};

	return positions;
}

vector<double> cloud::getvel(int i, cloud*)						// for getting velocities
{
	vector<double> velocities = {particleArray[i].vex, particleArray[i].vey, particleArray[i].vez};

	return velocities;
}

double cloud::getkin(int i, cloud*)								// for getting kinetic energy	
{
	return particleArray[i].kin;
}

double cloud::getpot(int i, cloud*)								// for getting potenial energy
{
	return particleArray[i].pot;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
// Initialization


void cloud::initPositions(cloud* c)								// setting initial positions with face centered packing
{		
	int k = 0;
	for (int kz = 0; kz < natomsside; kz++)
		for (int ky = 0; ky < natomsside; ky++)
			for (int kx = 0; kx < natomsside; kx++, k += 4)		// placing four particles a time
			{   
				setpos(k,kx*gridsize, ky*gridsize, kz*gridsize, c);
				setpos(k + 1, (kx + 0.5)*gridsize, (ky + 0.5)*gridsize, kz*gridsize, c);
				setpos(k + 2, (kx + 0.5)*gridsize, ky*gridsize, (kz + 0.5)*gridsize, c);
				setpos(k + 3, kx*gridsize, (ky + 0.5)*gridsize, (kz + 0.5)*gridsize, c);
			}   
}

void cloud::initVelocities(cloud* c)							// setting initial velocities, as gaussian in the components
{
	double mean = 0;											// mean velocity
	double std = sqrt(kB*T/argonmass);							// standard deviation

	default_random_engine generator;							// random number generator to be used
	normal_distribution<double> veldist(mean,std);				// probability distribution to be used

	for (int k = 0; k < natoms; k++)							// setting velocites
	{
		setvel(k,veldist(generator), veldist(generator), veldist(generator), c);
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
		
	x -= static_cast<int> (x >= 0 ? z*resboxlength + 0.5 : x*resboxlength - 0.5)*boxlength;		// periodic boundary conditions
	y -= static_cast<int> (y >= 0 ? z*resboxlength + 0.5 : y*resboxlength - 0.5)*boxlength;
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

void cloud::particle::speed (particle* p1)								// calculating speeds
{
	double s = sqrt(pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2));

	 p1 -> spd = s;
}

void cloud::particle::speedsquare (particle* p1)						// caclulating square of speed to save computer time
{
	double s = pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2);
	
	p1 -> spdsqr = s;
}

void cloud::force (particle* p1, particle* p2)					// force between two particles
{
	vector<double> dir = direction(p1,p2);

	double r = distancesquare(p1, p2);
	double F = 4*epsilon*(-6*pow(sigma,6)*pow(r,-4) + 12*pow(sigma,12)*pow(r,-7) );		// Lennard-Jones Force with squared distances to save compter time
	
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

void cloud::particle::posverlet(double timestep,particle* p1)			// taking a verlet position step
{
	p1->pox += timestep*p1->vex;
	p1->poy += timestep*p1->vey;
	p1->poz += timestep*p1->vez;
}

void cloud::particle::velverlet(double timestep, particle* p1)			// taking a verlet velocity step, also calculating accelerations
{	
	p1->vex += timestep*p1->fox*resargonmass*0.5;
	p1->vey += timestep*p1->foy*resargonmass*0.5;
	p1->vez += timestep*p1->foz*resargonmass*0.5;
}

void cloud::particle::langevin(double langfric, double langacc, double timestep, particle* p1)
{
	p1->vex += timestep*p1->fox*resargonmass*0.5 - langfric + langacc;
	p1->vey += timestep*p1->foy*resargonmass*0.5 - langfric + langacc;
	p1->vez += timestep*p1->foz*resargonmass*0.5 - langfric + langacc;
}


void cloud::particle::kinetic(particle* p1)							// calculates the kinetic energy of a particle
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

void cloud::speedcloud(cloud*)									// calculates all speeds
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speed(&particleArray[k]);
}

void cloud::speedsquarecloud(cloud*)							// calculates all speed squares
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speedsquare(&particleArray[k]);
}

void cloud::forcecloud(cloud* c)								// calculates all forcews
{
	for(int k = 0; k < natoms; k++)								// resets forces
	{
		setforce(k,0,0,0,c);
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)						// calculates the force on two particles
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

void cloud::kincloud(cloud* c)								// calculating kinetic energies
{
	c->speedsquarecloud(c);

	for(int k = 0; k < natoms; k++)
	{
		particleArray[k].kinetic(&particleArray[k]);
	}
	
	c->totkin = 0;
	for(int k = 0; k < natoms; k++)								// a bit cumbersome but clear?
	{
		c->totkin += getkin(k ,c);
	}
	
	//cout << totkin << endl;
}

void cloud::potcloud(cloud* c)								// calculating potenital energies
{
	for(int k = 0; k < natoms; k++)								// resetting
	{
		setpot(k,0,c);
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)
	{
		for (int p2 = p1 + 1; p2 < natoms; p2++)
		{   
			if (distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)		// calculates E if r < rc
			{
				potential(&particleArray[p1],&particleArray[p2]);
			}
		} 
	}

	c->totpot = 0;
	for(int k = 0; k < natoms; k++)
	{
		c->totpot += getpot(k,c);
	}
}

void cloud::temp(cloud* c)
{
	c->T = 2*totkin/(3*natoms*kB);
}

void cloud::posverletcloud(double timestep, cloud*)				// new positions using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].posverlet(timestep, &particleArray[k]);
}

void cloud::velverletcloud(double timestep, cloud*)				// new velocities using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].velverlet(timestep, &particleArray[k]);
}

void cloud::langevincloud(double timestep, cloud*)			
{
	double mean = 0;											// mean velocity
	double std = sqrt(2*langfric*argonmass*kB*T);				// standard deviation
	
	default_random_engine generator;							// random number generator to be used
	normal_distribution<double> langacc(mean,std);				// probability distribution to be used
	
	for(int k = 0; k < natoms; k++)
		particleArray[k].langevin(langfric, langacc(generator), timestep, &particleArray[k]);
}

void cloud::radialdistfunc(int densbins, vector<double> *densArr, cloud*)		// calculates the radial distribution
{
	vector<double> distances (pairs);											// length is number of pairs, sum (i=1->n-1) (i)
	double resdensbins = 1/(static_cast<double> (densbins));					// reciprocal densbins
	
	unsigned int m = 0;
	for(int k = 0; k < natoms - 1; k++)
	{   
		for(int l = k + 1; l < natoms; l++)
		{   
			distances[m] = distance(&particleArray[k],&particleArray[l]);		// build a vector of all distances
			m++;
		}
	}

	vector<int>hist(densbins);
	double shell1 = 0, shell2;
	
	for(int k = 0; k < densbins; k++)
	{   
		shell2 = ((double) k + 1)*resdensbins*boxlength*1.5;
		hist[k] = std::count_if(distances.begin(), distances.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;}); // count distances into bins
		shell1 = shell2;
	}

	shell1 = 0, shell2 = 0;
	double vol;
	for (double k = 0; k < densbins; k++)
	{   
		shell2 = (4/3*pi*pow((k + 1)*resdensbins*boxlength,3));		// calculade densities
		//shell2 = (k + 1)*resdensbins*boxlength;		// calculade densities
		vol = shell2 - shell1;
		//vol = 4*pi*shell2*shell2;
		(*densArr)[k] = hist[k]/(partdens*vol);
		shell1 = shell2;
	}
}

void cloud::veldist(double maxspeed, int velbins, vector<double> *velArr, cloud*)		// calculates the velocity distrobution
{
	vector<double> velocities(natoms);
	double resvelbins = 1/(static_cast<double> (velbins));						// reciprocal velbins

	for(int k = 0; k < natoms; k++)
	{   
		velocities[k] = particleArray[k].spd;
	}

	double shell1 = 0, shell2;
	for(double k = 0; k < velbins; k++)
	{   
		shell2 = (k + 1)*resvelbins*maxspeed;
		(*velArr)[k] = std::count_if(velocities.begin(), velocities.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;}); // count velocities into bins
		shell1 = shell2;
	}
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//
// Printing functions

void cloud::printPositions(cloud*)									// printing positions
{
	ofstream filestream;
	filestream.open("positions.txt");
	for(int k = 0; k < natoms; k++)
		filestream << particleArray[k].pox << " "  << particleArray[k].poy << " " << particleArray[k].poz << "\n";
	filestream.close();
}

void cloud::printSpeeds(cloud*)										// printing speeds
{
	ofstream filestream;
	filestream.open("velocities.txt");
	double vel;
	for(int k = 0; k < natoms; k++)
	{
		vel = particleArray[k].spd;
		filestream << vel << "\n";
	}
	filestream.close();
}

void cloud::printDensArr(int densbins, vector<double> *densArr)		// print velocity distribution
{
	ofstream filestream;
	filestream.open("radialdistfunc.txt");
	for(int k = 0; k < densbins; k++)
		filestream << (*densArr)[k]  << "\n";
	filestream.close();
}

void cloud::printSpeedDist(int velbins, vector<double> *velArr)		// print velocity distribution
{
	ofstream filestream;
	filestream.open("speeddist.txt");
	for(int k = 0; k < velbins; k++)
	{
		filestream << (*velArr)[k]  << "\n";
	}
	
	filestream.close();
}

void cloud::display(int k, cloud*)
{
	//cout << "pos: " << particleArray[k].pox << " "  << particleArray[k].poy << " " << particleArray[k].poz << "\n"; 
	cout << "vel: " << particleArray[k].vex << " "  << particleArray[k].vey << " " << particleArray[k].vez << "\n"; 
	//cout << "force: " << particleArray[k].fox << " "  << particleArray[k].foy << " " << particleArray[k].foz << "\n"; 
	//cout << "Ekin: " << particleArray[k].kin << "\n"; 
	//cout << "Epot: " << particleArray[k].pot << "\n";
	//cout << "Epot: " << particleArray[k].pot << "\n";
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
// Main program



int main()
{
	// simulation specifics
	double timestep = 1e-14;										// time step to be used
	int duration = 20;											// duration of simulation

	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();		//for measuring execution time
	
	int side = 6;													// box side in number of atoms
	int atoms = 4*pow(side,3);										// total number of atoms
	double length = 10.229*sigma;									// box length
	double grid = length/side;										// initial separation of atoms

	// creating particle cloud
	cloud particleCloud(side, atoms, length, grid, 10e-5);
	
	// set initial conditions
	srand(time(NULL));												// using current time as random seed
	particleCloud.initPositions(&particleCloud);					// initializing positions
	particleCloud.initVelocities(&particleCloud);					// initializing velocities
	particleCloud.forcecloud(&particleCloud);						// calculating initial forces
	
	ofstream filestream;
	filestream.open("loopoutput.txt");

	for(int k = 0; k< duration; k++)								// main simulation loop
	{
		//particleCloud.velverletcloud(timestep, &particleCloud);		// calculating new velocities
		particleCloud.langevincloud(timestep, &particleCloud);
		particleCloud.posverletcloud(timestep, &particleCloud);		// caluclating new positions
		particleCloud.forcecloud(&particleCloud);					// calculating new forces
		//particleCloud.velverletcloud(timestep, &particleCloud);		// calculating new velocities
		particleCloud.langevincloud(timestep, &particleCloud);

		particleCloud.kincloud(&particleCloud);			// calculating kinetic energies
		particleCloud.potcloud(&particleCloud);			// calculating potential energies
		particleCloud.temp(&particleCloud);

		filestream << particleCloud.totkin << " " << particleCloud.totpot << " " << particleCloud.T << endl;
	}

	filestream.close();	

	particleCloud.speedcloud(&particleCloud);						// calculating speeds
	particleCloud.speedsquarecloud(&particleCloud);					// calculating square speeds

	int velbins = 100;												// number of bins in speed distribution
	double maxspeed = 800.0;										// maximum speed
	vector<double> velArr (velbins);								// vector of all speeds
	particleCloud.veldist(maxspeed, velbins, &velArr, &particleCloud);		// calculating speed distribution
	
	int densbins = 500;												// number of bins in density distribution
	vector<double> densArr (densbins);								// vector of density over distance
	particleCloud.radialdistfunc(densbins, &densArr, &particleCloud); // calculating radial distribution
	
	particleCloud.printPositions(&particleCloud);					// print positions
	particleCloud.printSpeeds(&particleCloud);						// print raw speeds
	particleCloud.printSpeedDist(velbins, &velArr);					// print speed distribution
	particleCloud.printDensArr(densbins, &densArr);					// print density distribution

	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	cout << benchtime << "\n" ;

	//while(true) // interface for simple debugging
	//{
	//	int k;
	//	cout << "give k\n";
	//	cin >> k;
	//	particleCloud.display(k,&particleCloud);
	//}

	return 0;
}









