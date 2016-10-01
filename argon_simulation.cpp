
// Molecular dynamics dimulation of argon liquid n = 864, T = 120 K.
// Originally performed by A. Rahman, Phys letters, vol 136:2A

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//
// Preliminaries

// includes

#include <iostream>			// basic io stuffi
#include <fstream>			// writing to file
#include <cmath>			// for math functions
#include <random>			// for random generation
#include <algorithm>		// for count_if
#include <vector>			// for using the vector type members
#include <chrono>			// for measuring execution times

using namespace std;		// tidy up

// constants

const double umass = 1.660539040e-27;				// atomic mass
const double argonmass = umass*39.948;				// argon mass
const int T = 120;									// temperature
const double pi = atan(1)*4;						// for using pi
const double kB = 1.38064852e-23;					// Boltzmann's constant
const double sigma = 3.4e-10;						// interaction term
const double epsilon = 120*kB;						// energy term
const double dcutoff = pow(2,1/6)*sigma;			// interaction cutoff 
const double dcutoffsquare = pow(dcutoff,2);		// interaction cutoff 

// classes

class particle //class to hold particle data and functions
{
	// data members
	double pox;		// positions
	double poy;
	double poz;
	double vex;		// velocities
	double vey;
	double vez;
	double fox;
	double foy;
	double foz;
	double spd;
	double spdsqr;
	double kin;		// potential energy
	double pot;		// kinetic energy

	// member functions (none)

	// freinds with bennifits
	friend class cloud;						// a class of particles

	// 
	public:
	void posverlet(double,particle*);		// verlet velocity function
	void velverlet(double,particle*);		// vertet position function
	void kinetic(particle*);				// calculating kinetic energy
	void speed(particle*);
	void speedsquare(particle*);
	
	// accessors
	void setpos(double,double,double);		// for setting parameters
	void setvel(double,double,double);		
	vector<double> getpos();				// for getting parameters			
	vector<double> getvel();	
};

class cloud // class for the set of all particles and functions on these
{
	int natomsside;
	int natoms;
	double boxlength;			// from paper
	double gridsize;	// initial separation

	unsigned int ugly = natoms;		// why does it need to be unsigned?
	
	// object member
	vector<particle> particleArray{ugly};						// vector of all particles
	
	public:
	cloud(int a, unsigned int b, double c, double d)
		: natomsside(a), natoms(b), boxlength(c), gridsize(d)
	{}
	
	double distance(particle*,particle*);				// distance between particles
	double distancesquare(particle*,particle*);			// distance between particles
	vector<double> direction(particle*,particle*);
	void force(particle*,particle*);						// force between particles
	void potential(particle*,particle*);					// calculating potential energy

	// function on object
	void forcecloud(cloud*);							// total force on all particls
	void velverletcloud(double,cloud*);						// verlet velocity step
	void posverletcloud(double,cloud*);						// verlet position step
	void speedcloud(cloud*);
	void speedsquarecloud(cloud*);
	void potcloud(cloud*);							// potential energy
	void kincloud(cloud*);									// kinetic energy
	void radialdistfunc(int,vector<double>*,cloud*);	// print density over distanc
	void veldist(double,int,vector<double>*,cloud*);		// print velocity distribution

	// data collection
	void initPositions(cloud*);						// initial positions
	void initVelocities(cloud*);							// initial velocities
	void printPositions(cloud*);							// print positions
	void printSpeeds(cloud*);								// print speeds
	void printDensArr(int,vector<double>*);
	void printSpeedDist(int,vector<double>*);
	void display(int,cloud*);									// print stuff to terminal
};


// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
//
// Accessors

void particle::setpos(double k, double l, double m)		// for setting positons
{
	pox = k;
	poy = l;
	poz = m;
}

void particle::setvel(double k, double l, double m)		// for setting velocities
{
	vex = k;
	vey = l;
	vez = m;
}

vector<double> particle::getpos()		// for getting positions
{
	std::vector<double> returnvals;
	returnvals.push_back(pox);
	returnvals.push_back(poy);
	returnvals.push_back(poz);

	return returnvals;
}

vector<double> particle::getvel()		// for getting velocities
{
	std::vector<double> returnvals;
	returnvals.push_back(vex);
	returnvals.push_back(vey);
	returnvals.push_back(vez);

	return returnvals;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
//Initializing


void cloud::initPositions(cloud*) // setting initial positions with face centered packing
{		
	int k = 0;
	for (int kz = 0; kz < natomsside; kz++)
		for (int ky = 0; ky < natomsside; ky++)
			for (int kx = 0; kx < natomsside; kx++, k += 4) // placing four particles a time
			{   
				double displace = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				displace = 1; // uncomment to add white noise

				particleArray[k].setpos(kx*gridsize*displace, ky*gridsize*displace, kz*gridsize*displace);
				particleArray[k+1].setpos((kx + 0.5)*gridsize*displace, (ky + 0.5)*gridsize*displace, kz*gridsize*displace);
				particleArray[k+2].setpos((kx + 0.5)*gridsize*displace, ky*gridsize*displace, (kz + 0.5)*gridsize*displace);
				particleArray[k+3].setpos(kx*gridsize*displace, (ky + 0.5)*gridsize*displace, (kz + 0.5)*gridsize*displace);
			}   
}

void cloud::initVelocities(cloud*) // setting initial velocities, assuming approx gaussian
{
	double mean = 0; // centre of 1d velocities
	double std = sqrt(kB*T/argonmass); // standard deviation

	default_random_engine generator; // random number generator to be used
	normal_distribution<double> veldist(mean,std); // probability distribution

	for (int k = 0; k < natoms; k++) // setting velocites
	{
		particleArray[k].setvel(veldist(generator), veldist(generator), veldist(generator));
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
// Functions on the particles

vector<double> cloud::direction (particle* p1, particle* p2)
{
	vector<double> direction;
	
	double dx,dy,dz;
		
	dx = fabs(p1->pox - p2->pox); // 1d distances
	dy = fabs(p1->poy - p2->poy);
	dz = fabs(p1->poz - p2->poz);

	dx -= static_cast<int> (dx/boxlength + 0.5) * boxlength; //periodic boundary
	dy -= static_cast<int> (dy/boxlength + 0.5) * boxlength;
	dz -= static_cast<int> (dz/boxlength + 0.5) * boxlength;

	direction.push_back(dx);
	direction.push_back(dy);
	direction.push_back(dz);

	return direction;
}	

double cloud::distance (particle* p1, particle* p2) // distance between particles with periodic boundaries
{
	double r, dx, dy, dz;
	
	vector<double> dir = direction(p1,p2);

	dx = dir[0];
	dy = dir[1];
	dz = dir[2];
	
	r = sqrt(dx*dx + dy*dy + dz*dz); // euklidian distance

	return r;
}

double cloud::distancesquare (particle* p1, particle* p2) // calculating square of distance to save computer time
{
	double r, dx, dy, dz;
	
	vector<double> dir = direction(p1,p2);

	dx = dir[0];
	dy = dir[1];
	dz = dir[2];
	
	r = dx*dx + dy*dy + dz*dz; // euklidian distance squared

	return r;
}

void particle::speed (particle* p1) // calculating speeds
{
	double s = sqrt(pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2));

	 p1->spd = s;
}

void particle::speedsquare (particle* p1) // caclulating square of speed without sqrt to save computer time
{
	double s = pow(p1->vex,2)+pow(p1->vey,2)+pow(p1->vez,2);

	p1->spdsqr = s;
}

void cloud::force (particle* p1, particle* p2) // force between two particles
{
	vector<double> dir = direction(p1,p2);

	double r = distancesquare(p1, p2); // see distance function
	double F = 24*epsilon*(2*pow(sigma,12)/pow(r,7)-pow(sigma,6)/pow(r,4)); // LJ-force with halved exponents to save computer time
	//double r = distance(boxlength, p1, p2); // see distance function
	//double F = 24*epsilon*(pow(sigma,6)*pow(r,8) - 2*pow(sigma,12)/pow(r,14) );

	double dx,dy,dz;
	
	dx = dir[0];
	dy = dir[1];
	dz = dir[2];
	
	p1->fox += dx*F;	
	p1->foy += dy*F;	
	p1->foz += dz*F;

	p2->fox -= dx*F;
	p2->foy -= dy*F;
	p2->foz -= dz*F;
}


void particle::posverlet(double timestep,particle* p1) // taking a verlet position step
{
	p1->pox += timestep*p1->vex;
	p1->poy += timestep*p1->vey;
	p1->poz += timestep*p1->vez;
}

void particle::velverlet(double timestep, particle* p1) // taking a verlet velocity step, also calculating accelerations
{	
	double argonmas = 1;

	p1->vex += timestep*p1->fox/argonmas*0.5;
	p1->vey += timestep*p1->foy/argonmas*0.5;
	p1->vez += timestep*p1->foz/argonmas*0.5;
}

void particle::kinetic(particle* p1) // calculates the kinetic energy of a particle
{
	p1->kin = 0.5*argonmass*p1->spdsqr;
}

void cloud::potential(particle* p1, particle* p2)
{
	double r = distancesquare(p1, p2);
	double pot = 4*epsilon*(pow(sigma/r,6)-pow(sigma/r,3)); // LJ-potential with squared ristance

	p1->pot += pot;
	p2->pot += pot;
}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//
// Functions on the cloud

void cloud::speedcloud(cloud*)
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speed(&particleArray[k]);
}

void cloud::speedsquarecloud(cloud*)
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].speedsquare(&particleArray[k]);
}

void cloud::forcecloud(cloud*)
{
	for(int k = 0; k < natoms; k++)
	{
		particleArray[k].fox = 0;
		particleArray[k].foy = 0;
		particleArray[k].foz = 0;
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)
	{
		for(int p2 = p1 + 1; p2 < natoms - 1; p2++)
		{
			if (p1 != p2 && distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)
			{
				force(&particleArray[p1],&particleArray[p2]);
			}
		}	
	}
}

void cloud::kincloud(cloud*) // calculating kinetic energies
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].kinetic(&particleArray[k]);
}

void cloud::potcloud(cloud*) // calculating potenital energies
{
	for(int k = 0; k < natoms; k++)
	{
		particleArray[k].pot = 0;
	}

	for(int p1 = 0; p1 < natoms - 1; p1++)
	{
		for (int p2 = p1 + 1; p2 < natoms; p2++)
		{   
			if (p1 != p2 && distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)
			{
				potential(&particleArray[p1],&particleArray[p2]);
			}
		} 
	}
}

void cloud::posverletcloud(double timestep, cloud*) // new positions using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].posverlet(timestep, &particleArray[k]);
}

void cloud::velverletcloud(double timestep, cloud*) // new velocities using verlet
{
	for(int k = 0; k < natoms; k++)
		particleArray[k].velverlet(timestep, &particleArray[k]);
}

void cloud::radialdistfunc(int densbins, vector<double> *densArr, cloud*) // finds the density over distance
{
	vector<double> distances;

	for(int k = 0; k < natoms; k++)
	{   
		for(int l = 0; l < natoms; l++)
		{   
			if(l !=  k){distances.push_back(distance(&particleArray[k],&particleArray[l]));} // build a vector of all distances
		}
	}

	vector<int>hist(densbins);
	double shell1 = 0, shell2;
	for(int k = 0; k < densbins; k++)
	{   
		shell2 = ((double) k + 1)/densbins*boxlength;
		hist[k] = std::count_if(distances.begin(), distances.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;}); // count distances into bins
		shell1 = shell2;
	}

	shell1 = 0, shell2 = 0;
	double vol;
	for (int k = 0; k < densbins; k++)
	{   
		shell2 = 4/3*pi*pow(((double) k + 1)/densbins*boxlength,3); // calculade densities
		vol = shell2 - shell1;
		(*densArr)[k] = hist[k]*argonmass/vol;
		shell1 = shell2;
	}
}

void cloud::veldist(double maxspeed, int velbins, vector<double> *velArr, cloud*) // finds the velocity distrobution
{
	vector<double> velocities(natoms);

	for(int k = 0; k < natoms; k++)
	{   
		velocities[k] = particleArray[k].spd;
	}
	double shell1 = 0, shell2;
	for(int k = 0; k < velbins; k++)
	{   
		shell2 = ((double) k + 1)/velbins*maxspeed;
		(*velArr)[k] = std::count_if(velocities.begin(), velocities.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;}); // count velocities into bins
		shell1 = shell2;
	}
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//
// Printing functions

void cloud::printPositions(cloud*) // printing positions
{
	ofstream filestream;
	filestream.open("positions.txt");
	for(int k = 0; k < natoms; k++)
		filestream << particleArray[k].pox << " "  << particleArray[k].poy << " " << particleArray[k].poz << "\n";
	filestream.close();
}

void cloud::printSpeeds(cloud*) // printing speeds
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

void cloud::printDensArr(int densbins, vector<double> *densArr) // print velocity distribution
{
	ofstream filestream;
	filestream.open("radialdistfunc.txt");
	for(int k = 0; k < densbins; k++)
		filestream << (*densArr)[k]  << "\n";
	filestream.close();
}

void cloud::printSpeedDist(int velbins, vector<double> *velArr) // print velocity distribution
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
	cout << particleArray[k].pox << " " << particleArray[k].vex << " " << particleArray[k].fox << "\n";
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
// Main program



int main()
{

	// simulation specifics
	const double timestep = 1e-14;
	const int duration = 1000;

	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	int side = 6;
	int atoms = 4*pow(side,3);
	double length = 10.229*sigma;			// from paper
	double grid = length/side;				// initial separation
	
	// creating particle cloud
	cloud particleCloud(side, atoms, length, grid);
	
	// set initial conditions

	srand(time(NULL));	// setting random seed
	particleCloud.initPositions(&particleCloud);	// initializing positions
	particleCloud.initVelocities(&particleCloud);				// initializing velocities

	for(int k = 0; k< duration; k++)	// main simulation loop
	{
		particleCloud.forcecloud(&particleCloud);	// calculating force
		particleCloud.velverletcloud(timestep, &particleCloud);	// calculating new velocities
		particleCloud.posverletcloud(timestep, &particleCloud);	// caluclating new positions
		//particleCloud.display(275,&particleCloud); 
		particleCloud.forcecloud(&particleCloud);	// calculating forces
		particleCloud.velverletcloud(timestep, &particleCloud);	// calculating new velocities
	}

	particleCloud.speedcloud(&particleCloud);
	particleCloud.speedsquarecloud(&particleCloud);
	particleCloud.potcloud(&particleCloud);	// calculating potential energies
	particleCloud.kincloud(&particleCloud);				// calculating kinetic energies

	//	particleCloud.display(1,&particleCloud);	

	int velbins = 100;					// number of bins in speed distribution
	double maxspeed = 800.0;			// maximum speed
	vector<double> velArr (velbins);	// vector of all speeds
	particleCloud.veldist(maxspeed, velbins, &velArr, &particleCloud); // calculating speed distribution
	
	int densbins = 500;					// number of bins in density distribution
	vector<double> densArr (densbins);	// vector of density over distance
	particleCloud.radialdistfunc(densbins, &densArr, &particleCloud); // calculating radial distribution func
	
	particleCloud.printPositions(&particleCloud);	// print positions
	particleCloud.printSpeeds(&particleCloud);		// print raw speeds
	particleCloud.printSpeedDist(velbins, &velArr);			// print speed distribution
	particleCloud.printDensArr(densbins, &densArr);			// print density distribution

	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	cout << benchtime << "\n" ;

	// -----------------------------------------------------------------------
	// testing stuff
	
		




























	while(true) // interface for simple debugging
	{
		int k;
		cout << "give k\n";
		cin >> k;
		particleCloud.display(k,&particleCloud);
	}

	return 0;
}









