
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
#include <ctime>				// for random seed time

using namespace std;		// tidy up

// constants

const double umass = 1.660539040e-27;				// atomic mass
const double argonmass = umass*39.948;				// argon mass
const double T = 94.4;								// temperature
const double pi = atan(1)*4;						// for using pi
const double kB = 1.38064852e-23;					// Boltzmann's constant
const double sigma = 3.4e-10;						// interaction term
//const double sigma = 0;						// interaction term
const double epsilon = 120*kB;						// energy term
const double dcutoff = 2.25*sigma;			// interaction cutoff 
//const double dcutoff = 1;  
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
	
	
};

class cloud // class for the set of all particles and functions on these
{
	int natomsside;
	int natoms;
	double boxlength;			// from paper
	double gridsize;	// initial separation

	unsigned int ugly = natoms;		// why does it need to be unsigned?
	
	// object member
	
	public:
	vector<particle> particleArray{ugly};						// vector of all particles

	cloud(int a, unsigned int b, double c, double d)
		: natomsside(a), natoms(b), boxlength(c), gridsize(d)
	{}	
	
	// accessors
	void setpos(int,double,double,double,cloud*);		// for setting parameters
	void setvel(int,double,double,double,cloud*);		
	void setforce(int,double,double,double,cloud*);
	void setpot(int,double,cloud*);
	vector<double> getpos(int,cloud*);				// for getting parameters			
	vector<double> getvel(int,cloud*);	
	double getpot(int,cloud*);
	double getkin(int,cloud*);

	//	
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

void cloud::setpos(int i, double k, double l, double m, cloud*)		// for setting positons
{
	particleArray[i].pox = k;
	particleArray[i].poy = l;
	particleArray[i].poz = m;
}

void cloud::setvel(int i, double k, double l, double m, cloud*)		// for setting velocities
{
	particleArray[i].vex = k;
	particleArray[i].vey = l;
	particleArray[i].vez = m;
}

void cloud::setforce(int i, double k, double l, double m, cloud*)		// for setting positons
{
	particleArray[i].fox = k;
	particleArray[i].foy = l;
	particleArray[i].foz = m;
}

void cloud::setpot(int i, double k, cloud*)		// for setting positons
{
	particleArray[i].pot = k;
}

vector<double> cloud::getpos(int i, cloud*)		// for getting positions
{
	vector<double> returnvals = {particleArray[i].pox, particleArray[i].poy, particleArray[i].poz};

	return returnvals;
}

vector<double> cloud::getvel(int i, cloud*)		// for getting velocities
{
	vector<double> returnvals = {particleArray[i].vex, particleArray[i].vey, particleArray[i].vez};

	return returnvals;
}

double cloud::getkin(int i, cloud*)
{
	return particleArray[i].kin;
}

double cloud::getpot(int i, cloud*)
{
	return particleArray[i].pot;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
//Initializing


void cloud::initPositions(cloud* c) // setting initial positions with face centered packing
{		
	int k = 0;
	for (int kz = 0; kz < natomsside; kz++)
		for (int ky = 0; ky < natomsside; ky++)
			for (int kx = 0; kx < natomsside; kx++, k += 4) // placing four particles a time
			{   
				double displace = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				displace = 1; // uncomment to add white noise

				setpos(k,kx*gridsize*displace, ky*gridsize*displace, kz*gridsize*displace, c);
				setpos(k + 1, (kx + 0.5)*gridsize*displace, (ky + 0.5)*gridsize*displace, kz*gridsize*displace, c);
				setpos(k + 2, (kx + 0.5)*gridsize*displace, ky*gridsize*displace, (kz + 0.5)*gridsize*displace, c);
				setpos(k + 3, kx*gridsize*displace, (ky + 0.5)*gridsize*displace, (kz + 0.5)*gridsize*displace, c);
			}   
}

void cloud::initVelocities(cloud* c) // setting initial velocities, assuming approx gaussian
{
	double mean = 0; // centre of 1d velocities
	double std = sqrt(kB*T/argonmass); // standard deviation

	default_random_engine generator; // random number generator to be used
	normal_distribution<double> veldist(mean,std); // probability distribution

	for (int k = 0; k < natoms; k++) // setting velocites
	{
		setvel(k,veldist(generator), veldist(generator), veldist(generator), c);
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
// Functions on the particles

vector<double> cloud::direction (particle* p1, particle* p2)
{
	vector<double> direction;
	
	double x,y,z;
	
	x = p1->pox - p2->pox;			// 1d distances
	y = p1->poy - p2->poy;
	z = p1->poz - p2->poz;
		
	x -= static_cast<int> (x >= 0 ? z/boxlength + 0.5 : x/boxlength - 0.5)*boxlength;
	y -= static_cast<int> (y >= 0 ? z/boxlength + 0.5 : y/boxlength - 0.5)*boxlength;
	z -= static_cast<int> (z >= 0 ? z/boxlength + 0.5 : z/boxlength - 0.5)*boxlength;

	//cout << "dist1 " << dx << " " << dy << " " << dz << "\n";

	direction = {x,y,z};

	//cout << "dist2 " << direction[0] << " " << direction[1] << " " << direction[2]	<< "\n";

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
	double F = 4*epsilon*(-6*pow(sigma,6)/pow(r,4) + 12*pow(sigma,12)/pow(r,7) );
	
	//double r = distance(p1, p2); // see distance function
	//double F = 24*epsilon*(-pow(sigma,6)*pow(r,8) + 2*pow(sigma,12)/pow(r,14) );

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

	//cout << p1->foz << "\n";
}


void particle::posverlet(double timestep,particle* p1) // taking a verlet position step
{
	p1->pox += timestep*p1->vex;
	p1->poy += timestep*p1->vey;
	p1->poz += timestep*p1->vez;
}

void particle::velverlet(double timestep, particle* p1) // taking a verlet velocity step, also calculating accelerations
{	
	p1->vex += timestep*p1->fox/argonmass*0.5;
	p1->vey += timestep*p1->foy/argonmass*0.5;
	p1->vez += timestep*p1->foz/argonmass*0.5;
}

void particle::kinetic(particle* p1) // calculates the kinetic energy of a particle
{
	p1->kin = 0.5*argonmass*p1->spdsqr;
}

void cloud::potential(particle* p1, particle* p2)
{
	double r = distancesquare(p1, p2);
	double pot = 4*epsilon*(pow(sigma,12)/pow(r,6)-pow(sigma,6)/pow(r,3)); // LJ-potential with squared ristance
	
	//cout << p1-> pot << "\n";
	p1->pot += pot;
	p2->pot += pot;
	//cout << p1-> pot << "\n";
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

	//cout << particleArray[0].fox << "\n";
	
	for(int p1 = 0; p1 < natoms - 1; p1++)
	{
		for(int p2 = p1 + 1; p2 < natoms; p2++)
		{
			//cout << p1 << " " << p2 << "\n";
			if (distancesquare(&particleArray[p1],&particleArray[p2]) <= dcutoffsquare)
			{
				force(&particleArray[p1],&particleArray[p2]);
				//cout << distancesquare(&particleArray[p1],&particleArray[p2]) << " " << dcutoffsquare << "\n";
				//cout << p1 << " " << p2 << "\n";
				//cout << particleArray[p2].fox << " " << particleArray[p2].foy << " " << particleArray[p2].foz << "\n";
			}
		}
	}

	//cout << particleArray[0].fox << "\n";
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
	double timestep = 5e-14;
	int duration = 1000;

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

	//particleCloud.speedcloud(&particleCloud);
	//particleCloud.speedsquarecloud(&particleCloud);
	//particleCloud.potcloud(&particleCloud);	// calculating potential energies
	//particleCloud.kincloud(&particleCloud);				// calculating kinetic energies

	//	particleCloud.display(1,&particleCloud);	

	//int velbins = 100;					// number of bins in speed distribution
	//double maxspeed = 800.0;			// maximum speed
	//vector<double> velArr (velbins);	// vector of all speeds
	//particleCloud.veldist(maxspeed, velbins, &velArr, &particleCloud); // calculating speed distribution
	
	//int densbins = 500;					// number of bins in density distribution
	//vector<double> densArr (densbins);	// vector of density over distance
	//particleCloud.radialdistfunc(densbins, &densArr, &particleCloud); // calculating radial distribution func
	
	particleCloud.printPositions(&particleCloud);	// print positions
	//particleCloud.printSpeeds(&particleCloud);		// print raw speeds
	//particleCloud.printSpeedDist(velbins, &velArr);			// print speed distribution
	//particleCloud.printDensArr(densbins, &densArr);			// print density distribution

	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	cout << benchtime << "\n" ;

	// -----------------------------------------------------------------------
	// testing stuff XXXXXXXXXXXXXXXXXXXXXXXX REMEMBER TO MAKE PARICLE ARRAY PRIVATE AGAIN 
	
	
	//int side = 6;
	//int atoms = 2;
	//double length = 10.229*sigma;			// from paper
	//double grid = length/side;				// initial separation
	//
	//// creating particle cloud
	//cloud testcloud(side, atoms, length, grid);
	//
	////checked and works
	////testcloud.setpos(0,length-15e-10,length-15e-10,length-15e-10,&testcloud);
	////testcloud.setpos(0,length-15e-10,length-15e-10,length-15e-10,&testcloud);
	//
	////checked and works
	//testcloud.setpos(0,0,0,0,&testcloud);
	//testcloud.setpos(1,0,0,length-2*sigma,&testcloud);

	////testcloud.setpos(0,0,0,0,&testcloud);
	////testcloud.setpos(1,length-sigma,0,0,&testcloud);

	//testcloud.setvel(0,0,0,0,&testcloud);
	//testcloud.setvel(1,0,0,0,&testcloud);

	//int duration = 1000;
	//double timestep = 5e-14;

	//vector <double> ppos(3);
	//vector <double> vvel(3);
	//double kin, pot;

	//ofstream filestream;
	//filestream.open("testpositions.txt");

	//testcloud.setforce(0,0,0,0,&testcloud);
	//testcloud.setforce(1,0,0,0,&testcloud);
	//if(testcloud.distancesquare(&testcloud.particleArray[0], &testcloud.particleArray[1]) <= dcutoffsquare)
	//	testcloud.force(&testcloud.particleArray[0], &testcloud.particleArray[1]);
	//
	//for(int k = 0; k< duration; k++)	// main simulation loop
	//{			
	//	testcloud.velverletcloud(timestep, &testcloud);	// calculating new velocities
	//	testcloud.posverletcloud(timestep, &testcloud);	// caluclating new positions
	//	
	//	testcloud.setforce(0,0,0,0,&testcloud);
	//	testcloud.setforce(1,0,0,0,&testcloud);
	//	
	//	if(testcloud.distancesquare(&testcloud.particleArray[0], &testcloud.particleArray[1]) <= dcutoffsquare)
	//		testcloud.force(&testcloud.particleArray[0], &testcloud.particleArray[1]);
	//	
	//	testcloud.velverletcloud(timestep, &testcloud);	// calculating new velocities
	//	
	//	ppos = testcloud.getpos(0,&testcloud);
	//	vvel = testcloud.getvel(1,&testcloud);

	//	testcloud.setpot(0,0,&testcloud);
	//	testcloud.setpot(1,0,&testcloud);
	//	testcloud.potential(&testcloud.particleArray[0], &testcloud.particleArray[1]);

	//	testcloud.speedsquarecloud(&testcloud);		
	//	testcloud.kincloud(&testcloud);

	//	
	//	ppos = testcloud.getpos(0,&testcloud);
	//	vvel = testcloud.getvel(0,&testcloud);
	//	kin = testcloud.getkin(0,&testcloud);
	//	pot = testcloud.getpot(0,&testcloud);
	//	filestream << ppos[0] << " "  << ppos[1] << " " << ppos[2] << "\n";
	//	filestream << vvel[0] << " "  << vvel[1] << " " << vvel[2] << "\n";
	//	filestream << kin << " "  << pot << "\n";
	//	
	//	ppos = testcloud.getpos(1,&testcloud);
	//	vvel = testcloud.getvel(1,&testcloud);
	//	kin = testcloud.getkin(0,&testcloud);
	//	pot = testcloud.getpot(0,&testcloud);
	//	filestream << ppos[0] << " "  << ppos[1] << " " << ppos[2] << "\n";
	//	filestream << vvel[0] << " "  << vvel[1] << " " << vvel[2] << "\n";
	//	filestream << kin << " "  << pot << "\n";
	//	
	//	testcloud.display(0,&testcloud);
	//}

	//filestream.close();
	//
	//testcloud.speedcloud(&testcloud);


	//----------------------------------------------------------------------------------------------
	// Testing more stuff -------- 

	//int side = 1;
	//int atoms = 3;;
	////int atoms = 4*pow(side,3);;
	//double length = 20.229*sigma;			// from paper
	//double grid = length/side;				// initial separation
	//
	//// creating particle cloud
	//cloud testcloud(side, atoms, length, grid);

	//// 1295 1e-15 make badness!	
	//int duration = 1000;
	//double timestep = 1e-14;

	//srand(time(NULL));	// setting random seed
	//
	////testcloud.initPositions(&testcloud);
	////testcloud.initVelocities(&testcloud);

	//testcloud.setpos(0,length*0.5,length*0.5,length*0.5,&testcloud);
	////testcloud.setpos(0,0,0,5*sigma,&testcloud);
	//testcloud.setpos(1,0,0,0,&testcloud);
	//testcloud.setpos(2,0,0,1*sigma,&testcloud);

	//testcloud.setvel(0,0,0,0,&testcloud);
	//testcloud.setvel(1,0,0,0,&testcloud);
	//testcloud.setvel(2,0,0,0,&testcloud);

	//testcloud.forcecloud(&testcloud);	// calculating force
	//for(int k = 0; k< duration; k++)	// main simulation loop
	//{	
	//	testcloud.velverletcloud(timestep, &testcloud);	// calculating new velocities
	//	testcloud.posverletcloud(timestep, &testcloud);	// caluclating new positions
	//	testcloud.display(1,&testcloud);
	//	testcloud.forcecloud(&testcloud);	// calculating forces
	//	testcloud.velverletcloud(timestep, &testcloud);	// calculating new velocities

	//	//cout << k << "\n";
	//}

	//testcloud.printPositions(&testcloud);



	//while(true) // interface for simple debugging
	//{
	//	int k;
	//	cout << "give k\n";
	//	cin >> k;
	//	particleCloud.display(k,&particleCloud);
	//}

	return 0;
}









