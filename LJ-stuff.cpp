
#include <iostream> 
#include <fstream> // file stream
#include <cmath> 
#include <random> 

// Constants

const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;
const int T = 120;
const double pi = 3.14159265;
const double kB = 1.38064852e-23;
const double sigma = 3.4e-10; //?
const double epsilon = 120*kB; //?
const double dcuttoff =10*pow(2,1/6)*sigma; // 10 * cutoff for interactions

//Box

const int boxside = 6; // number of atoms on each side
const int boxsize = 4*6*6*6; // total number of atoms (face centred packing)
const double boxlength = 10.229*sigma; //from paper
const double gridsize = boxlength/boxside; // spatial separation in initialization


struct particle // all relevant particle properties goes here
{
	double posx;
	double posy;
	double posz;
	double vx;
	double vy;
	double vz;
};

particle particleArray[boxsize];

// structures to be used for returning stuff

struct threevals // used to return three values
{
	double val1, val2, val3;
};

struct fourvals // used to return four values
{
	double val1, val2, val3, val4;
};

// functions

// initializes the positions, separation is specified by gridsize above.
// each loop places 4 particles corresponding to face centered crystal packing

void initPositions() 
{
	int k = 0;
	
	for (int kz = 0; kz < boxside; kz++)
		for (int ky = 0; ky < boxside; ky++)
			for (int kx = 0; kx < boxside; kx++, k += 4)
			{
				double displace = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				displace *= gridsize;
				displace = 0; // uncomment to remove noise

				particleArray[k].posx = kx*gridsize + displace;
				particleArray[k].posy = ky*gridsize + displace;
				particleArray[k].posz = kz*gridsize + displace;
				particleArray[k+1].posx = (kx + 0.5)*gridsize + displace;
				particleArray[k+1].posy = (ky + 0.5)*gridsize + displace;
				particleArray[k+1].posz = kz*gridsize + displace;
				particleArray[k+2].posx = (kx + 0.5)*gridsize + displace;
				particleArray[k+2].posy = ky*gridsize + displace;
				particleArray[k+2].posz = (kz + 0.5)*gridsize + displace;
				particleArray[k+3].posx = kx*gridsize + displace;
				particleArray[k+3].posy = (ky + 0.5)*gridsize + displace;
				particleArray[k+3].posz = (kz + 0.5)*gridsize + displace;
			}
}

// assigning random velicities to the particles, doing this by assigning random
// speeds in x, y, z directions. these speeds are normal distributed with
// centre and standard deviation as given below.

void initVelocities()
{
	double mean = 0; // centre of 1d velocities, not same as post prob speed!
	double std = sqrt(kB*T/argonmass); // standard deviation
	
	std::default_random_engine generator; // specifiying random number generator
	std::normal_distribution<double> veldist(mean,std); // distribution

	for (int k = 0; k < boxsize; k++)
		{
			particleArray[k].vx = veldist(generator);
			particleArray[k].vy = veldist(generator);
			particleArray[k].vz = veldist(generator);
		}
}

double distance (int p1, int p2)
{
	double x,y,z;

	x = fabs(particleArray[p1].posx - particleArray[p2].posx); 
	y = fabs(particleArray[p1].posy - particleArray[p2].posy); 
	z = fabs(particleArray[p1].posz - particleArray[p2].posz);

	x -= static_cast<int> (x/boxlength + 0.5) * boxlength;
	y -= static_cast<int> (y/boxlength + 0.5) * boxlength;
	z -= static_cast<int> (z/boxlength + 0.5) * boxlength;

	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

	return r;
}

threevals direction (int p1, int p2)
{
	threevals returnvals;
	double xdir,ydir,zdir;

	xdir = particleArray[p1].posx - particleArray[p2].posx; 
	ydir = particleArray[p1].posy - particleArray[p2].posy; 
	zdir = particleArray[p1].posz - particleArray[p2].posz;
	
	returnvals.val1 = xdir;
	returnvals.val2 = ydir;
	returnvals.val3 = zdir;

	return returnvals;
}

double norm(threevals v)
{
	double norm = sqrt(v.val1*v.val1 + v.val2*v.val2 + v.val3*v.val3);
	
	return norm;
}

double LJpot(int p1, int p2)
{
	double r = distance(p1, p2);
    double V = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
			 
    return V;
}

double LJforce(int p1, int p2) // returns the magnitude of the LJ-force
{
	double r = distance(p1, p2);
    double F = 24*epsilon*pow(sigma,6)*(pow(r,6)-2*pow(sigma,6))/pow(r,13);
	
	return F;
}

fourvals LJforceTot (int p1)
{
	double getforce;
	fourvals totforce = {};
	threevals getdir;
	double getnorm;

	for (int p2 = 0; p2 < boxsize; p2++)
	{	
		if (p1 != p2 && distance(p1,p2) <= dcuttoff)
		{
			getforce = LJforce(p1,p2);
			getdir = direction(p1,p2);
			getnorm = norm(getdir);
		
			getdir.val1 /= getnorm;
			getdir.val2 /= getnorm;
			getdir.val3 /= getnorm;
			
			totforce.val1 += getforce*getdir.val1;
			totforce.val2 += getforce*getdir.val2;
			totforce.val3 += getforce*getdir.val3;
		}
		totforce.val4 = sqrt(pow(totforce.val1,2) + pow(totforce.val2,2) + pow(totforce.val3,2));
	}

	return totforce;
}


int main()
{
	int k;
	fourvals F;

	initPositions();
	initVelocities();

	while(true)
	{
		std::cout << "Give particle index:\n";
		std::cin >> k;

		F = LJforceTot(k);

		std::cout << "Lennard-Jones Force:\n" << F.val4 << "\n";
		std::cout << "Lennard-Jones x-direction:\n" << F.val1 << "\n";
		std::cout << "Lennard-Jones y-diirection:\n" << F.val2 << "\n";
		std::cout << "Lennard-Jones z-direction:\n" << F.val3 << "\n";	
	}

	return 0;
}























