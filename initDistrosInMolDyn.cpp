
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <random> 

// Constants

const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;
const int T = 120;
const double pi = 3.14159265;
const double kB = 1.38064852e-23;
const double sigma = 3.4e-10; //?

//Box

const int boxside = 6;
const int boxsize = 4*6*6*6;
double gridsize = 10.229*sigma/boxside; // spatial separation in initialization

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

// initializes the positions, separation is specified by gridsize above.
// each loop plaves 4 particles corresponding to face centered crystal packing

void initPositions() 
{
	int k = 0;
	
	for (int kz = 0; kz < boxside; kz++)
		for (int ky = 0; ky < boxside; ky++)
			for (int kx = 0; kx < boxside; kx++, k += 4)
			{
				particleArray[k].posx = kx*gridsize;
				particleArray[k].posy = ky*gridsize;
				particleArray[k].posz = kz*gridsize;
				particleArray[k+1].posx = (kx + 0.5)*gridsize;
				particleArray[k+1].posy = (ky + 0.5)*gridsize;
				particleArray[k+1].posz = kz*gridsize;
				particleArray[k+2].posx = (kx + 0.5)*gridsize;
				particleArray[k+2].posy = ky*gridsize;
				particleArray[k+2].posz = (kz + 0.5)*gridsize;
				particleArray[k+3].posx = kx*gridsize;
				particleArray[k+3].posy = (ky + 0.5)*gridsize;
				particleArray[k+3].posz = (kz + 0.5)*gridsize;
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

void printInitPositions()
{
	std::ofstream initPosFile;
	initPosFile.open("initPositions.txt");
	for(int k = 0; k < boxsize; k++)
		initPosFile << particleArray[k].posx << " "  << particleArray[k].posy << " " << particleArray[k].posz << " \n";
	initPosFile.close();
}

int main()
{
	int k;

	initPositions();
	initVelocities();
	printInitPositions();

	while(true)
	{
		std::cout << "Give particle index:\n";
		std::cin >> k;
		std::cout << "x position:\n" << particleArray[k].posx << "\n";
		std::cout << "y position:\n" << particleArray[k].posy << "\n";
		std::cout << "z position:\n" << particleArray[k].posz << "\n";
		std::cout << "x velocity:\n" << particleArray[k].vx << "\n";
		std::cout << "y velocity:\n" << particleArray[k].vy << "\n";
		std::cout << "z velocity:\n" << particleArray[k].vz << "\n";
	}

	return 0;
}























