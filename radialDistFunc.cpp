
#include <iostream> // for io
#include <fstream> //for file
#include <cmath> //for time
#include <random> // for random
#include <ctime> //for measuring time?
#include <cstdlib> //for srand (random seed)?
#include <algorithm> // for count
#include <vector> //for well, vector

// Constants

const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;
const int T = 120;
const double kB = 1.38064852e-23;
const double sigma = 3.4e-10; // interaction cut off
const double pi = atan(1)*4;

//Box

const int boxside = 6;
const int boxsize = 4*6*6*6; // also number of particles, should change name?
const double boxlength = 10.229*sigma; //from paper
const double gridsize = boxlength/boxside; // spatial separation in initialization

struct particle // all relevant particle properties goes here
{
	double posx;
	double posy;
	double posz;
};

particle particleArray[boxsize];

// initializes the positions, separation is specified by gridsize above.
// each loop plaves 4 particles corresponding to face centered crystal packing

void initPositions() 
{
	int k = 0;
	std::srand(std::time(0)); // user current time as random seed

	for (int kz = 0; kz < boxside; kz++)
		for (int ky = 0; ky < boxside; ky++)
			for (int kx = 0; kx < boxside; kx++, k += 4)
			{
				
				double displace = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);	
				displace *= gridsize;
				//displace = 0; // uncomment to remove noise
				
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

double dist (int p1, int p2)
{
	double x = fabs(particleArray[p1].posx - particleArray[p2].posx); 
	double y = fabs(particleArray[p1].posy - particleArray[p2].posy); 
	double z = fabs(particleArray[p1].posz - particleArray[p2].posz);

	x -= static_cast<int> (x/boxlength + 0.5) * boxlength;
	y -= static_cast<int> (y/boxlength + 0.5) * boxlength;
	z -= static_cast<int> (z/boxlength + 0.5) * boxlength;

	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

	return r;
}

const int bins = 1000;
double densArr[bins];
std::vector<double> distances;

void radialdistfunc() //counts the number of particles in "shells" around x_i
{
	for(int k = 0; k < boxsize; k++)
	{
		for(int l = 0; l < boxsize; l++)
		{
			if(l !=  k){distances.push_back(dist(k,l));}
		}
	}
	
	int hist[bins];
	double shell1 = 0, shell2;
	for(int k = 0; k < bins; k++)
	{
		shell2 = ((double) k + 1)/bins*boxlength;
		hist[k] = std::count_if(distances.begin(), distances.end(), [&shell1, &shell2](double i){return shell1 < i && i < shell2;});
		shell1 = shell2;
		std::cout << hist[k] << '\n';
	}

	shell1 = 0, shell2 = 0;
	double vol;
	for (int k = 0; k < bins; k++)
	{
		shell2 = 4/3*pi*pow(((double) k + 1)/bins*boxlength,3);
		vol = shell2 - shell1;
		densArr[k] = hist[k]*argonmass/vol;
		shell1 = shell2;
	}
}


void printInitPositions()
{
	std::ofstream filestream;
	filestream.open("radialdistpos.txt");
	for(int k = 0; k < boxsize; k++)
		filestream << particleArray[k].posx << " "  << particleArray[k].posy << " " << particleArray[k].posz << " \n";
	filestream.close();
}

void printDensArr()
{
	std::ofstream filestream;
	filestream.open("radialdistfunc.txt");
	for(int k = 0; k < bins; k++)
		filestream << k << " " << densArr[k] << "\n";
	filestream.close();
}

int main()
{
	int k;

	initPositions();
	printInitPositions();
	radialdistfunc();
	printDensArr();

	while(true)
	{
		std::cout << "Give bin index:\n";
		std::cin >> k;
		std::cout << "density in this bin is: \n" << distances[k] << "\n";
	}

	return 0;
}
















