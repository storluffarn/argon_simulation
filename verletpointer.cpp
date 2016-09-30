// for benchmarking, does not perform better

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

//definitions

const double timestep = 5e-15;
const int duration = 1000; // in 1e-14s, number of timesteps

const int T = 120;
const double sigma = 3.4e-10;
const double kB = 1.38064852e-23;
const double epsilon = T*kB;
const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;


struct particle // all relevant particle properties goes here
{
    double posx;
    double vx;
	double kinE;
	double potE;
};

//initial conditions

double x1 = 0;
double x2 = 1.25e-11;
double v1 = 1;
double v2 = -2;

//1-d LJ-force

double LJforce(particle *p1, particle *p2)
{
	double r = fabs(p1 -> posx - p2 -> posx);
	double F = 24*epsilon*pow(sigma,6)*(pow(r,6)-2*pow(sigma,6))/pow(r,13);

	return F;
}

double LJpot(particle *p1, particle *p2) 
{
	double r = fabs(p1 -> posx - p2 -> posx);
	double V = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));

	return V;
}

// make sub rutines for printing, and energy calculations
//
//

void verlet (particle *p1, particle *p2)
{
	double forcenew, forceold;

	double Etot1;
	double Etot2;
	double Etot;

	std::ofstream verlet1dfile;
	verlet1dfile.open("verlet1dpointer.txt");

	for(int k = 0; k < duration; k++)
	{
		forceold = LJforce(p1,p2);

		p1 -> vx += timestep*forceold/2;
		p2 -> vx += timestep*(-forceold)/2; // - because no vector formulation
		
		p1 -> posx += timestep*p1 -> vx;
		p2 -> posx += timestep*p2 -> vx;
		
		forcenew = LJforce(p1,p2);

		p1 -> vx += timestep*forcenew/2;
		p2 -> vx += timestep*(-forcenew)/2;
		
		p1 -> kinE = argonmass * p1 -> vx * p1 -> vx/2;
		p1 -> potE = LJpot(p1,p2);
		
		p2 -> kinE = argonmass * p2 -> vx * p2 -> vx/2;
		p2 -> potE = LJpot(p1,p2);
		
		Etot1 = p1 -> kinE + p1 -> potE;
		Etot2 = p2 -> kinE + p2 -> potE;
		Etot = p1 -> kinE + p1 -> potE + p2 -> kinE + p2 -> potE;
		
		verlet1dfile << p1 -> posx << " " << p1 -> vx << " " << p1 -> kinE << " " << p1 -> potE << "\n";
		verlet1dfile << p2 -> posx << " " << p2 -> vx << " " << p2 -> kinE << " " << p2 -> potE << "\n";
		verlet1dfile << Etot1 << " " << Etot2 << " " << Etot << "\n";
	}
	verlet1dfile.close();
}

int main()
{	
	particle p1 = {x1,v1,0,0};
	particle p2 = {x2,v2,0,0};
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	verlet(&p1, &p2);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	
	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	std::cout << "execution in microseconds: " << benchtime  << "\n";

	return 0;
}












