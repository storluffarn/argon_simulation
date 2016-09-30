
// Assignment for Astrid's simulation course. Molecular dynamics of an argon liquid.

// See also notes.txt 

// Preliminaries

#include <cmath>
#include <iostream>
#include <cstdlib>

#define  k_B = 1.380649e-23 // Boltzmann's constant
#define sigma = 3.4e-10 // distance at which the potential is zero
#define epsilon = 120*k_B // depth of potentail wall
#define boxlength = 10.229*sigma
// #define umass = 1.660539E-27 // atomic mass unit

// Funcions

// Setting up box

void buildbox()
{
	
}

// Stepping function

// Lennard Jones Potential, takes particle coordinates gives potential

double LJpot(double r_x, double r_y, double r_z)
{
	double r = cbrt(r_x*r_x + r_y*r_y + r_z*r_z);
	double V = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));

	return V;
}

// Lennard Jones Force, D[V(r),r] evaluated and simplified with mathematica

double LJforce(double r_x, double r_y, double r_z)
{
	double r = cbrt(r_x*r_x + r_y*r_y + r_z*r_z);
	double F = 24*epsilon*pow(sigma,6)*(pow(r,6)-2*pow(sigma,6))/pow(r,13);

	return F;
}

// Main program

int main()
{
	// Defining stuff
	double k,l,V,F;
	
	// Initilizing
	buildbox
	
	// Calculating stuff
	std::cout << "Give position for first particle\n";
	std::cin >> k;
	std::cout << "Give position for second particle\n";
	std::cin >> l;
	V=LJpot(k,l);
	F=LJforce(k,l);
	std::cout << "Lennard Jones potential:\n" << V << "\n" << "Lennard Jones force:\n" << F << "\n";

	return 0;
}






















