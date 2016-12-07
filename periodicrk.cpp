
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>

using namespace std;

//definitions

const double timestep = 7e-14;
const int duration = 1000; // in 1e-14s, number of timesteps

const int T = 120;
const double sigma = 3.4e-10;
const double kB = 1.38064852e-23;
const double epsilon = T*kB;
const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;
const int boxside = 6;
const double boxside2 = 10.299*sigma;
const double gridsize = boxside2/boxside;
const double npart = 2;

struct particle // all relevant particle properties goes here
{
    double pox;
    double vex;
	double acx;
	double kinE;
	double potE;
};

struct virtualpart
{
	double pox;
	double vex;
	double acx;

	public:
	virtualpart(double a, double b, double c)
		: pox(a), vex(b), acx(c)
	{}
};

//initial conditions

double x1 = boxside2-10.0e-10;
double x2 = boxside2-0e-14;
double v1 = 100;
double v2 = -100;

particle p1 = {x1,v1,0,0,0};
particle p2 = {x2,v2,0,0,0};

//1-d LJ-force

double LJforce()
{
	double r = fabs(p2.pox - p1.pox);
	r -= static_cast<int> (r/boxside2 + 0.5) * boxside2; // +0.5 for cast to work // static_cast<int>(n >= 0 ? n + 0.5 : n - 0.5) for nearest int
	
	double F = - 24*epsilon*(-pow(sigma,6)/pow(r,7) + 2*pow(sigma,12)/pow(r,13) ); // works and correct

	return F;
}

double LJpot() 
{
	double r = fabs(p2.pox - p1.pox);
	r -= static_cast<int> (r/boxside2 + 0.5) * boxside2; 

	double V = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));

	return V;
}

double calcacc(double r)
{
	r -= static_cast<int> (r >= 0 ? r/boxside2 + 0.5 : r/boxside2 - 0.5)*boxside2;

	// r -= static_cast<int> (r/boxside2 + 0.5) * boxside2; // +0.5 for cast to work // static_cast<int>(n >= 0 ? n + 0.5 : n - 0.5) for nearest int
	
	double F = - 24*epsilon*(-pow(sigma,6)/pow(r,7) + 2*pow(sigma,12)/pow(r,13) ); // works and correct
	//double F = - 4*epsilon*(-6*pow(sigma,6)*pow(r,-7) + 12*pow(sigma,12)*pow(r,-13) );
	
	return F/argonmass;
}

//pair <double,double> rkstep(double pos, double vel, double acc, double time) // might want to not return vector if we want to use to do pos and vel steps separetly?
//{
//	pair <double, double> rkstep;
//
//	rkstep.first = pos + vel * time;
//	rkstep.second = vel + acc * time;
//
//	return rkstep;
//}

double qstep(double q, double qdot, double time) // might want to not return vector if we want to use to do pos and vel steps separetly?
{
	return q + qdot * time;
}

void rk4(particle* p1, particle* p2, double timestep, double time)
{
	double kpos, kvel, kacc;

	kpos = qstep(p1->pox,p1->vex,time);
	kvel = qstep(p1->vex,p1->acx,time);
	
	virtualpart k11 (p1->pox, p1->vex, p1->acx);
	
	kpos = qstep(p2->pox,p2->vex,time);
	kvel = qstep(p2->vex,p2->acx,time);
	
	virtualpart k12 (p2->pox, p2->vex, p2->acx);
	
	kacc = calcacc(k11.pox-k12.pox);
	k11.acx = kacc;
	k12.acx = -kacc;
	
	kpos = qstep(p1->pox,0.5*timestep*k11.vex,time+0.5*timestep);
	kvel = qstep(p1->vex,0.5*timestep*k11.acx,time+0.5*timestep);
	
	virtualpart k21 = {kpos,kvel,0};
	
	kpos = qstep(p2->pox,0.5*timestep*k12.vex,time+0.5*timestep);
	kvel = qstep(p2->vex,0.5*timestep*k12.acx,time+0.5*timestep);
	
	virtualpart k22 = {kpos,kvel,0};

	kacc = calcacc(k21.pox-k22.pox);
	k21.acx = kacc;
	k22.acx = -kacc;

	kpos = qstep(p1->pox,0.5*timestep*k21.vex,time+0.5*timestep);
	kvel = qstep(p1->vex,0.5*timestep*k21.acx,time+0.5*timestep);
	
	virtualpart k31 = {kpos,kvel,0};
	
	kpos = qstep(p2->pox,0.5*timestep*k22.vex,time+0.5*timestep);
	kvel = qstep(p2->vex,0.5*timestep*k22.acx,time+0.5*timestep);
	
	virtualpart k32 = {kpos,kvel,0};
		
	kacc = calcacc(k31.pox-k32.pox);
	k31.acx = kacc;
	k32.acx = -kacc;	
	
	kpos = qstep(p1->pox,timestep*k31.vex,time+timestep);
	kvel = qstep(p1->vex,timestep*k31.acx,time+timestep);
	
	virtualpart k41 = {kpos,kvel,0};
	
	kpos = qstep(p2->pox,timestep*k32.vex,time+timestep);
	kvel = qstep(p2->vex,timestep*k32.acx,time+timestep);
	
	virtualpart k42 = {kpos,kvel,0};
	
	kacc = calcacc(k41.pox-k42.pox);
	k41.acx = kacc;
	k42.acx = -kacc;

	p1->pox += timestep/6.0 * (k11.vex + k21.vex + k31.vex + k41.vex);
	p1->vex += timestep/6.0 * (k11.acx + k21.acx + k31.acx + k41.acx);
	p1->acx = timestep/6.0*(k11.acx + k21.acx + k31.acx + k41.acx);
	p2->pox += timestep/6.0 * (k12.vex + k22.vex + k32.vex + k42.vex);
	p2->vex += timestep/6.0 * (k12.acx + k22.acx + k32.acx + k42.acx);
	p2->acx = timestep/6.0*(k12.acx + k22.acx + k32.acx + k42.acx);


	//cout << timestep/6.0 * (k11.vex + k21.vex + k31.vex + k41.vex) << endl;
	
	//cout << p1->pox << endl;
}


// make sub rutines for printing, and energy calculations
//
//

void verlet ()
{
	//double force;
	double acc;
	double time;

	double Etot1;
	double Etot2;
	double Etot;
	
	acc = calcacc(p1.pox - p2.pox);
	p1.acx = acc;
	p1.acx = -acc;

	std::ofstream verlet1dfile;
	verlet1dfile.open("periodicrk.txt");

	for(int k = 0; k < duration; k++)
	{
		time = timestep*k;
		
		//cout << p1.pox << endl;
		
		rk4(&p1, &p2, timestep, time);

		//cout << p1.pox << endl;
				
		//p1.vx = rk4(p1.vx, acc, timestep, time);
		//p2.vx = rk4(p2.vx, -acc, timestep, time);

		//p1.posx = rk4(p1.posx,p1.vx,timestep, timestep);
		//p2.posx = rk4(p2.posx,p2.vx,timestep, time);
		
		//force = LJforce(); // values 1e-19 - 1e-11
		
		//p1.vx += timestep*force/argonmass*0.5;
		//p2.vx += -timestep*force/argonmass*0.5;

		//p1.posx += timestep*p1.vx;
		//p2.posx += timestep*p2.vx;

		//force = LJforce();
		
		//p1.vx += timestep*force/argonmass*0.5;
		//p2.vx += -timestep*force/argonmass*0.5;
		
		p1.kinE = argonmass*p1.vex*p1.vex/2;
		p1.potE = LJpot();
		
		p2.kinE = argonmass*p2.vex*p2.vex/2;
		p2.potE = LJpot();
		
		Etot1 = p1.kinE + p1.potE;
		Etot2 = p2.kinE + p2.potE;
		Etot = p1.kinE + p1.potE - (p2.kinE + p2.potE);
		
		verlet1dfile << p1.pox << " " << p1.vex << " " << p1.kinE << " " << p1.potE << "\n";
		verlet1dfile << p2.pox << " " << p2.vex << " " << p2.kinE << " " << p2.potE << "\n";
		verlet1dfile << Etot1 << " " << Etot2 << " " << Etot << "\n";
		verlet1dfile << p1.acx << "\n";
	}
verlet1dfile.close();
}

int main()
{	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	verlet();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	std::cout << "execution in microseconds: " << benchtime << "\n";

	return 0;
}













