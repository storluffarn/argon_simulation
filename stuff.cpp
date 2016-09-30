#include <iostream>

struct particle // all relevant particle properties goes here
{
    double posx;
    double vx; 
    double kinE;
    double potE;
};

//initial conditions

double x1 = 0;
double x2 = 1.0e-11;
double v1 = 0;
double v2 = -1; 

particle p1 = {x1,v1,0,0};
particle p2 = {x2,v2,0,0};

int main()
{
	p2.posx += 1.0e-11;
	std::cout << p2.posx << " " << p2.vx <<"\n";
	return 0;
}
