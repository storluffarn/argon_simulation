// Stuff in this file:
// *Constructing 2d arrays and 2d vectors
// *Using pointers to 2d arrays and 2d vectors
// *Assigning random numbers to elements in matrices
// *Measuring time to construct and assign values to 2d vectors and 2d arrays
//
// *Conclusion: Arrays are faster.
//

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib> //for random number seed?
#include <chrono> //for measuring times

using namespace std;
using namespace std::chrono;

//// Accessing the array with a pointer WORKING
//
//int main()
//{
//	int array[5];
//	for (int i = 0; i < 5; i++)
//	{
//		array[i]=i*2;
//	}
//
//	int k;
//
//	cout << "Which element of the 1-d array?\n";
//	cin >> k;
//	cout << *(array+k) << "\n";
//
//	return 0;
//}

double sigma = 3.4e-10;

double boxlength = 10.229*sigma;
const int hi = 10000, wi = 3;
double arrbox[hi][wi];
vector< vector<double> > vecbox(hi, vector<double>(wi));

void buildarrbox()
{
     int m=0, k, l;
//	 srand(static_cast<unsigned int>(time(0)));
	 srand(2);
     
     for (k=0; k < hi; k++)
     {
         for (l=0; l < wi; l++, m++)
         {
             *(*(arrbox + k) + l) = rand()*boxlength;
         }
     }
}

void buildvecbox()
{
     int m=0, k, l;
//	 srand(static_cast<unsigned int>(time(0)));
	 srand(1);

	for (k=0; k < hi; k++)
 	{
//	double *ptrvecbox = &vecbox[k][0];
        for (l=0; l < wi; l++, m++)
        {
			vecbox[k][l] = rand()*boxlength;
		}
	}
}

int main()
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	buildvecbox();
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	double vecduration = duration_cast<microseconds>( t2 - t1 ).count();


	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	buildarrbox();
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	double arrduration = duration_cast<microseconds>( t4 - t3 ).count();

//	int k, l;

	cout << "Duration vector matrix\n";
	cout << vecduration << " 10e-6 s\n";
	cout << "Duration array matrix\n";
	cout << arrduration << " 10e-6 s\n";
//	double *ptrvecbox = &vecbox[k][0]; 	  //use these to access vector by reference
//	cout <<  *(ptrvecbox + l) << "\n";  //no need to do that though?
	cout << "Sanity check \n";
	cout << arrbox[1][1] << "\n";
	cout << vecbox[1][1] << "\n";
	cout << arrbox[2][2] << "\n";
	cout << vecbox[2][2] << "\n";

	return 0;
}








