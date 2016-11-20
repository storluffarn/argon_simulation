#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector <double> data = {1,10,100};

double mean = accumulate(data.begin(),data.end(),0.0)/data.size();

double corrfunc()
{
	//auto corr = accumulate(data.begin(),data.end(),0.0,[&](int k, int l){return (k - mean)*(l-mean);})/data.size();

	double corr = 0;

	for ( unsigned int k = 0; k < data.size() - 1; k++)
	{
		corr += (data[k] - mean)*(data[k+1] - mean);
	}

	return corr/data.size();
}

int main()
{
	double output = corrfunc();

	cout << output << endl;

	double sanity = ((1-mean)*(10-mean)+(10-mean)*(100-mean))/3;

	cout << sanity << " " << mean << endl;
}
