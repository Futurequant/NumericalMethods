// CW2 Task 2 and Task 3 
// Matthew Morris 
// Last Updated: 27/3/2017
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

// Define Normal Distribution using Polar-Marsaglia 
double normal(double mean, double std)
{
	static int iset = 0;
	static double gset; 
	double fac, r, v1, v2;
	// Create two normally-distributed numbers 
	if (iset == 0)
	{
		r = 0;
		// Compute two possibilities 
		do 
			{
				v1 = 2.0 * rand() / RAND_MAX - 1.0;
				v2 = 2.0 * rand() / RAND_MAX - 1.0;
				// define the radius 
				r = v1 * v1 + v2 * v2;
			}
			while (r >= 1.0 || r ==0.0);
			
			// Inside the unit circle? if not try again..
			fac = sqrt((-2 * log(r)) / r); // Box-Muller Transform
			gset = (v1 * fac);
			iset = 1; // save one
			v2 = v2 * fac * std + mean; // Scale and return on
			return v2;
	}
	else 
	{
		iset = 0;
		return (gset * std) + mean;
			// scale and return the saved one
	return 0;
	}
}
// Estimate the Guassian Integral
void gauss_integral(vector<int>& nums, vector<double>& est)
{
	for (auto i=nums.begin(); i != nums.end(); i++)
	{
		double temp = 0.0;
		double sum = 0.0;
		for(auto j=0; j<=*i; j++)
		{
			double u = 0.0;
			u = normal(0,1);
			sum += pow(u,2); 
		}
		temp = sum / *i;
		est.push_back(temp);
	}
}
// Estimate Monte Carlo Integral
void monte_carlo(vector<int>& nums, vector<double>& mcEst)
{
	for(auto i=nums.begin(); i != nums.end(); i++)
	{
		double sum = 0.0;
		double result = 0.0;
		double brac = 0.0;
		double nom = 0.0;
		double den = 0.0;
		double power = 0.0;

		for(auto j=0; j<=*i; j++)
		{
			double u = 0.0;
			u = double(rand())/RAND_MAX;
			brac = ((1/u) - 1 );
			power = -1 * (brac * brac);
			nom = exp(power);
			den = pow(u,2);
			sum += nom / den;
		}
		result = sum / *i;
		mcEst.push_back(result);
	}
}

int main(){
	// House Keeping 
	ofstream print;
	print.open("mc.csv");

	// Define Vectors
	vector<int> numbers = {10, 100, 1000, 10000, 100000, 1000000, 10000000};
	vector<double> est;
	vector<double> mcEst;

	// Call Functions
	gauss_integral(numbers, est);
	monte_carlo(numbers, mcEst);

	// Print Difference in Errors between vectors to the terminal
	cout << "********************************************" << endl;
	cout << "Task 2 Normal Random Numbers" << endl;
	auto numit = numbers.begin();	
	for(auto i : est)
	{
		// Print out comma separated value file 
		//print << *numit << "," << i << endl;
		// Terminal print 
		cout << "N " << *numit << ": " << i << endl;
		numit++;
	}
	// Print Monte Carlo estimate results 
	cout << "********************************************" << endl;
	cout << "Task 3 Uniform Random Numbers" << endl;
	numit = numbers.begin();
	for(auto i : mcEst)
	{	
		// Print out comma separated value file 
		print << *numit << "," << i << endl;
		// Terminal print 
		cout << "N " << *numit << ": " << i << endl;
		numit++;
	}

}

