// EUROPEAN OPTION PRICER 
// MATTHEW MORRIS 
// LAST UPDATED: 3/3/17
//
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// Define CDF Function 
double CDF(double X)
{
	const double pi=4.0*atan(1.0);
	const double a1=0.319381530, a2=-0.356563782; 
	const double a3=1.781477937, a4=-1.821255978, a5=1.330274429;

	double N, CDF, n;
	double x=0, k=0; 
	x=fabs(X);
	k=1/(1+0.2316419*x);
	n=(1/sqrt(2*pi))*exp(-0.5*x*x); 
	N=1-n*(a1*k+a2*k*k+a3*pow(k,3)+a4*pow(k,4)+a5*pow(k,5)); 

	if (X>=0) 
		CDF = N;
	else 
		CDF = 1-N;

	return CDF;
}

double polar_mars(double std, double mean)
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
			v2 = v2 * fac * std + mean; // Scale and return one

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
// Calculate the price of a Call Option using the closed form solution
double analytical_call_price(const double& K, const double& S, const double& r,
					  const double& sigma, const double& T)
{
	double sigma_sqrt_T = sigma * sqrt(T);
	double d_1 = (log(S/K) + (r + sigma * sigma * 0.5) * T) / sigma_sqrt_T;
	double d_2 = d_1 - sigma_sqrt_T;

	return S * CDF(d_1) - K * exp(-r*T) * CDF(d_2);
}
// Calculate the price of a Put Option using the closed form solution
double analytical_put_price(const double& K, const double& S, const double& r,
					  const double& sigma, const double& T)
{
	double sigma_sqrt_T = sigma * sqrt(T);
	double d_1 = (log(S/K) + (r + sigma * sigma * 0.5) * T) / sigma_sqrt_T;
	double d_2 = d_1 - sigma_sqrt_T;
	return K * exp(-r*T) * CDF(-d_2) - S * CDF(-d_1);
}

// Pricing a European call option with a Monte Carlo method
double monte_carlo_call_price(const int& num_sims, const double& K,
							  const double& S, const double& r, 
							  const double& sigma, const double& T) 
{
	// Define variables for computation
	double payoff = 0.0;
	// Define the increment in the path
	double S_inc = S * exp( T * (r - (0.5 * sigma * sigma)));
	// Current price of the path 
 	double S_cur = 0.0;

	for (int i=0; i<num_sims; i++) 
	{
    	double random = polar_mars(1.0, 0.0);
    	S_cur = S_inc * exp(sqrt(sigma * sigma * T) * random);
    	payoff += max(S_cur - K, 0.0);
	}	

	// make use of static_cast to cast const in to double for calc
	return (payoff / static_cast<double>(num_sims)) * exp(-r*T);
}

// Pricing a European put option with a Monte Carlo method
double monte_carlo_put_price(const int& num_sims, const double& K,
							  const double& S, const double& r, 
							  const double& sigma, const double& T) 
{
	double payoff = 0.0;
	double S_inc = S * exp( T * (r - 0.5 * sigma * sigma));
 	double S_cur = 0.0;

	for (int i=0; i<num_sims; i++) 
	{
    	double random = polar_mars(1.0, 0.0);
    	S_cur = S_inc * exp(sqrt(sigma * sigma * T) * random);
    	payoff += max(K - S_cur, 0.0);
	}

	// make use of static_cast to cast const in to double for calc
	return (payoff / static_cast<double>(num_sims)) * exp(-r*T);
}

int main() 
{
	// Housekeeping 
	clock_t timet;
	timet = clock();
	srand((unsigned)time(NULL));
	ofstream print;
	print.open("option.csv");

	// Define Params 
	int num_sims = 10000000;
	double K = 100.0;
	double S = 100.0;
	double r = 0.05;
	double sigma = 0.2;
	double T = 1.0;

	// Get prices 
	double call_price = monte_carlo_call_price(num_sims, K, S, r, sigma, T);
	double put_price = monte_carlo_put_price(num_sims, K, S, r, sigma, T);
	double closed_call_price = analytical_call_price(K, S, r, sigma, T);
	double closed_put_price = analytical_put_price(K, S, r, sigma, T);

	// Terminal Ouput for user 
	cout << "\nNumber of Simulations: " << num_sims << endl;
	cout << "Strike Price: " << K << endl;
	cout << "Underlying Price: " << S << endl;
	cout << "Risk free rate: " << r << endl;
	cout << "Time to Maturity: " << T << endl;
	cout << "Call Price: " << call_price << endl;
	cout << "Put Price: " << put_price << endl;
	cout << "Closed Form Call Price: " << closed_call_price << endl;
	cout << "Closed Form Put Price: " << closed_put_price << endl;
	cout << "********************************" << endl;

	vector<int> numbers = {1000, 10000, 100000, 1000000, 10000000};

	// Demonstrate the Law of Large numbers through increased simulations 
	for (auto i = numbers.begin(); i != numbers.end(); i++)
	{	
		double est_call = monte_carlo_call_price(*i, K, S, r, sigma, T);
		double difference = closed_call_price - est_call;
		cout << "N: " << *i << " " << est_call << " " << difference << endl;
		print << *i << "," << difference << endl;
	}
	cout << "********************************" << endl;
	timet = clock() - timet;
	cout << "Program Run Time: " << ((float)timet) /CLOCKS_PER_SEC << " second(s)." << endl;
	cout << "*************************************" << endl;

	return 0;

}
