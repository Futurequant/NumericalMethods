// Asian Option Pricing Calculator 
// Matthew Morris 
// Last Updated: 27/3/2017
//
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <fstream>
using namespace std;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Function Delcarations 
double CDF(double X);

double normal(double mean, double std);

void calc_path(vector<double>& spot_prices, const double& r, const double& v, 
			   const double& T);

// Averaging Functions 
vector<double> simple_moving_avg(vector<double>& prices, int window);

double arithmetic_avg(vector<double>& spot_prices);

double geometric_avg(vector<double>& spot_prices);

// Pricing Functions
double discrete_asian_call(const int num_paths, const int num_intervals, 
						   const int window, const double K, const double S, 
						   const double r, const double v, const double T);

double discrete_asian_put(const int num_paths, const int num_intervals, 
						   const int window, const double K, const double S, 
						   const double r, const double v, const double T);

double arithmetic_asian_call(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T);

double arithmetic_asian_put(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T);

double geometric_asian_call(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T);

double geometric_asian_put(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Main Program
int main()
{
	// Housekeeping 
	ofstream print;
	print.open("sma.csv");
	clock_t timet;
	timet = clock();
	srand((unsigned)time(NULL));

	// Define Params 
	int num_paths = 100000;
	int num_intervals = 252;
	int window = 20;
	double K = 100;
	double S = 100;
	double r = 0.05;
	double v = 0.2;
	double T = 1.0;

	// Test Pricing functions
	double arith_call_price = arithmetic_asian_call(num_paths, num_intervals, 
													K, S, r, v, T);
	double arith_put_price = arithmetic_asian_put(num_paths, num_intervals, 
													K, S, r, v, T);
	double geo_call_price = geometric_asian_call(num_paths, num_intervals, 
													K, S, r, v, T);
	double geo_put_price = geometric_asian_put(num_paths, num_intervals, 
													K, S, r, v, T);

	double dis_call_price = discrete_asian_call(num_paths, num_intervals, window, 
													K, S, r, v, T);

	double dis_put_price = discrete_asian_put(num_paths, num_intervals, window, 
													K, S, r, v, T);

	// Output to terminal
	cout << "Arithmetic Asian Call: "<< arith_call_price << endl;
	cout << "Arithmetic Asian Put: "<< arith_put_price << endl;
	cout << "Geometric Asian Call: "<< geo_call_price << endl;
	cout << "Geometric Asian Put: "<< geo_put_price << endl;
	cout << "Discrete Asian Call: " << dis_call_price << endl;
	cout << "Discrete Asian Put: " << dis_put_price << endl;

	// Test Simple Moving Average Function
	vector<double> spot_prices(num_intervals, S);
	vector<double> sma_prices(num_intervals, S);
	// Simulate asset path
	calc_path(spot_prices, r, v, T);
	sma_prices = simple_moving_avg(spot_prices, window);
	for (int i = 0; i < spot_prices.size(); i++) 
	{
		// Print path and SMA to csv file
		print << sma_prices[i] << ",";
		print << spot_prices[i] << ",";
		print << " " << endl;
	}
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Function Definitions
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

// Calculate Asset Path Function 
void calc_path(vector<double>& spot_prices, const double& r, const double& v, 
			   const double& T)
{
	// Define parameters for stochastic process
	double dt = T/static_cast<double>(spot_prices.size());
	double drift = (1 + r * dt);
	double diff = v * sqrt(dt);
	double milstein = 0.5 * v * v * dt;

	// Simulate asset path using Milstein method 
	for(int i=1; i<spot_prices.size(); i++){
		double random = normal(0,1);
		spot_prices[i] = spot_prices[i-1] * 
						 (drift + (random * diff) + (milstein * 
						 ((random * random) - 1)));
	}
}

////////////////////////////////////////////////////////////////////////////////
// Averaging Functions 
double arithmetic_avg(vector<double>& spot_prices) 
{
	unsigned length = spot_prices.size();
	double sum = accumulate(spot_prices.begin(), spot_prices.end(), 0);
  	double arith_mean = sum / static_cast<double>(length);

  	return arith_mean;
}

double geometric_avg(vector<double>& spot_prices)
{
	unsigned length = spot_prices.size();
	double log_sum = 0.0;
		for (int i=0; i<length; i++) {
			log_sum += log(spot_prices[i]);
		}
	double geom_mean = exp(log_sum / static_cast<double>(length) );
	return geom_mean;
}

vector<double> discrete_avg(vector<double>& prices, int window)
{
	//ofstream print;
	//print.open("Paths.csv")
	int cnt = 0;
	double sum = 0.0;
	double temp = 100.0;

	unsigned length = prices.size();
	vector<double> discrete_avg(length, prices[0]); 
	//
	for (int i = 0; i < prices.size(); i++) 
	{
		discrete_avg[i] = temp;
		//print << discrete_avg[i] << ",";
		//print << prices[i] << ",";
    	sum += prices[i];
    	cnt++;

    if (cnt >= window) 
    {
        temp = (sum / static_cast<double>(window));
        discrete_avg[i] = temp;
        sum = 0.0;
        cnt = 0.0;
    }
    }
	return discrete_avg;
}

vector<double> simple_moving_avg(vector<double>& prices, int window)
{
	int cnt = 0.0;
	double sum = 0.0;
	unsigned length = prices.size();
	vector<double> sma(length, prices[0]); 
	//
	for (int i = 0; i < length; i++) {
    sum += prices[i];
    cnt++;
    if (cnt >= window) {
        sma[i] = (sum / static_cast<double>(window));
        sum -= prices[cnt - window];           
        }
    }
    return sma;
}
////////////////////////////////////////////////////////////////////////////////
// Option Pricing Functions
double discrete_asian_call(const int num_paths, const int num_intervals, 
						   const int window, const double K, const double S, 
						   const double r, const double v, const double T)
{
	// Define pricing parameters
	double payoff = 0.0;
	double spot_price = 0.0;
	vector<double> avg_prices;
	vector<double> spot_prices(num_intervals, S);
	// Calculate payoff for each simulated asset path
	for (int i=0; i<num_paths; i++) {

		// calculate the path and discrete averages
    	calc_path(spot_prices, r, v, T);
    	avg_prices = discrete_avg(spot_prices, window);

    	// Get spot price from previous average window
   		int length = avg_prices.size();
   		spot_price = avg_prices[length-1];

   		// Discount Payoff
    	payoff += max(spot_price - K, 0.0);
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }
    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

double discrete_asian_put(const int num_paths, const int num_intervals, 
						   const int window, const double K, const double S, 
						   const double r, const double v, const double T)
{
	double payoff = 0.0;
	double spot_price = 0.0;
	int length = 0.0;
	vector<double> avg_prices;

	for (int i=0; i<num_paths; i++) {

		// calculate the path and discrete averages
		vector<double> spot_prices(num_intervals, S); 
    	calc_path(spot_prices, r, v, T);
    	avg_prices = discrete_avg(spot_prices, window);

    	// Get spot price from previous average window
   		length = avg_prices.size();
   		spot_price = avg_prices[length-1];

   		// Discount Payoff
    	payoff += max(K - spot_price, 0.0);
    	// Reset vector for next path
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }
    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

double arithmetic_asian_call(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T)
{
	vector<double> spot_prices(num_intervals, S); 
	double payoff = 0.0;
	int window = 100;
	for (int i=0; i<num_paths; i++){
    	calc_path(spot_prices, r, v, T);
    	payoff += max(arithmetic_avg(spot_prices) - K, 0.0);
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }

    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

double arithmetic_asian_put(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T)
{
	vector<double> spot_prices(num_intervals, S); 
	double payoff = 0.0;
	for (int i=0; i<num_paths; i++) {
    	calc_path(spot_prices, r, v, T);
    	payoff += max(K - arithmetic_avg(spot_prices), 0.0);
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }
    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

double geometric_asian_call(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T)
{
	vector<double> spot_prices(num_intervals, S); 
	double payoff = 0.0;
	for (int i=0; i<num_paths; i++) {
    	calc_path(spot_prices, r, v, T);
    	payoff += max(geometric_avg(spot_prices) - K, 0.0);
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }
    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

double geometric_asian_put(const int num_paths, const int num_intervals,  
							 const double K, const double S, const double r, 
						     const double v, const double T)
{
	vector<double> spot_prices(num_intervals, S); 
	double payoff = 0.0;
	for (int i=0; i<num_paths; i++) {
    	calc_path(spot_prices, r, v, T);
    	payoff += max(K - geometric_avg(spot_prices), 0.0);
    	fill(spot_prices.begin(), spot_prices.end(), S);
    }
    return(payoff / static_cast<double>(num_paths)) * exp(-r*T); 
}

