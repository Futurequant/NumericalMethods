// European Option Pricing Via Monte Carlo 
// 
// Matthew Morris 
// 
#include <algorithm>
#include <iostream>
#include <cmath>

// Implementation of the Box-Muller method, used to generate Guassian 
// random numbers. In C++11 we can use std::normal_distribution form 
// the <random> library 

double gaussian_box_muller(){
	double x = 0.0;
	double y = 0.0; 
	double euclid_sq;
	
	// Continute generating two uniform random variable until
	// the square of their "Euclidean Distance" is less than 
	// unity 
	do {
		x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
		euclid_sq = x*x + y*y;
	}while (euclid_sq >= 1.0);

	return x * sqrt(-2 * log(euclid_sq)/euclid_sq);
}

// Pricing a European Option using the Monte Carlo Method 

double monte_carlo_call_price(const &int nums_sims, const double& S,
	const double& K, const double& r, const double& v, const double& T){
	
	double S_adjust = S * exp(T * (r-0.5*v*v));
	double S_cur = 0.0;
	double payoff_sum = 0.0;

	for (int i=0; i<nums_sims; i++){
		double gauss_bm = gaussian_box_muller();
		S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
		payoff_sum += std::max(S_cur - K, 0.0);
	}

	return (payoff_sum / static_cast<double>(nums_sims)) * exp(-r*T);
}

double monte_carlo_put_price(const &int nums_sims, const double& S,
	const double& K, const double& r, const double& v, const double& T){

	double S_adjust = S * exp(T * (r-0.5*v*v));
	double S_cur = 0.0;
	double payoff_sum = 0.0;

	for (int i=0; i<nums_sims; i++){
		double gauss_bm = gaussian_box_muller();
		S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
		payoff_sum += std::max(K - S_cur, 0.0);
	}
}

int main(){

	//Option Parameter List 
	int num_sims = 1000000;
	double S = 100.0;
	double K = 100.0;
	double r = 0.05;
	double v = 0.2;
	double T = 1.0;

	// Calculate the Call Put Values via Monte Carlo 
	double monte_carlo_call_price(nums_sims, S, K, r, v, T);
	double monte_carlo_put_price(nums_sims, S, K, r, v, T);

	std::cout << "Number of Paths: " << num_sims << std::endl;
  	std::cout << "Underlying:      " << S << std::endl;
  	std::cout << "Strike:          " << K << std::endl;
  	std::cout << "Risk-Free Rate:  " << r << std::endl;
  	std::cout << "Volatility:      " << v << std::endl;
  	std::cout << "Maturity:        " << T << std::endl;

  	std::cout << "Call Price:      " << call << std::endl;
  	std::cout << "Put Price:       " << put << std::endl;

  	return 0;	
}



