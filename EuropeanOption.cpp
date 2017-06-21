////////////////////////////////////////////////////////////////////////////////
//European Option Pricing Class
// Implementation File
// 
// Matthew Morris 
// 
#include "EuropeanOption.hpp" // Declarations of functions
#include <cmath> 			  // For mathematical functions

// An approximation to the cumulative distribution function
// for the standard normal distribution
// Note: This is a recursive function
double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + 
    	k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

double EuropeanOption::CallPrice() const{ 

	double temp = sig * sqrt(T);

	double d1 = (log(U/K) + (b + (sig * sig) * 0.5) * T) / temp;
	double d2 = d1 - temp;

	return (U * exp(b-r) * T) * norm_cdf(d1) - (K * exp(-r * T)* norm_cdf(-d2));
}

double EuropeanOption::PutPrice() const{

	double temp = sig * sqrt(T);
	
	double d1 = (log(U/K) + (b + (sig * sig) * 0.5) * T)/temp;
	double d2 = d1 - temp;

	return (K * exp(-r * T) * norm_cdf(-d2)) - (U * exp(b - r) * T) * norm_cdf(-d1);
}

double EuropeanOption::CallDelta() const{

	double temp = sig * sqrt(T);

	double d1 = (log(U/K) + (b + (sig * sig) * 0.5) * T) /temp;
	
	return exp((b-r) * T) * norm_cdf(d1);
}

double EuropeanOption::PutDelta() const{

	double temp = sig * sqrt(T);

	double d1 = (log(U/K) + (b + (sig * sig) * 0.5) * T)/ temp;

	return exp((b-r) * T) * (norm_cdf(d1) - 1.0);
}

void EuropeanOption::init(){

//Initilise deafault values 
	//DEfault Values 
	r = 0.08;
	sig = 0.30;
	K = 65.0;
	T = 0.25;
	U = 65.0;
	b = r;
	optType = "C"; 
}
void EuropeanOption::copy(const EuropeanOption& obj){
	r = obj.r;
	sig = obj.sig;
	K = obj.K;
	T = obj.T;
	U = obj.U;
	b = obj.b;

	optType = obj.optType;
}

// Define Constrcutors
EuropeanOption::EuropeanOption(){
	// Defualt Call Option
	init();
}

EuropeanOption::EuropeanOption(const EuropeanOption& obj){
// Copy Constructor
	copy(obj);
}

EuropeanOption::EuropeanOption (const std::string& optionType){
	//Create Option Type
	init();
	optType = optionType;
	if (optType == "c")
			optType = "C";
}
// Destructor for memory clean up
EuropeanOption::~EuropeanOption(){
	// Destructor
}

EuropeanOption& EuropeanOption::operator = (const EuropeanOption& option2) {
	// Assignment Operator (deep copy)
		if (this == &option2) return *this;
		copy (option2);

		return *this;
}

// Functions that caclulate the option price and sensitivites 
double EuropeanOption::Price() const{

	if(optType == "C"){
		return CallPrice();
	}
	else
		return PutPrice();
}

double EuropeanOption::Delta() const{
	if (optType == "C")
		return CallDelta();
	else
		return PutDelta();
}

// Modifier Functions 
// Toggle between Put and Call
void EuropeanOption::toggle(){
	if(optType == "C")
		optType = "P";
	else
		optType = "C";
}