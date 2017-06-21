////////////////////////////////////////////////////////////////////////////////
// European Option Main Pricing Function 
// 
// Matthew Morris 
// 
#include "EuropeanOption.hpp"
#include <iostream>
int main(){
	// All Options are European 

	// Call option on a stock 
	EuropeanOption callOption;
	cout << "Put option on on stock: " << callOption.Price() << endl;

	// Call option on a stock index 
	EuropeanOption indexOption;
	indexOption.optType = "P";
	indexOption.U = 100.0;
	indexOption.K = 95.0;
	indexOption.T = 0.5;
	indexOption.r = 0.10;
	indexOption.sig = 0.2;

	double q = 0.05; // Dividend Yeild

	cout << "Put option on on index: " << indexOption.Price() << endl;

	// Call and put options on a future
	EuropeanOption futureOption;

	futureOption.optType = "P";
	futureOption.U = 19.0;
	futureOption.K = 19.0;
	futureOption.T = 0.75;
	futureOption.r = 0.10;
	futureOption.sig = 0.28;
	futureOption.b = 0.0;

	cout << "Put option on on future: " << futureOption.Price() << endl;


	// Now change over to a call on the option
	futureOption.toggle();
	cout << "Call on future: " << futureOption.Price() << endl;

}