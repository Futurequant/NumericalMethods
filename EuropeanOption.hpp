// European Options Pricing Class 
// Header File
// 
// Matthew Morris 
// 
#include <string>
using namespace std;

class EuropeanOption
{
private: 
	void init();
	void copy(const EuropeanOption& obj);
	// "Kernel" Functions for options calculations
	// Const Functions are illegal to write inside 
	double CallPrice() const;
	double PutPrice() const;
	double CallDelta() const;
	double PutDelta() const;

public: 
	//Public member data for convience only
	double r;   	// Interest rate
	double sig; 	// Volatility
	double K;		// Strike Price
	double T; 		// Expiry Date
	double U; 		// Current underlying price
	double b; 		// Cost of carry 

	std::string optType; // Option name (call, put)

public: 
	// Constructors 
	EuropeanOption();		//Default Constructor
	EuropeanOption(const EuropeanOption& option2); //Copy Constructor
	EuropeanOption(const std::string& optionType);	   //Create option type

	// Destructor 
	virtual ~EuropeanOption();

	// Assignment operator 
	EuropeanOption& operator = (const EuropeanOption& option2);

	// Functions that calculate option price and some sensitivities 
	double Price() const;
	double Delta() const;

	// Modifier functions 
	void toggle();	//Change put to call (C/P P/C)
};

