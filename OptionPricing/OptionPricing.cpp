// OptionPricing.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#pragma warning(disable: 4819)
#include<ql/quantlib.hpp>
#include <iostream>
#include "business248.hpp"
using namespace QuantLib;

int main()
{
    //Specify	1) calendar for specific market
	//			2) date count convention
	//			3) current date
	//			4) maturity date
	Calendar calendar = TARGET();
	DayCounter dayCounter = Business248();
	Date t0(29, Dec, 2011);
	Settings::instance().evaluationDate() = t0;
	Date T(26, Mar, 2012);

	//Specify	5) option type
	//			6) current stock price
	//			7) strike price
	Option::Type type(Option::Call);
	Real S0 = 23.44;
	Real K = 25;

	//Specify	8) dividend yield
	//			9) interest rate
	//			10) volatility
	Spread q = 0.00;
	Rate r = 0.00602;
	Volatility sigma = 0.3287;

	//Construct the stock quote object, the pointer and the handle of the object
	Handle<Quote> underlyingH(
		boost::shared_ptr<Quote>(
			new SimpleQuote(S0)
			)
	);

	//Construct the divident TS object, the pointer and the handle of the object
	Handle<YieldTermStructure> flatDividentTS(
		boost::shared_ptr<YieldTermStructure>(
			new FlatForward(t0, q, dayCounter)
			)
	);

	//Construct the interest rate TS object, the pointer and the handle of the object
	Handle<YieldTermStructure> flatTermStructure(
		boost::shared_ptr<YieldTermStructure>(
			new FlatForward(t0, r, dayCounter)
			)
	);

	//Construct the volatility TS object, the pointer and the handle of the object
	Handle<BlackVolTermStructure> flatVolTS(
		boost::shared_ptr<BlackVolTermStructure>(
			new BlackConstantVol(t0, calendar, sigma, dayCounter)
			)
	);

	//Construct the stock movement process object and the pointer of the object
	boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
		new BlackScholesMertonProcess(
			underlyingH, flatDividentTS, flatTermStructure, flatVolTS
		)
	);

	//Construct the option exercise object and the pointer of the object
	boost::shared_ptr<Exercise> europeanExercise(
		new EuropeanExercise(T)
	);

	//Construct the option payoff object and the pointer of the object
	boost::shared_ptr<StrikedTypePayoff>  payoff(
		new PlainVanillaPayoff(type, K)
	);

	//Construct the option object
	VanillaOption europeanOption(payoff, europeanExercise);

	//Construct the pricing engine and assign it to the pricing engine of the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new AnalyticEuropeanEngine(bsmProcess)
			)
	);

	//Output the result to the output window
	Real OptionPrice = europeanOption.NPV();
	std::cout << "S0 = " << S0 << std::endl;
	std::cout << "K = " << K << std::endl;
	std::cout << "sigma = " << sigma << std::endl;
	std::cout << "T = " << dayCounter.yearFraction(t0, T) << std::endl;
	std::cout << "r = " << r << std::endl;
	std::cout << "q = " << q << std::endl;
	std::cout << "C = " << OptionPrice << std::endl;
	system("pause");

	//Specify the jump intensity, jump mean and jump volatility objects
	double jumpIntensity = 0.4736842;
	double jumpVolatility = 0.4647414;
	double jumpMean = log(0.9989906) - 0.5*jumpVolatility*jumpVolatility;

	//Construct the jump intensity, jump volatility, and jump mean quote
	//objects and the pointer of the objects
	Handle<Quote> jmpIntH(
		boost::shared_ptr<Quote>(
			new SimpleQuote(jumpIntensity)
			)
	);

	Handle<Quote> jmpVolH(
		boost::shared_ptr<Quote>(
			new SimpleQuote(jumpVolatility)
			)
	);

	Handle<Quote> jmpMeanH(
		boost::shared_ptr<Quote>(
			new SimpleQuote(jumpMean)
			)
	);

	//Construct the stock movement process object and the pointer of the object
	boost::shared_ptr<Merton76Process> jumpProcess(
		new Merton76Process(
			underlyingH, flatDividentTS, flatTermStructure, flatVolTS,
			jmpIntH, jmpMeanH, jmpVolH)
	);

	//Construct the pricing engine object and assign it to the pricing engine 
	//for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new JumpDiffusionEngine(jumpProcess)
			)
	);

	//Output the results to the output window
	OptionPrice = europeanOption.NPV();
	std::cout << "C (with jump) = " << OptionPrice << std::endl;
	system("pause");

}