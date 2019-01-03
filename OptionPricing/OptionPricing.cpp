// OptionPricing.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#pragma warning(disable: 4819)
#include<ql/quantlib.hpp>
#include <iostream>
#include "business248.hpp"
using namespace QuantLib;

class BSMJProcess : public BlackScholesMertonProcess {

public:
	//The contructor for initializing the parameters
	BSMJProcess(
		const Handle<Quote>& stateVariable,
		const Handle<YieldTermStructure>& dividendTS,
		const Handle<YieldTermStructure>& riskFreeTS,
		const Handle<BlackVolTermStructure>& blackVolTS,
		//allow the main program to pass the jump paramters
		const double jumpInt,
		const double logJMean,
		const double logJVol,
		const double kappa,
		//Initilize the random number generator by the system clock
		const MersenneTwisterUniformRng& URng_Pt = MersenneTwisterUniformRng(SeedGenerator::instance().get()),
		const MersenneTwisterUniformRng& URng_Ji = MersenneTwisterUniformRng(SeedGenerator::instance().get()),
		const boost::shared_ptr<discretization>& d =
		boost::shared_ptr<discretization>(new EulerDiscretization)
	)
		: BlackScholesMertonProcess(stateVariable, dividendTS, riskFreeTS, blackVolTS, d),
		//Initialize these parameters accordingly
		jumpIntensity_(jumpInt), logMeanJump_(logJMean), logJumpVolatility_(logJVol),
		kappa_(kappa), URng_Pt_(URng_Pt), URng_Ji_(URng_Ji) {}

	//Implement the adjustment for the drift term
	double drift_adj(Time dt) const {
		return jumpIntensity_ * kappa_*dt;
	}

	//Implement the adjustment for the diffusion term
	double diffusion_adj(Time dt) const {
		InverseCumulativePoisson invP(jumpIntensity_*dt);
		InverseCumulativeNormal invN(logMeanJump_, logJumpVolatility_);

		int Pt = (int)invP(URng_Pt_.next().value);

		double total_jump = 0;
		double Ji;

		for (int i = 1; i <= Pt; i++) {
			Ji = invN(URng_Ji_.next().value);
			total_jump += Ji;
		}

		return total_jump;
	}

	//Override the original evovle function by the modified func with jump
	Real evolve(Time t0, Real x0, Time dt, Real dw) const {
		return apply(x0, discretization_->drift(*this, t0, x0, dt) - drift_adj(dt) +
			stdDeviation(t0, x0, dt)*dw + diffusion_adj(dt));
	}

private:
	//Additional parameters for the jump diffussion modelling
	double jumpIntensity_;
	double logMeanJump_;
	double logJumpVolatility_;
	double kappa_;
	//Uniform random number generateors for the Monte Carlo Simulation
	MersenneTwisterUniformRng URng_Pt_;
	MersenneTwisterUniformRng URng_Ji_;
};

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


	//Specify the jump intensity, jump mean and jump volatility objects
	double jumpIntensity = 0.4736842;
	double jumpVolatility = 0.4647414;
	double jumpMean = log(0.9989906) - 0.5*jumpVolatility*jumpVolatility;
	double kappa = 0.9989906 - 1; //The average jump size - 1

	//Construct the stock movement process object and the pointer to the object
	boost::shared_ptr<BSMJProcess> bsmjProcess(
		new BSMJProcess(
			underlyingH, flatDividentTS, flatTermStructure, flatVolTS,
			jumpIntensity, jumpMean, jumpVolatility, kappa
		)
	);

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

	// PseudoRandom
	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<PseudoRandom>(
				bsmProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				10000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);
	
	//write column headings
	Size widths[] = { 25, 25, 25 };
	std::cout << std::setw(widths[0]) << std::left << "MC Samples"
		<< std::setw(widths[1]) << std::left << "Call price estimate"
		<< std::setw(widths[2]) << std::left << "Err of estimate"
		<< std::endl;


	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 10000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << europeanOption.errorEstimate()
		<< std::endl;
	
	
	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<PseudoRandom>(
				bsmProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				100000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);

	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 100000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << europeanOption.errorEstimate()
		<< std::endl;

	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<PseudoRandom>(
				bsmProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				1000000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);

	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 1000000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << europeanOption.errorEstimate()
		<< std::endl;
	

	//Jump diffusion
	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<LowDiscrepancy>(
				bsmjProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				10000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);

	//write column headings
	std::cout << std::setw(widths[0]) << std::left << "MC Samples"
		<< std::setw(widths[1]) << std::left << "Call price estimate"
		<< std::endl;

	
	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 10000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::endl;
	
	
	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<LowDiscrepancy>(
				bsmjProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				100000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);
	
	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 100000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::endl;

	//Contruct the pricing engin object and assign it to the pricing engine for the option
	europeanOption.setPricingEngine(
		boost::shared_ptr<PricingEngine>(
			new MCEuropeanEngine<LowDiscrepancy>(
				bsmjProcess, //process
				10, //timeSteps
				Null<Size>(),//timeStepsPerYear
				false, //brownianBridge
				false, //antitheticVariate
				1000000, //requriedSamples
				Null<Real>(), //requiredTolerance
				Null<Size>(), //maxSamples
				SeedGenerator::instance().get() //seed
				)
			)
	);

	//Output the result to the output window
	std::cout << std::setw(widths[0]) << std::left << 1000000
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::endl;
	
}