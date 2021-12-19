/*
** Author: Wu Chen
** 2017/08/07
*/

#include "PricingRoutine.h"

double GetQCTrigger(
	double S0,                            /* Current Stcok Prcie               */
	double CTrigger,                      /* Callable Trigger Value            */
	double T,                             /* Time to Maturity                  */
	double r,                             /* Constant interest rate            */
	double sigma,                         /* Stock Volatility                  */
	int n,                                /* n out of m days                   */
	int m,                                /* n out of m days                   */
	int N,                                /* Simulation Times                  */
	std::vector<double> S_history,        /* Historical Stock Price            */
	std::vector<double> CTrigger_history, /* Historical CTrigger               */
	std::vector<double> noCallTime,       /* No allow call time period         */
	std::string windcode,                 /* windcode for current bond         */
	std::string date                      /* calculation date for current bond */
) {
	std::wstring FilePath = L".\\logs_GetQCTrigger\\";
	CreateDirectoryW(FilePath.c_str(), NULL);
	std::ofstream out(".\\logs_GetQCTrigger\\" + windcode + " " + date + ".txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	std::cout << "**************************************************" << std::endl;
	std::cout << "***         start GetQCTrigger(...)            ***" << std::endl;
	std::cout << "**************************************************" << std::endl;

	return calculateQCTrigger(S0, CTrigger, T, r, sigma, n, m, N, S_history, CTrigger_history, noCallTime);
}

double GetQPTrigger(
	double S0,                            /* Current Stcok Prcie               */
	double PTrigger,                      /* Callable Trigger Value            */
	double T,                             /* Time to Maturity                  */
	double r,                             /* Constant interest rate            */
	double sigma,                         /* Stock Volatility                  */
	int n,                                /* n out of m days                   */
	int m,                                /* n out of m days                   */
	int N,                                /* Simulation Times                  */
	std::vector<double> S_history,        /* Historical Stock Price            */
	std::vector<double> PTrigger_history, /* Historical PTrigger               */
	std::vector<double> noPutTime,        /* No allow put time period          */
	std::string windcode,                 /* windcode for current bond         */
	std::string date                      /* calculation date for current bond */
) {
	std::wstring FilePath = L".\\logs_GetQPTrigger\\";
	CreateDirectoryW(FilePath.c_str(), NULL);
	std::ofstream out(".\\logs_GetQPTrigger\\" + windcode + " " + date + ".txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	std::cout << "**************************************************" << std::endl;
	std::cout << "***         start GetQPTrigger(...)            ***" << std::endl;
	std::cout << "**************************************************" << std::endl;

	return calculateQPTrigger(S0, PTrigger, T, r, sigma, n, m, N, S_history, PTrigger_history, noPutTime);
}

GetImpSpreadOutput GetImpSpread(
	double faceValue,                    /* Convertible bond face value              */
	double conversionPrice,              /* Conversion price                         */
	double bonusRate,                    /* Bonus Rate for holding to expiry         */
	double S0,                           /* Current stock price                      */
	double sigma,                        /* Historical Volatility                    */
	double T,                            /* Time to maturity                         */
	double q,                            /* Stock dividend yield                     */
	double p,                            /* Default Probability                      */
	double R,                            /* Recovery rate                            */
	double eta,                          /* When default, stock price = (1-eta) * S  */
	double QCTrigger,                    /* Quivalent Callable Trigger               */
	double QPTrigger,                    /* Quivalent Puttable Trigger               */
	double theta,                        /* Scheme Parameter, theta = 0, 0.5, 1 for  */
										 /* explicit, CN and implicit respectively   */
	double r,                            /* Risk free rate                           */
	double annuFactor,                   /* Annualized Factor                        */
	double callRedeemPrice,              /* Call option redemption price             */
	double putResellPrice,               /* Put option reselling price               */
	double mktPrice,                     /* Market Price                             */
	std::vector<double> noCallTime,      /* No allow call time period                */
	std::vector<double> noPutTime,       /* No allow put time period                 */
	std::vector<double> noConvertTime,   /* No allow conversion time period          */
	std::vector<double> couponTimeRate,  /* Coupon rate for different times          */
	std::string windcode,                /* windcode for current bond                */
	std::string date                     /* calculation date for current bond        */
) {
	std::wstring FilePath = L".\\logs_GetImpSpread\\";
	CreateDirectoryW(FilePath.c_str(), NULL);
	std::ofstream out(".\\logs_GetImpSpread\\" + windcode + " " + date + ".txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	GetImpSpreadOutput result;
	std::cout << "**************************************************" << std::endl;
	std::cout << "***         start GetImpSpread(...)            ***" << std::endl;
	std::cout << "**************************************************" << std::endl;
	// accuracy level
	int Nt = ceil(T * 252.0), Ns = 2000;
	double ds = QCTrigger / Ns;
	double dt = T / Nt;
	double tol = _MIN(ds * ds, pow(dt, 1.5));
	std::cout << "   Current tol = " << tol << std::endl;

	std::cout.precision(18);
	std::cout << "   Define objective function ..." << std::endl;
	// objective function
	auto f0 = [Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, 
		eta, QCTrigger, QPTrigger, theta, r, annuFactor, noCallTime, callRedeemPrice, putResellPrice,
		noPutTime, noConvertTime, couponTimeRate, mktPrice](double oas) 
	{return ConvertibleBondPDEPrice(Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma,
		T, q, p, R, eta, QCTrigger, QPTrigger, theta, r, oas, annuFactor, callRedeemPrice, putResellPrice,
		noCallTime, noPutTime, noConvertTime, couponTimeRate) - mktPrice; };
	std::cout << "   Define objective function success ..." << std::endl;

	double impSpread;
	std::cout << "   Start root bracketing ..." << std::endl;
	// initial interval
	double a = -1, b = 1;
	std::cout << "   Initial interval: (" << a << "," << b << ")" << std::endl;

	if (RootBracketing(f0, a, b, -1e-16, 1e16)) {
		std::cout << "   Root bracketing successed!" << std::endl;
		impSpread = rfbrent(f0, a, b, tol);
		if (ALMOST_EQUAL(mktPrice, S0 * faceValue / conversionPrice)) {
			std::cout << "   Market Price = nS, start minimum root bisection searching ..." << std::endl;
			double lower = impSpread / 2;
			double upper = impSpread;
			double minimum_impSpread = impSpread;
			double accuracy_level = 1e-6;
			int count = 0;
			double tested_max_lower = a;
			while(upper - lower >= accuracy_level) {
				count++;
				if (f0(lower) < tol && f0(lower) > - tol) {
					minimum_impSpread = lower;
					upper = lower;
					lower = lower / 2;
					if (lower < tested_max_lower) {
						lower = (upper + tested_max_lower) / 2;
					}
				}
				else {
					if (lower > tested_max_lower) {
						tested_max_lower = lower;
					}
					lower = (upper + lower) / 2;
				}
				std::cout << "   Current (lower, upper) = (" << lower << "," << upper << ")" << std::endl;
				std::cout << "   Current minimum_impSpread = " << minimum_impSpread << std::endl;
			}
			std::cout << "   Total loop times = " << count << std::endl;
			std::cout << "   Finding impSpread successed! impSpread = " << minimum_impSpread << std::endl;
			result.impliedSpread = minimum_impSpread;
			result.indicator = true;
			std::cout.rdbuf(coutbuf); //reset to standard output again
			return result;
		}
		else {
			std::cout << "   Finding impSpread successed! impSpread = " << impSpread << std::endl;
			result.impliedSpread = impSpread;
			result.indicator = true;
			std::cout.rdbuf(coutbuf); //reset to standard output again
			return result;
		}
	}
	else {
		std::cout << "   Cannot Find the implied spread, will return -1 !!" << std::endl;
		result.impliedSpread = -1;
		result.indicator = false;
		std::cout.rdbuf(coutbuf); //reset to standard output again
		return result;
	}
}

int test_routine()
{
    /*
	imp spread
	*/
	/*
    double faceValue = 100.0, conversionPrice = 43.28, bonusRate = 0.015, S0 = 49.56, sigma = 0.8547012569581931, T_ACT = 2.9397260273972603;
	double q = 0, p = 0, R = 1, eta = 0;
	double QCTrigger = 10000000.0, QPTrigger = 0.0;
	int theta = 1;
	double annuFactor = 365.25;
	double r = 0.0343161763731807;

	double value;
	double noCallTime1[2] = { 0, 0 };
	double noPutTime1[2] = { 0, 0 };
	double noConvertTime1[2] = { 1.9917808219178084, 3.0027397260273974 };
	double couponTimeRate1[3] = { 0.015, 0.015, 0.015 };
	std::vector<double> noCallTime(noCallTime1, noCallTime1 + 2);
	std::vector<double> noPutTime(noPutTime1, noPutTime1 + 2);
	std::vector<double> noConvertTime(noConvertTime1, noConvertTime1 + 2);
	std::vector<double> couponTimeRate(couponTimeRate1, couponTimeRate1 + 3);
	double callRedeemPrice = -1;
	double putResellPrice = -1;

	double mktPrice = 133.38;

	clock_t start;
	double duration;
	start = clock();

	GetImpSpreadOutput result;
	result = GetImpSpread(faceValue, conversionPrice, bonusRate, S0, sigma, T_ACT, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, annuFactor, callRedeemPrice, putResellPrice, mktPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate, windcode, date);

	if (result.indicator) {
		std::cout << result.impliedSpread << std::endl;
	}
	else
	{
		std::cout << "     other treatment should be done in python." << std::endl;
	}

	std::cout.precision(17);
	//std::cout << value << std::endl;
	
	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "The running time is " << duration << std::endl;
	*/
	return 0;
}

int main_test2() {
	test_routine();
	system("pause");
	return 0.0;
}