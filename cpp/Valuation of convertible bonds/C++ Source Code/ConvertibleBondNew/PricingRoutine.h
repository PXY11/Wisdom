/*
** Author: Wu Chen
** 2017/07/07
*/
#pragma once
#include "CalculateQCTrigger.h"
#include "CalculateQPTrigger.h"
#include "ConvertibleBondPDE.h"
#include "RootFinding.h"
#include <fstream>
#include <string>
#include <windows.h>

#define ERRORLEVEL     1e-6
#define IS_ALMOST_ZERO(x) ((x) <= ERRORLEVEL && (x) >= -ERRORLEVEL)

#define MACHINE_EPSILON 1e-9//2.220446049250313e-16
#define ALMOST_EQUAL(x, y) ((x) <= y * (1 + MACHINE_EPSILON) && (x) >= y * (1 - MACHINE_EPSILON))

struct GetImpSpreadOutput {
	double impliedSpread;
	bool indicator;
};

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
);

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
);

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
);