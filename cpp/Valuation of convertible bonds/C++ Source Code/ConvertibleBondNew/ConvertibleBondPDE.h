/*
** Author: Wu Chen
** 2017/07/03
*/
#pragma once
#include <iostream>
#include <vector>
#include <ctime>
#include "pde_util.h"
#include "PricingRoutine.h"

struct PDEOutput {
	double value;
	double delta;
	double gamma;
	double theta;
};

struct FullPDEOutput {
	double value;
	double delta;
	double gamma;
	double theta;
	double vega;
	double rho;
};

FullPDEOutput ConvertibleBondPDEFullResult(
	int Nt,                              /* Time direction discretization number     */
	int Ns,                              /* Stock direction discretization number    */
	double faceValue,                    /* Convertible bond face value              */
	double conversionPrice,              /* Conversion price                         */
	double bonusRate,                    /* Bonus Rate for holding to expiry         */
	double S0,                           /* Current stock price                      */
	double sigma,                        /* Implied Volatility                       */
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
	double oas,                          /* Option-adjusted spread                   */
	double annuFactor,                   /* Annualized Factor                        */
	double callRedeemPrice,              /* Call option redemption price             */
	double putResellPrice,               /* Put option reselling price               */
	std::vector<double> noCallTime,      /* No allow call time period                */
	std::vector<double> noPutTime,       /* No allow put time period                 */
	std::vector<double> noConvertTime,   /* No allow conversion time period          */
	std::vector<double> couponTimeRate,  /* Coupon rate for different times          */
	std::string windcode,                /* windcode for current bond                */
	std::string date                     /* calculation date for current bond        */
);

PDEOutput ConvertibleBondPDE(
	int Nt,                              /* Time direction discretization number     */
	int Ns,                              /* Stock direction discretization number    */
	double faceValue,                    /* Convertible bond face value              */
	double conversionPrice,              /* Conversion price                         */
	double bonusRate,                    /* Bonus Rate for holding to expiry         */
	double S0,                           /* Current stock price                      */
	double sigma,                        /* Implied Volatility                       */
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
	double oas,                          /* Option-adjusted spread                   */
	double annuFactor,                   /* Annualized Factor                        */
	double callRedeemPrice,              /* Call option redemption price             */
	double putResellPrice,               /* Put option reselling price               */
	std::vector<double> noCallTime,      /* No allow call time period                */
	std::vector<double> noPutTime,       /* No allow put time period                 */
	std::vector<double> noConvertTime,   /* No allow conversion time period          */
	std::vector<double> couponTimeRate   /* Coupon rate for different times          */
);

double ConvertibleBondPDEPrice(
	int Nt,                              /* Time direction discretization number     */
	int Ns,                              /* Stock direction discretization number    */
	double faceValue,                    /* Convertible bond face value              */
	double conversionPrice,              /* Conversion price                         */
	double bonusRate,                    /* Bonus Rate for holding to expiry         */
	double S0,                           /* Current stock price                      */
	double sigma,                        /* Implied Volatility                       */
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
	double oas,                          /* Option-adjusted spread                   */
	double annuFactor,                   /* Annualized Factor                        */
	double callRedeemPrice,              /* Call option redemption price             */
	double putResellPrice,               /* Put option reselling price               */
	std::vector<double> noCallTime,      /* No allow call time period                */
	std::vector<double> noPutTime,       /* No allow put time period                 */
	std::vector<double> noConvertTime,   /* No allow conversion time period          */
	std::vector<double> couponTimeRate   /* Coupon rate for different times          */
);