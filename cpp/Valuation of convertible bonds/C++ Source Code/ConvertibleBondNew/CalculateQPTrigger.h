/*
** Author: Wu Chen
** 2017/07/05
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include "random.h"
#include <assert.h>

#define CONVERTIBLE_BOND_MAY_BE_PUTTED -2;

double calculateQPTrigger(
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
	std::vector<double> noPutTime         /* No allow put time period          */);