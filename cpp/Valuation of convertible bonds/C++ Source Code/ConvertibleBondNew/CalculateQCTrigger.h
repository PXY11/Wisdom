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

#define CONVERTIBLE_BOND_WAS_CALLED -1;

double calculateQCTrigger(
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
	std::vector<double> noCallTime        /* No allow call time period         */);