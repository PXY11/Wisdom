/*
** Author: Wu Chen
** 2017/07/05
*/

#include "CalculateQPTrigger.h"

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
	std::vector<double> noPutTime         /* No allow put time period          */)
{
	double d_t = 1.0 / 252.0;
	double d_t_sqrt = sqrt(d_t);
	int M = round(T / d_t);

	// Ziggurat algorithm of Marsaglia and Tsang
	float fn[128];
	uint32_t kn[128];
	uint32_t seed;
	float value;
	float wn[128];
	r4_nor_setup(kn, fn, wn);
	seed = rand() % 10000;

	// size shouble be equal
	assert(S_history.size() == PTrigger_history.size());

	// Calculate startWindowSum
	int startWindowSum = 0;
	for (int i = 0; i < S_history.size(); i++) {
		if (S_history[i] < PTrigger_history[i]) {
			startWindowSum++;
		}
	}

	if (startWindowSum >= n) {
		return -1;
	}

	// initialization
	std::vector<double> S(M);
	double windowSum, temp, windowSizeBeforeS;
	bool windowSumStop = false;
	bool triggerSumStop = false;
	int windowSize = m;
	double S_history_sum = 0.0;
	int path_count = 0;

	// Case 1: T > noPutTime[0]
	// In this case, S_history = null
	if (T > noPutTime[0]) {
		std::cout << "        Monte Carlo is running: QPTrigger - Case 1" << std::endl;
		// S_history = null
		double time_gap = T - noPutTime[0];
		M = ceil((T - time_gap) / d_t);

		for (int i = 1; i <= N; i++) { // path simulation
			// Initial value for each path
			windowSumStop = false;
			triggerSumStop = false;
			windowSum = 0.0;
			temp = r4_nor(seed, kn, fn, wn) * sqrt(time_gap);
			S[0] = S0 * exp((r - 0.5 * sigma * sigma) * time_gap + sigma * temp);
			// PTrigger
			if (S[0] < PTrigger) {
				windowSum++;
			}

			// j = 1 to m - 1
			for (int j = 1; j < m; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j * d_t + time_gap) + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}

			// j = m to M - 1
			for (int j = m; j < M; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j * d_t + time_gap) + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS < PTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}
		}
	}

	// Case 2: S_history.size() < m
	// In this case, S_history is truncated. The length of S_history < m
	if (S_history.size() < m && T <= noPutTime[0]) {
		std::cout << "        Monte Carlo is running: QPTrigger - Case 2" << std::endl;

		for (int i = 1; i <= N; i++) { // path simulation
									  // Initial value for each path
			temp = 0;
			windowSumStop = false;
			triggerSumStop = false;
			windowSum = startWindowSum;

			// j = 0 to m - 1
			for (int j = 0; j < m; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					if (j >= m - S_history.size()) {
						windowSizeBeforeS = S_history[j - m + S_history.size()];
						if (windowSizeBeforeS < PTrigger && windowSum > 0) {
							windowSum--;
						}
					}
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}
			// j = m to M - 1
			for (int j = m; j < M; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS < PTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}
		}
	}

	// Case 3: S_history.size() = m
	// In this case, S_history is not truncated. The length of S_history = m
	if (S_history.size() == m && T <= noPutTime[0]) {
		std::cout << "        Monte Carlo is running: QPTrigger - Case 3" << std::endl;

		for (int i = 1; i <= N; i++) { // path simulation
									  // Initial value for each path
			temp = 0;
			windowSumStop = false;
			triggerSumStop = false;
			windowSum = startWindowSum;

			// j = 0 to m
			for (int j = 0; j < m; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S_history[j];
					if (windowSizeBeforeS < PTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}
			// j = m to M - 1
			for (int j = m; j < M; j++) {
				if (triggerSumStop && windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// PTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS < PTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] < PTrigger) {
						windowSum++;
					}
					if (windowSum >= n) {
						windowSumStop = true;
						S_history_sum += S[j];
						path_count++;
					}
				}
			}
		}
	}
	std::cout << "total path: " << path_count << std::endl;

	if (path_count < N * 1e-5)
		return S_history_sum / path_count * 0.7;
	else
		return S_history_sum / path_count;
}