/*
** Author: Wu Chen
** 2017/07/05
*/

#include "CalculateQCTrigger.h"

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
	std::vector<double> noCallTime        /* No allow call time period         */)
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
	assert(S_history.size() == CTrigger_history.size());

	// Calculate startWindowSum
	int startWindowSum = 0;
	for (int i = 0; i < S_history.size(); i++) {
		if (S_history[i] >= CTrigger_history[i]) {
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
	int windowSize = m;
	double S_history_sum = 0.0;
	int path_count = 0;

	// Case 1: T > noCallTime[0]
	// In this case, S_history = null
	if (T > noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 1" << std::endl;
		// S_history = null
		assert(S_history.empty());
		double time_gap = T - noCallTime[0];
		M = ceil((T - time_gap) / d_t);

		for (int i = 1; i <= N; i++) { // path simulation
			// Initial value for each path
			windowSumStop = false;
			windowSum = 0.0;
			temp = r4_nor(seed, kn, fn, wn) * sqrt(time_gap);

			// j = 0
			S[0] = S0 * exp((r - 0.5 * sigma * sigma) * time_gap + sigma * temp);
			// CTrigger
			if (S[0] >= CTrigger) {
				windowSum++;
			}


			// j = 1 to m - 1
			for (int j = 1; j < m; j++) {
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j * d_t + time_gap) + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					if (S[j] >= CTrigger) {
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
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j * d_t + time_gap) + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS >= CTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] >= CTrigger) {
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

	// Case 2: 0 < S_history.size() < m && T <= noCallTime[0]
	// In this case, S_history is truncated. The length of S_history < m
	if (S_history.size() < m && T <= noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 2" << std::endl;

		for (int i = 1; i <= N; i++) { // path simulation
			// Initial value for each path
			temp = 0;
			windowSumStop = false;
			windowSum = startWindowSum;

			// j = 0 to m - 1
			for (int j = 0; j < m; j++) {
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					if (j >= m - S_history.size()) {
						windowSizeBeforeS = S_history[j - m + S_history.size()];
						if (windowSizeBeforeS >= CTrigger && windowSum > 0) {
							windowSum--;
						}
					}
					if (S[j] >= CTrigger) {
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
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS >= CTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] >= CTrigger) {
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

	// Case 3: S_history.size() = m && T <= noCallTime[0]
	// In this case, S_history is not truncated. The length of S_history = m
	if (S_history.size() == m && T <= noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 3" << std::endl;

		for (int i = 1; i <= N; i++) { // path simulation
									  // Initial value for each path
			temp = 0;
			windowSumStop = false;
			windowSum = startWindowSum;

			// j = 0 to m
			for (int j = 0; j < m; j++) {
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S_history[j];
					if (windowSizeBeforeS >= CTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] >= CTrigger) {
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
				if (windowSumStop) {
					break;
				}
				temp += r4_nor(seed, kn, fn, wn) * sqrt(d_t);
				S[j] = S0 * exp((r - 0.5 * sigma * sigma) * (j + 1) * d_t + sigma * temp);
				// CTrigger
				if (!windowSumStop) {
					windowSizeBeforeS = S[j - windowSize];
					if (windowSizeBeforeS >= CTrigger && windowSum > 0) {
						windowSum--;
					}
					if (S[j] >= CTrigger) {
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
		return S_history_sum / path_count * 2.0;
	else
		return S_history_sum / path_count;
}


