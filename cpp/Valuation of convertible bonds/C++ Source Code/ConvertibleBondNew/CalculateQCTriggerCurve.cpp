/*
** Author: Wu Chen
** 2017/10/25
*/

#include "CalculateQTriggerCurve.h"

QTriggerCurveReturnType calculateQCTriggerCurve(
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
	std::ofstream out("1.txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	double d_t = 1.0 / 252.0;
	double d_t_sqrt = sqrt(d_t);
	int M;

	if (T > noCallTime[0]) {
		double time_gap = T - noCallTime[0];
		M = ceil((T - time_gap) / d_t);
	}
	else {
		M = round(T / d_t);
	}

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

	// initialization
	std::vector<double> S(M);
	double windowSum, temp, windowSizeBeforeS;
	bool windowSumStop = false;
	int windowSize = m;
	double S_history_sum = 0.0;
	int path_count = 0;

	// initialize the TriggerCurve
	QTriggerCurve zero;
	zero.count = 0;
	zero.stockValueSum = 0.0;
	zero.time = 0.0;
	std::vector<QTriggerCurve> QCTriggerCurve(M, zero);

	// Case 1: T > noCallTime[0]
	// In this case, S_history = null
	if (T > noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 1" << std::endl;
		// S_history = null
		assert(S_history.empty());
		double time_gap = T - noCallTime[0];

		while (path_count < N) { // path simulation
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + j * d_t;
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + j * d_t;
					}
				}
			}
		}
	}

	// Case 2: 0 < S_history.size() < m && T <= noCallTime[0]
	// In this case, S_history is truncated. The length of S_history < m
	if (S_history.size() < m && T <= noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 2" << std::endl;

		while (path_count < N) { // path simulation
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + (j + 1) * d_t;
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + (j + 1) * d_t;
					}
				}
			}
		}
	}

	// Case 3: S_history.size() = m && T <= noCallTime[0]
	// In this case, S_history is not truncated. The length of S_history = m
	if (S_history.size() == m && T <= noCallTime[0]) {
		std::cout << "        Monte Carlo is running: QCTrigger - Case 3" << std::endl;

		while (path_count < N) { // path simulation
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + (j + 1) * d_t;
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
						QCTriggerCurve[j].count++;
						QCTriggerCurve[j].stockValueSum += S[j];
						QCTriggerCurve[j].time = noCallTime[1] - noCallTime[0] + (j + 1) * d_t;
					}
				}
			}
		}
	}
	std::cout << "total path: " << path_count << std::endl;
	QTriggerCurveReturnType result;
	std::vector<double> zeros(M, 0);
	result.time = zeros;
	result.QTrigger = zeros;
	for (int i = 0; i < QCTriggerCurve.size(); i++)
	{
		result.time[i] = noCallTime[1] - QCTriggerCurve[i].time;
		result.QTrigger[i] = QCTriggerCurve[i].stockValueSum / QCTriggerCurve[i].count;
	}
	return result;
}

void printTriggerCurve(std::vector<QTriggerCurve> QCTriggerCurve, double T)
{
	for (int i = 0; i < QCTriggerCurve.size(); i++)
	{
		std::cout << "i = " << i << std::endl;
		std::cout << "Time to Maturity:" << T - QCTriggerCurve[i].time << std::endl;
		std::cout << "QCTrigger = " << QCTriggerCurve[i].stockValueSum / QCTriggerCurve[i].count << std::endl;
		std::cout << "Path_number = " << QCTriggerCurve[i].count << std::endl;
	}
}

void printTriggerCurve2(QTriggerCurveReturnType QCTriggerCurve)
{
	for (int i = 0; i < QCTriggerCurve.time.size(); i++)
	{
		std::cout << "Time to Maturity:" << QCTriggerCurve.time[i] << std::endl;
		std::cout << "QCTrigger = " << QCTriggerCurve.QTrigger[i] << std::endl;
	}
}

int test_QCNew() {

	double S0 = 7.68, conversionPrice = 5.73, T = 5.6, sigma = 0.4;
	double faceValue = 100;
	int nCall = 15, mCall = 30, nPut = 30, mPut = 30;
	double q = 0, p = 0, R = 1, eta = 0;
	double PTrigger = 0.7 * conversionPrice;
	double CTrigger = 1.3 * conversionPrice;
	int theta = 1;
	double annuFactor = 365.25;
	std::string windcode = "132001.SH";
	std::string date = "2015-01-01";

	double r = 0.035;
	double value;
	double noCallTime1[2] = { 5.5, 6.0 };
	double noPutTime1[2] = { 2, 6.0 };
	double noConvertTime1[2] = { 5.5, 6.0 };
	double couponTimeRate1[6] = { 0.018, 0.015, 0.015, 0.01, 0.007, 0.005 };
	std::vector<double> noCallTime(noCallTime1, noCallTime1 + 2);
	std::vector<double> noPutTime(noPutTime1, noPutTime1 + 2);
	std::vector<double> noConvertTime(noConvertTime1, noConvertTime1 + 2);
	std::vector<double> couponTimeRate(couponTimeRate1, couponTimeRate1 + 6);

	double S_history1[30] = { 6.25, 6.13, 6.23, 6.03, 6.19, 6.32, 6.39, 6.33, 6.20, 6.04, 6.13, 6.12, 6.11, 6.25, 6.19, 6.09, 6.01, 6.13, 6.13, 6.24, 6.33, 6.80, 7.48, 7.70, 8.30, 7.75, 7.77, 7.81, 7.61, 7.68 };
	double CTrigger_history1[30] = { CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger };
	double PTrigger_history1[30] = { PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger };
	std::vector<double> S_history = {};
	std::vector<double> CTrigger_history = { CTrigger_history1 , CTrigger_history1 + 30 };
	std::vector<double> PTrigger_history = { PTrigger_history1 , PTrigger_history1 + 30 };

	double bonusRate = 0.042;

	double mktPrice = 133.722;

	clock_t start;
	double duration;
	start = clock();

	double QCTrigger, QPTrigger;
	int N = 500;

	QTriggerCurveReturnType result = calculateQCTriggerCurve(S0, CTrigger, T, r, sigma, nCall, mCall, N, S_history, CTrigger_history, noCallTime);
	//std::cout << result << std::endl;
	printTriggerCurve2(result);

	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "The running time is " << duration << std::endl;

	system("pause");
	return 0.0;
}

int main_Curve()
{
	test_QCNew();
	return 0;
}




