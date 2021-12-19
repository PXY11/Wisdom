/*
** Author: Wu Chen
** 2017/07/03
*/

#include "ConvertibleBondPDE.h"

double ConvertibleBondPDEPriceNew(
	int Nt,                              /* Time direction discretization number     */
	int Ns,                              /* Stock direction discretization number    */
	double faceValue,                    /* Convertible bond face value              */
	double conversionPrice,              /* Conversion price                         */
	double Smax,                         /* Upper bound of stock price               */
	double S0,                           /* Current stock price                      */
	double sigma,                        /* Implied Volatility                       */
	double T,                            /* Time to maturity                         */
	double finalCouponRate,              /* Coupon rate for last payment             */
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
	std::vector<double> noCallTime,      /* No allow call time period                */
	std::vector<double> noPutTime,       /* No allow put time period                 */
	std::vector<double> noConvertTime,   /* No allow conversion time period          */
	std::vector<double> couponTimeRate   /* Coupon rate for different times          */
)
{
	// constant interest rate cases
	r = r + oas;
	double D = 0.0; // fixed dividend
	double dtime = 0.33; // dividend at time t = X.33
	double ds = Smax / Ns;
	double dt = T / Nt;
	double rhoPenaltyCall = 1000000 / (dt * dt);
	double rhoPenaltyPut = 1000000 / (dt * dt);

	// Boundary Condition
	double Bc_T = faceValue + couponTimeRate[0] * faceValue; // Bc in the final time, using dirty price
	double Bp_T;

	if (QPTrigger > 0) {
		Bp_T = faceValue + couponTimeRate[0] * faceValue; // Bp in the final time, using dirty price
	}
	else {
		Bp_T = 0.0;
	}

	double k_T = faceValue / (conversionPrice - D * ceil(T - dtime)); // dividend influence
	std::vector<double> I(Ns + 1);
	std::vector<double> S(Ns + 1);
	std::vector<double> u(Ns + 1);
	std::vector<double> B(Ns + 1);
	std::vector<double> dsVec(Ns + 1);
	for (int i = 0; i < Ns + 1; i++) {
		I[i] = i + 1;  // stock price, [1, 2, ..., Ns + 1], u(x, y) = u(T - xdt, (y - 1)dh)
		S[i] = ds * i; // stock price, [0, h, 2h, ..., Ns*h]
		double temp = _MIN(_MAX(k_T * ds * i, faceValue + couponTimeRate[0] * faceValue + finalCouponRate * faceValue), _MAX(Bc_T, k_T * QCTrigger));
		u[i] = temp;
		B[i] = faceValue + couponTimeRate[0] * faceValue + finalCouponRate * faceValue;
	}

	// Solve PDE by finite difference method, combined with penalty method
	int count = 0;
	int count1 = 0;
	int count2 = 0;
	int count3 = 0;
	int count4 = 0;
	int count5 = 0;
	int flag = 0; // fi not converge, break the whole loop

				  // parameters declaration
	std::vector<double> alpha(Ns + 1);
	std::vector<double> beta(Ns + 1);
	std::vector<double> gammaMu(Ns + 1);
	std::vector<double> gammaMb(Ns + 1);
	std::vector<double> rhs(Ns + 1);
	std::vector<double> u_n;
	std::vector<double> Muu;
	std::vector<double> P1_old;
	std::vector<double> P2_old;
	std::vector<double> u_old;
	std::vector<double> P1;
	std::vector<double> P2;
	std::vector<double> n_upper(Ns + 1);
	std::vector<double> n_lower(Ns + 1);

	double Bc_n, Bp_n, k_n, t, AccI;

	/*
	***************************************************
	*   Constant Interest Rate Cases
	***************************************************
	*/

	// coefficients of the discretization
	// Mu = spdiags([alpha, gammaMu, beta],(-1:1),Ns+1,Ns+1)
	// Mb = spdiags([alpha, gammaMb, beta],(-1:1),Ns+1,Ns+1)
	//               lower, middle , upper
	for (int j = 0; j < Ns + 1; j++) {
		if (j < Ns - 1) {
			alpha[j] = dt * (pow(sigma, 2.0) * pow((j + 1) * ds, 2.0) / 2.0 / pow(ds, 2.0)
				- (r + p * eta - q) * (j + 1) * ds / 2.0 / ds);
		}
		else {
			alpha[j] = 0.0;
		}
		if (j < 2) {
			beta[j] = 0.0;
		}
		else {
			beta[j] = dt * (pow(sigma, 2.0) * pow((j - 1) * ds, 2.0) / 2.0 / pow(ds, 2.0)
				+ (r + p * eta - q) * (j - 1) * ds / 2.0 / ds);
		}
		if (j == 0) {
			gammaMu[j] = -(r + p) * dt;
			gammaMb[j] = -(r + p * (1 - R)) * dt;

		}
		else if (j == Ns) {
			gammaMu[j] = 0.0;
			gammaMb[j] = 0.0;
		}
		else {
			double beta = dt * (pow(sigma, 2.0) * pow(j * ds, 2.0) / 2.0 / pow(ds, 2.0)
				+ (r + p * eta - q) * j * ds / 2.0 / ds);
			gammaMu[j] = -(alpha[j - 1] + beta + (r + p) * dt);
			gammaMb[j] = -(alpha[j - 1] + beta + (r + p * (1 - R)) * dt);
		}
	}

	for (int i = 0; i < Nt; i++) {

		if (flag == 1) {
			break;
		}
		t = dt * (i + 1);

		// dividend influence to stock price
		for (int j = 0; j < S.size(); j++) {
			S[j] += D * ceil(t - dtime);
		}
		// coupon rate payment
		if ((ceil(t + dt) - ceil(t)) == 1.0) {
			count5++;
			for (int j = 0; j < u.size(); j++) {
				u[j] += couponTimeRate[ceil(t + dt + 1e-6) - 1] * faceValue;
				B[j] += couponTimeRate[ceil(t + dt + 1e-6) - 1] * faceValue;
			}
		}

		AccI = (ceil(t) - t) * couponTimeRate[ceil(t) - 1] * faceValue;


		if (noCallTime[0] <= t && t <= noCallTime[1] || QCTrigger == 10000000) {
			Bc_n = 10000000; // when n value of Bc, take a large number means no call
		}
		else {
			Bc_n = faceValue + AccI; // dirty price
		}
		if (noPutTime[0] <= t && t <= noPutTime[1] || QPTrigger == 0.0) {
			Bp_n = 0.0; // small number means no put
		}
		else {
			Bp_n = faceValue + AccI; // dirty price
		}
		if (noConvertTime[0] <= t && t <= noConvertTime[1]) {
			k_n = 0.0; // small number means no convert
		}
		else {
			k_n = faceValue / (conversionPrice - D * (ceil(T - dtime) - ceil(t - dtime)));
		}

		// now value of u is n-1 step, will convert to n step following
		u_n = u;
		// linear system: Ax = rhs, where A = (I - theta * Mb), b = (I + (1 - theta) * Mb) and  rhs = b * B
		my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMb, (1 - theta) * beta, B, rhs);
		Thomas_Algorithm(-theta * alpha, 1 + -theta * gammaMb, -theta * beta, rhs, B);
		for (int j = 0; j < B.size(); j++) {
			if (B[j] > Bc_n) { // explicit B <= Bc
				B[j] = Bc_n;
			}
		}

		// [Muu(1:Ns,:);0]
		Muu = p * dt * vecMAX(k_n * (1 - eta) * S, R * B);
		Muu[Ns] = 0.0;
		my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMu, (1 - theta) * beta, u_n, rhs);
		rhs = rhs + Muu;

		// solve the linear system
		Thomas_Algorithm(-theta * alpha, 1 + -theta *gammaMu, -theta * beta, rhs, u);

		// check upper bound and lower bound
		// new version
		for (int j = 0; j < S.size(); j++) {
			if (S[j] < QPTrigger) {
				u[j] = Bp_n;
				n_lower[j] = Bp_n;
				n_upper[j] = Bp_n;
			}
			else if (S[j] >= QCTrigger) {
				u[j] = k_n * QCTrigger;
				n_lower[j] = k_n * QCTrigger;
				n_upper[j] = k_n * QCTrigger;
			}
			else {
				u[j] = _MAX(u[j], k_n * S[j]);
				n_lower[j] = k_n * S[j];
				n_upper[j] = 10000000;
			}
		}
		// MAX(k_n * QCTrigger, Bc_n) for non callable cases
		// following is the old version: wrong!
		// u = numvecMIN(_MAX(k_n * QCTrigger, Bc_n), numvecMAX(Bp_n, vecMAX(u, k_n * S)));
		P1 = rhoPenaltyPut * compare(u, n_lower);
		P2 = -rhoPenaltyCall * compare(n_upper, u);

		// loop for penalty method
		while (true) {
			count++;
			u_old = u;
			P1_old = P1;
			P2_old = P2;
			// solve the linear system
			// speye(Ns+1) + (1-theta)*Mu)*u_n
			my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMu, (1 - theta) * beta, u_n, rhs);
			// + [Muu(1:Ns,:);0]
			rhs = rhs + Muu;
			// + P1_old*max(Bp_n,k_n*S)-P2_old*max(Bc_n,k_n*qctrigger)
			rhs = rhs + element_prod(P1_old, n_lower) - element_prod(P2_old, n_upper);
			Thomas_Algorithm(-theta * alpha, 1 + (-theta * gammaMu + P1_old - P2_old), -theta * beta, rhs, u);
			// new panelty
			P1 = rhoPenaltyPut * compare(u, n_lower);
			P2 = -rhoPenaltyCall * compare(n_upper, u);
			if (P1 == P1_old && P2 == P2_old) {
				count1 = count1 + 1;
				break;
			}
			else if (maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) < (1 / 0.0001)) {
				count2++;
				break;
			}
			else if (count > Nt * 40) {
				count4++;
				std::cout << "Probably does not converge.\n" << std::endl;
				flag = 1;
				break;
			}
			else {
				count3++;
			}
		}
	}
	double Uindex = Ns * (S0 - S[0]) / (S[Ns] - S[0]) + 1;
	// convertible bond value
	double value = u[ceil(Uindex) - 1] * (Uindex - floor(Uindex + 0.001)) + u[floor(Uindex) - 1] * (ceil(Uindex + 0.001) - Uindex);
	return value;
}


int TestPDE()
{
	int Nt = ceil(252 * 3.4904); // Time direction
	int Ns = 500; // Price direction
	double T = 3.4904; // Time to maturaty
	double faceValue = 100.0; // Face value
	double conversionPrice = 9.26; // Conversion price
	double S0 = 7.95;// Current stock price
	double Smax = 4 * S0; // Upbound of stock price
	double sigma = 1.3176; // Stock volatility
	double finalCouponRate = 0.04; // Coupon rate for last payment
	double q = 0; // Dividend
	double p = 0; // Default probobility
	double R = 1; // recovery rate
	double eta = 0; // when default the stock price = (1 - eta)*S
	double QCTrigger = 15.360488000000002; // callable price
	double QPTrigger = 10; // callable price
	double theta = 1; // implicitness parameter, 0.5 for CN method, 1 for fully implicit
	double r = 0.0345; // risk free rate
	double annuFactor = 252.0; // annualized factor
	double oas = 0.1; // option-adjusted spread

					  // for later usage
					  // std::vector<double> discFactor;      /* Discount Factor curve data              */
					  // std::vector<double> discTerm;        /* Discount Factor Term Structure          */

	double noCallTime1[2] = { 0, 5.0 };
	double noPutTime1[2] = { 0, 5.0 };
	double noConvertTime1[2] = { 3.0, 5.0 };
	double couponTimeRate1[6] = { 0.02, 0.015, 0.01, 0.008, 0.006 };
	std::vector<double> noCallTime(noCallTime1, noCallTime1 + 2);               /* No allow call time period               */
	std::vector<double> noPutTime(noPutTime1, noPutTime1 + 2);                  /* No allow put time period                */
	std::vector<double> noConvertTime(noConvertTime1, noConvertTime1 + 2);      /* No allow conversion time period         */
	std::vector<double> couponTimeRate(couponTimeRate1, couponTimeRate1 + 12);  /* Coupon rate for different times         */

	clock_t start;
	double duration;
	start = clock();

	FullPDEOutput result;

	result = ConvertibleBondPDEFullResult(Nt, Ns, faceValue, conversionPrice, Smax, S0, sigma, T, finalCouponRate, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, noCallTime, noPutTime, noConvertTime, couponTimeRate);

	std::cout << "bond value:" << result.value << std::endl;
	std::cout << "delta:" << result.delta << std::endl;
	std::cout << "gamma:" << result.gamma << std::endl;
	std::cout << "theta:" << result.theta << std::endl;
	std::cout << "rho:" << result.rho << std::endl;
	std::cout << "vega:" << result.vega << std::endl;

	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "The running time is " << duration << std::endl;

	return 0;
}