/*
** Author: Wu Chen
** 2017/07/03
*/

#include "ConvertibleBondPDE.h"
#include "CalculateQCTrigger.h"
#include "CalculateQPTrigger.h"

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
)
{
	std::wstring FilePath = L".\\logs_GetPDEFullResult\\";
	CreateDirectoryW(FilePath.c_str(), NULL);
	std::ofstream out(".\\logs_GetPDEFullResult\\" + windcode + " " + date + ".txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	FullPDEOutput result;
	PDEOutput output, output2, output_rho, output_rho2, output_vega, output_vega2;
	double rho, vega;

	// Richardson extrapolation
	// value, delta, gamma, theta
	output = ConvertibleBondPDE(Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);
	// value, delta, gamma, theta
	output2 = ConvertibleBondPDE(Nt * 4, Ns * 2, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);

	// rho
	output_rho = ConvertibleBondPDE(Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r + 0.0001, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);
	output_rho2 = ConvertibleBondPDE(Nt * 4, Ns * 2, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r + 0.0001, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);
	rho = (4 * (output_rho2.value - output2.value) / 0.0001 / 100 - (output_rho.value - output.value) / 0.0001 / 100) / 3.0;
	// vega
	output_vega = ConvertibleBondPDE(Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma + 0.0001, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);
	output_vega2 = ConvertibleBondPDE(Nt * 4, Ns * 2, faceValue, conversionPrice, bonusRate, S0, sigma + 0.0001, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);
	vega = (4 * (output_vega2.value - output2.value) / 0.0001 / 100 - (output_vega.value - output.value) / 0.0001 / 100) / 3.0;

	// results
	result.value = (4 * output2.value - output.value) / 3.0;
	result.delta = (4 * output2.delta - output.delta) / 3.0;
	result.gamma = (4 * output2.gamma - output.gamma) / 3.0;
	result.theta = (4 * output2.theta - output.theta) / 3.0;
	result.rho = rho;
	result.vega = vega;

	std::cout << "bond value:" << result.value << std::endl;
	std::cout << "delta:" << result.delta << std::endl;
	std::cout << "gamma:" << result.gamma << std::endl;
	std::cout << "theta:" << result.theta << std::endl;
	std::cout << "rho:" << result.rho << std::endl;
	std::cout << "vega:" << result.vega << std::endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
	return result;
}

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
)
{
	double Smax = _MAX(S0 * 5, QCTrigger * 2);
	double ds = QCTrigger / Ns;
	Ns = round(Smax / ds);
	double dt = T / Nt;
	double Large = 1 / _MIN(ds * ds, pow(dt, 1.5));
	double rhoPenaltyCall = Large;
	double rhoPenaltyPut = Large;
	double tol = 1 / Large;

	// constant interest rate cases
	r = r + oas;
	double D = 0.0; // fixed dividend
	double dtime = 0.33; // dividend at time t = X.33

						 // Boundary Condition
	double Bc_T = faceValue + couponTimeRate[0] * faceValue;// +bonusRate * faceValue; // Bc in the final time, using dirty price

															// if the put provision exists
	double Bp_T;
	if (QPTrigger > 0) {
		Bp_T = faceValue + couponTimeRate[0] * faceValue;// +bonusRate * faceValue; // Bp in the final time, using dirty price
	}
	else {
		Bp_T = 0.0;
	}

	double k_T = faceValue / (conversionPrice - D * ceil(T - dtime)); // dividend influence

	std::vector<double> S(Ns + 1);
	std::vector<double> u(Ns + 1);
	std::vector<double> B(Ns + 1);

	for (int i = 0; i < Ns + 1; i++) {
		S[i] = ds * i;
		//u[i] = _MAX(k_T * ds * i, faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue);
		if (S[i] < QCTrigger) {
			u[i] = _MAX(k_T * ds * i, faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue);
		}
		else {
			u[i] = k_T * ds * i;
		}
		B[i] = faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue;
	}

	// Solve PDE by finite difference method, combined with penalty method
	int count, count1, count2, count3, count4, count5;
	int flag = 0; // if not converge, break the whole loop

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

		count = 0;
		count1 = 0;
		count2 = 0;
		count3 = 0;
		count4 = 0;
		count5 = 0;

		if (flag == 1) {
			break;
		}
		t = dt * (i + 1);


		if (noCallTime[0] <= t && t < noCallTime[1]) {
			gammaMu[Ns] = -(q + p) * dt;
		}


		// dividend influence to stock price
		for (int j = 0; j < S.size(); j++) {
			S[j] += D * ceil(t - dtime);
		}

		// Accured Interest
		AccI = (ceil(t + dt - 1e-10) - t) * couponTimeRate[ceil(t + dt - 1e-10) - 1] * faceValue;

		// In no convert time period
		if (noConvertTime[0] <= t && t <= noConvertTime[1]) {
			k_n = 0.0; // small number means no convert
		}
		// In convert time period
		else {
			k_n = faceValue / (conversionPrice - D * (ceil(T - dtime) - ceil(t - dtime)));
		}
		// In no call time period
		if (noCallTime[0] <= t && t <= noCallTime[1] || QCTrigger == 10000000.0) {
			Bc_n = 10000000; // when n value of Bc, take a large number means no call
		}
		// In call time period
		else {
			if (callRedeemPrice == -1) {
				Bc_n = faceValue + AccI; // dirty price
			}
			else {
				Bc_n = callRedeemPrice;
			}
		}
		// In no put time period
		if (noPutTime[0] <= t && t <= noPutTime[1] || QPTrigger == 0.0) {
			Bp_n = 0.0; // small number means no put
		}
		// In put time period
		else {
			if (putResellPrice == -1) {
				Bp_n = faceValue + AccI; // dirty price
			}
			else {
				Bp_n = putResellPrice;
			}
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
		// 2017-07-31
		for (int j = 0; j < S.size(); j++) {
			if (S[j] < QPTrigger) {
				u[j] = _MAX(u[j], Bp_n);
				n_upper[j] = 100000;
				n_lower[j] = _MAX(Bp_n, k_n * S[j]);
			}
			else if (S[j] >= QCTrigger) {
				u[j] = _MAX(u[j], k_n * S[j]);
				u[j] = _MIN(u[j], _MAX(k_n * S[j], Bc_n));
				n_upper[j] = _MAX(k_n * S[j], Bc_n);
				n_lower[j] = k_n * S[j];
			}
			else {
				u[j] = _MAX(u[j], k_n * S[j]);
				n_lower[j] = k_n * S[j];
				n_upper[j] = 100000;
			}
		}

		// New Method
		if (!(noConvertTime[0] <= t && t < noConvertTime[1])) {
			std::vector<double> uu;
			uu = u;
			P1_old = rhoPenaltyPut * compare(u, n_lower);
			P2_old = rhoPenaltyCall * compare(n_upper, u);
			while (true) {
				count++;
				u_old = u;
				P1 = rhoPenaltyPut * compare(u + uu, 2 * n_lower);
				P2 = rhoPenaltyCall * compare(2 * n_upper, u + uu);
				// solve the linear system
				my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMu, (1 - theta) * beta, u_n, rhs);
				rhs = rhs + Muu;
				rhs = rhs - element_prod(P1, uu - 2 * n_lower) - element_prod(P2, uu - 2 * n_upper);
				Thomas_Algorithm(-theta * alpha, 1 + (-theta * gammaMu + P1 + P2), -theta * beta, rhs, u);
				// new panelty
				P1 = rhoPenaltyPut * compare_equal(u + uu, 2 * n_lower);
				P2 = rhoPenaltyCall * compare_equal(2 * n_upper, u + uu);
				if (P1 == P1_old && P2 == P2_old) {
					count1 = count1 + 1;
					break;
				}
				if (maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) < tol) {
					count2++;
					break;
				}
				else if (count > 1e4) {
					count4++;
					std::cout << "Probability does not converge.\n" << std::endl;
					std::cout << "t = " << t << " i = " << i << std::endl;
					std::cout << "err = " << maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) << std::endl;
					//printvector(u);
					flag = 1;
					break;
				}
				else {
					count3++;
				}
				P1_old = P1;
				P2_old = P2;
			}
		}

		B = vecMIN(B, u);
		if (floor(t + dt + 1e-10) - floor(t + 1e-10) == 1 && floor(t + dt + 1e-10) + 1 <= couponTimeRate.size()) {
			count5++;
			for (int i = 0; i < S.size(); i++) {
				if (S[i] < QCTrigger) {
					u[i] = u[i] + couponTimeRate[floor(t + dt + 1e-10)] * faceValue;
				}
			}
			B = couponTimeRate[floor(t + dt + 1e-10)] * faceValue + B;
		}
	}
	double Uindex = Ns * (S0 - S[0]) / (S[Ns] - S[0]) + 1;
	// convertible bond value
	double value = u[ceil(Uindex - 1e-10) - 1] * (Uindex - floor(Uindex + 1e-10)) + u[floor(Uindex + 1e-10) - 1] * (ceil(Uindex - 1e-10) - Uindex);
	// delta
	double delta = (u[ceil(Uindex - 1e-10) - 1] - u[floor(Uindex + 1e-10) - 1]) / ds / (100.0 / S0);
	double delta1 = (u[ceil(Uindex - 1e-10) - 2] - u[floor(Uindex + 1e-10) - 2]) / ds / (100.0 / S0);
	// gamma
	double gamma = (u[ceil(Uindex - 1e-10)] - 2 * u[ceil(Uindex - 1e-10) - 1] + u[floor(Uindex + 1e-10) - 1]) / pow(ds, 2) / (100.0 / S0) * 100;
	gamma = (delta - delta1) / ds / (100.0 / S0) * 100;
	// theta 
	double value2 = u_old[ceil(Uindex - 1e-10) - 1] * (Uindex - floor(Uindex + 1e-10)) + u_old[floor(Uindex + 1e-10) - 1] * (ceil(Uindex - 1e-10) - Uindex);
	double _theta = (value2 - value) / dt / annuFactor;

	PDEOutput result;
	result.value = value;
	result.delta = delta;
	result.gamma = gamma;
	result.theta = _theta;

	return result;
}

/*
************************************************************************************
** According to current chinese convertible bond contract, we assume that
** Call Start Time = Conversion Start Time  = t1 < Put Start Time = t2
** < Call End Time = Conversion End Time = Put End Time = Maturity Date = T
** Day Count Convention for calculation: ACT/ACT
************************************************************************************
*/
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
)
{

	double Smax = _MAX(S0 * 5, QCTrigger * 2);
	double ds = QCTrigger / Ns;
	Ns = round(Smax / ds);
	double dt = T / Nt;
	double Large = 1 / _MIN(ds * ds, pow(dt, 1.5));
	double rhoPenaltyCall = Large;
	double rhoPenaltyPut = Large;
	double tol = 1 / Large;

	// constant interest rate cases
	r = r + oas;
	double D = 0.0; // fixed dividend
	double dtime = 0.33; // dividend at time t = X.33

	// Boundary Condition
	double Bc_T = faceValue + couponTimeRate[0] * faceValue;// +bonusRate * faceValue; // Bc in the final time, using dirty price

	// if the put provision exists
	double Bp_T;
	if (QPTrigger > 0) {
		Bp_T = faceValue + couponTimeRate[0] * faceValue;// +bonusRate * faceValue; // Bp in the final time, using dirty price
	}
	else {
		Bp_T = 0.0;
	}

	double k_T = faceValue / (conversionPrice - D * ceil(T - dtime)); // dividend influence

	std::vector<double> S(Ns + 1);
	std::vector<double> u(Ns + 1);
	std::vector<double> B(Ns + 1);

	for (int i = 0; i < Ns + 1; i++) {
		S[i] = ds * i;
		//u[i] = _MAX(k_T * ds * i, faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue);
		if (S[i] < QCTrigger) {
			u[i] = _MAX(k_T * ds * i, faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue);
		}
		else {
			u[i] = k_T * ds * i;
		}
		B[i] = faceValue + couponTimeRate[0] * faceValue + bonusRate * faceValue;
	}

	// Solve PDE by finite difference method, combined with penalty method
	int count, count1, count2, count3, count4, count5;
	int flag = 0; // if not converge, break the whole loop

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

		count = 0;
		count1 = 0;
		count2 = 0;
		count3 = 0;
		count4 = 0;
		count5 = 0;

		if (flag == 1) {
			break;
		}
		t = dt * (i + 1);

		
		if (noCallTime[0] <= t && t < noCallTime[1]) {
			gammaMu[Ns] = -(q + p) * dt;
		}
		

		// dividend influence to stock price
		for (int j = 0; j < S.size(); j++) {
			S[j] += D * ceil(t - dtime);
		}

		// Accured Interest
		AccI = (ceil(t + dt - 1e-10) - t) * couponTimeRate[ceil(t + dt - 1e-10) - 1] * faceValue;

		// In no convert time period
		if (noConvertTime[0] <= t && t <= noConvertTime[1]) {
			k_n = 0.0; // small number means no convert
		}
		// In convert time period
		else {
			k_n = faceValue / (conversionPrice - D * (ceil(T - dtime) - ceil(t - dtime)));
		}
		// In no call time period
		if (noCallTime[0] <= t && t <= noCallTime[1] || QCTrigger == 10000000.0) {
			Bc_n = 10000000; // when n value of Bc, take a large number means no call
		}
		// In call time period
		else {
			if (callRedeemPrice == -1) {
				Bc_n = faceValue + AccI; // dirty price
			}
			else {
				Bc_n = callRedeemPrice;
			}
		}
		// In no put time period
		if (noPutTime[0] <= t && t <= noPutTime[1] || QPTrigger == 0.0) {
			Bp_n = 0.0; // small number means no put
		}
		// In put time period
		else {
			if (putResellPrice == -1) {
				Bp_n = faceValue + AccI; // dirty price
			}
			else {
				Bp_n = putResellPrice;
			}
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
		// 2017-07-31
		for (int j = 0; j < S.size(); j++) {
			if (S[j] < QPTrigger) {
				u[j] = _MAX(u[j], Bp_n);
				n_upper[j] = 100000;
				n_lower[j] = _MAX(Bp_n, k_n * S[j]);
			}
			else if (S[j] >= QCTrigger) {
				u[j] = _MAX(u[j], k_n * S[j]);
				u[j] = _MIN(u[j], _MAX(k_n * S[j], Bc_n));
				n_upper[j] =  _MAX(k_n * S[j], Bc_n);
				n_lower[j] = k_n * S[j];
			}
			else {
				u[j] = _MAX(u[j], k_n * S[j]);
				n_lower[j] = k_n * S[j];
				n_upper[j] = 100000;
			}
		}
		
		if (!(noConvertTime[0] <= t && t < noConvertTime[1])) {
			// New Method
			std::vector<double> uu;
			uu = u;
			P1_old = rhoPenaltyPut * compare(u, n_lower);
			P2_old = rhoPenaltyCall * compare(n_upper, u);
			while (true) {
				count++;
				u_old = u;
				P1 = rhoPenaltyPut * compare(u + uu, 2 * n_lower);
				P2 = rhoPenaltyCall * compare(2 * n_upper, u + uu);
				// solve the linear system
				my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMu, (1 - theta) * beta, u_n, rhs);
				rhs = rhs + Muu;
				rhs = rhs - element_prod(P1, uu - 2 * n_lower) - element_prod(P2, uu - 2 * n_upper);
				Thomas_Algorithm(-theta * alpha, 1 + (-theta * gammaMu + P1 + P2), -theta * beta, rhs, u);
				// new panelty
				P1 = rhoPenaltyPut * compare_equal(u + uu, 2 * n_lower);
				P2 = rhoPenaltyCall * compare_equal(2 * n_upper, u + uu);
				if (P1 == P1_old && P2 == P2_old) {
					count1 = count1 + 1;
					break;
				}
				if (maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) < tol) {
					count2++;
					break;
				}
				else if (count > 1e4) {
					count4++;
					std::cout << "Probability does not converge.\n" << std::endl;
					std::cout << "t = " << t << " i = " << i << std::endl;
					std::cout << "err = " << maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) << std::endl;
					//printvector(u);
					flag = 1;
					break;
				}
				else {
					count3++;
				}
				P1_old = P1;
				P2_old = P2;
			}
			/*
			// old method
			P1 = rhoPenaltyPut * compare(u, n_lower);
			P2 = rhoPenaltyCall * compare(n_upper, u);

			// loop for penalty method
			while (true) {
				count++;
				u_old = u;
				P1_old = P1;
				P2_old = P2;
				// solve the linear system
				my_dot_product((1 - theta) * alpha, 1 + (1 - theta) * gammaMu, (1 - theta) * beta, u_n, rhs);
				rhs = rhs + Muu;
				rhs = rhs + element_prod(P1_old, n_lower) + element_prod(P2_old, n_upper);
				Thomas_Algorithm(-theta * alpha, 1 + (-theta * gammaMu + P1_old + P2_old), -theta * beta, rhs, u);
				// new panelty
				P1 = rhoPenaltyPut * compare(u, n_lower);
				P2 = rhoPenaltyCall * compare(n_upper, u);
				if (P1 == P1_old && P2 == P2_old) {
					count1 = count1 + 1;
					break;
				}
				else if (maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) < 1e-4) {
					count2++;
					break;
				}
				else if (count > 1e4) {
					count4++;
					std::cout << "Probability does not converge.\n" << std::endl;
					std::cout << "t = " << t << " i = " << i << std::endl;
					std::cout << "err = " << maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) << std::endl;
					//printvector(u);
					flag = 1;
					break;
				}
				else {
					count3++;
				}
			}
			*/
		}
		
		//std::cout << "current n = " << i << " count = " << count << " err = " <<
		//	maxinvec(element_divide(abs_diff(u, u_old), numvecMAX(1.0, absvec(u)))) << std::endl;
	
		B = vecMIN(B, u);
		if (floor(t + dt + 1e-10) - floor(t + 1e-10) == 1 && floor(t + dt + 1e-10) + 1 <= couponTimeRate.size()) {
			count5++;
			for (int i = 0; i < S.size(); i++) {
				if (S[i] < QCTrigger) {
					u[i] = u[i] + couponTimeRate[floor(t + dt + 1e-10)] * faceValue;
				}
			}
			B = couponTimeRate[floor(t + dt + 1e-10)] * faceValue + B;
		}
	}
	double Uindex = Ns * (S0 - S[0]) / (S[Ns] - S[0]) + 1;
	// convertible bond value
	double value = u[ceil(Uindex - 1e-10) - 1] * (Uindex - floor(Uindex + 1e-10)) + u[floor(Uindex + 1e-10) - 1] * (ceil(Uindex - 1e-10) - Uindex);
	return value;
}

int test_PDE() {

	double S0 = 7.68, conversionPrice = 5.73, T = 1.474, sigma = 0.4;
	double faceValue = 100;
	int nCall = 15, mCall = 30, nPut = 30, mPut = 30;
	double q = 0.05, p = 0, R = 1, eta = 0;
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
	double CTrigger_history1[30] = { CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger, CTrigger};
	double PTrigger_history1[30] = { PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger, PTrigger};
	std::vector<double> S_history = { S_history1 , S_history1 + 30};
	std::vector<double> CTrigger_history = { CTrigger_history1 , CTrigger_history1 + 30};
	std::vector<double> PTrigger_history = { PTrigger_history1 , PTrigger_history1 + 30};

	double bonusRate = 0.042;

	double mktPrice = 133.722;

	clock_t start;
	double duration;
	start = clock();

	double QCTrigger, QPTrigger;
	int N = 50000;

	std::cout.precision(15);
	QCTrigger = 8.97664998628881;// calculateQCTrigger(S0, CTrigger, T, r, sigma, nCall, mCall, N, S_history, CTrigger_history, noCallTime);
	QPTrigger = 3.42962446334706;// calculateQPTrigger(S0, PTrigger, T, r, sigma, nPut, mPut, N, S_history, PTrigger_history, noPutTime);
	//std::cout << QCTrigger << std::endl;
	//std::cout << QPTrigger << std::endl;

	int Nt = 5000, Ns = 5000;
	double oas = 0;

	double callRedeemPrice = -1;
	double putResellPrice = -1;

	double result = ConvertibleBondPDEPrice(Nt, Ns, faceValue, conversionPrice, bonusRate, S0, sigma, T, q, p, R, eta, QCTrigger, QPTrigger,
		theta, r, oas, annuFactor, callRedeemPrice, putResellPrice, noCallTime, noPutTime, noConvertTime, couponTimeRate);

	std::cout.precision(18);
	std::cout << "bond value:" << result << std::endl;
	//std::cout << "delta:" << result.delta << std::endl;
	//std::cout << "gamma:" << result.gamma << std::endl;
	//std::cout << "theta:" << result.theta << std::endl;
	//std::cout << "rho:" << result.rho << std::endl;
	//std::cout << "vega:" << result.vega << std::endl;

	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "The running time is " << duration << std::endl;

	system("pause");
	return 0.0;
}

int main()
{
	test_PDE();
	return 0;
}