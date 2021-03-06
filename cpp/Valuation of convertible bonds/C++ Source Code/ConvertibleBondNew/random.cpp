#include "random.h"

float r4_uni(uint32_t &jsr)

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNI returns a uniformly distributed real value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed.
//
//    Output, float R4_UNI, a uniformly distributed random value in
//    the range [0,1].
//
{
	uint32_t jsr_input;
	float value;

	jsr_input = jsr;

	jsr = (jsr ^ (jsr << 13));
	jsr = (jsr ^ (jsr >> 17));
	jsr = (jsr ^ (jsr << 5));

	value = fmod(0.5
		+ (float)(jsr_input + jsr) / 65536.0 / 65536.0, 1.0);

	return value;
}

uint32_t shr3_seeded(uint32_t &jsr)

//****************************************************************************80
//
//  Purpose:
//
//    SHR3_SEEDED evaluates the SHR3 generator for integers.
//
//  Discussion:
//
//    Thanks to Dirk Eddelbuettel for pointing out that this code needed to
//    use the uint32_t data type in order to execute properly in 64 bit mode,
//    03 October 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed, which is updated 
//    on each call.
//
//    Output, uint32_t SHR3_SEEDED, the new value.
//
{
	uint32_t jsr_input;
	uint32_t value;

	jsr_input = jsr;

	jsr = (jsr ^ (jsr << 13));
	jsr = (jsr ^ (jsr >> 17));
	jsr = (jsr ^ (jsr << 5));

	value = jsr_input + jsr;

	return value;
}

float r4_nor(uint32_t &jsr, uint32_t kn[128], float fn[128], float wn[128])

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR returns a normally distributed single precision real value.
//
//  Discussion:
//
//    The value returned is generated from a distribution with mean 0 and 
//    variance 1.
//
//    The underlying algorithm is the ziggurat method.
//
//    Before the first call to this function, the user must call R4_NOR_SETUP
//    to determine the values of KN, FN and WN.
//
//    Thanks to Chad Wagner, 21 July 2014, for noticing a bug of the form
//      if ( x * x <= y * y );   <-- Stray semicolon!
//      {
//        break;
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed.
//
//    Input, uint32_t KN[128], data computed by R4_NOR_SETUP.
//
//    Input, float FN[128], WN[128], data computed by R4_NOR_SETUP.
//
//    Output, float R4_NOR, a normally distributed random value.
//
{
	int hz;
	uint32_t iz;
	const float r = 3.442620;
	float value;
	float x;
	float y;

	hz = (int)shr3_seeded(jsr);
	iz = (hz & 127);

	if (fabs(hz) < kn[iz])
	{
		value = (float)(hz)* wn[iz];
	}
	else
	{
		for (; ; )
		{
			if (iz == 0)
			{
				for (; ; )
				{
					x = -0.2904764 * log(r4_uni(jsr));
					y = -log(r4_uni(jsr));
					if (x * x <= y + y)
					{
						break;
					}
				}

				if (hz <= 0)
				{
					value = -r - x;
				}
				else
				{
					value = +r + x;
				}
				break;
			}

			x = (float)(hz)* wn[iz];

			if (fn[iz] + r4_uni(jsr) * (fn[iz - 1] - fn[iz])
				< exp(-0.5 * x * x))
			{
				value = x;
				break;
			}

			hz = (int)shr3_seeded(jsr);
			iz = (hz & 127);

			if (fabs(hz) < kn[iz])
			{
				value = (float)(hz)* wn[iz];
				break;
			}
		}
	}

	return value;
}
//****************************************************************************80

void r4_nor_setup(uint32_t kn[128], float fn[128], float wn[128])

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR_SETUP sets data needed by R4_NOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t KN[128], data needed by R4_NOR.
//
//    Output, float FN[128], WN[128], data needed by R4_NOR.
//
{
	double dn = 3.442619855899;
	int i;
	const double m1 = 2147483648.0;
	double q;
	double tn = 3.442619855899;
	const double vn = 9.91256303526217E-03;

	q = vn / exp(-0.5 * dn * dn);

	kn[0] = (uint32_t)((dn / q) * m1);
	kn[1] = 0;

	wn[0] = (float)(q / m1);
	wn[127] = (float)(dn / m1);

	fn[0] = 1.0;
	fn[127] = (float)(exp(-0.5 * dn * dn));

	for (i = 126; 1 <= i; i--)
	{
		dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
		kn[i + 1] = (uint32_t)((dn / tn) * m1);
		tn = dn;
		fn[i] = (float)(exp(-0.5 * dn * dn));
		wn[i] = (float)(dn / m1);
	}

	return;
}

/*
int main()
{
	float fn[128];
	uint32_t kn[128];
	uint32_t seed;
	float value;
	float wn[128];
	r4_nor_setup(kn, fn, wn);
	seed = 19937;

	clock_t start;
	double duration;
	start = clock();

	for (int i = 0; i < 100000 * 250 * 2; i++) {
		value = r4_nor(seed, kn, fn, wn);
	}

	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "The running time is " << duration << std::endl;

	system("pause");

	return 0.0;
}
*/
