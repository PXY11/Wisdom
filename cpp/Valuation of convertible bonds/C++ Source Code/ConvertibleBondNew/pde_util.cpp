/*
** This module is modified by Wu Chen 2017/07/03
*/

#include "pde_util.h"

void printvector(std::vector<double> a) {
	std::cout << "vector length:" << a.size() << std::endl;
	std::cout.precision(18);
	for (int i = 0; i < a.size(); i++) {
		std::cout << a[i] << std::endl;
	}
}

// modified by Wu Chen for MATLAB index
void my_dot_product(const std::vector<double>& a, const std::vector<double>& b,
	const std::vector<double>& c, const std::vector<double>& v, std::vector<double>& result)
{
	int n = a.size();
	result[0] = b[0] * v[0] + c[0] * v[1];
	for (int i = 1; i < n - 1; i++)
		result[i] = a[i - 1] * v[i - 1] + b[i] * v[i] + c[i + 1] * v[i + 1];
	result[n - 1] = a[n - 2] * v[n - 2] + b[n - 1] * v[n - 1];
}

// modified by Wu Chen for MATLAB index
void Thomas_Algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c,
	const std::vector<double>& d, std::vector<double>& result)
{
	int n = a.size();
	std::vector<double> c_bar(n - 1);
	std::vector<double> d_bar(n);
	c_bar[0] = c[0] / b[1];
	for (int i = 1; i < n - 1; i++)
		c_bar[i] = c[i + 1] / (b[i] - c_bar[i - 1] * a[i - 1]);
	d_bar[0] = d[0] / b[0];
	for (int i = 1; i<n; i++)
		d_bar[i] = (d[i] - d_bar[i - 1] * a[i - 1]) / (b[i] - c_bar[i - 1] * a[i - 1]);
	result[n - 1] = d_bar[n - 1];
	for (int i = n - 2; i >= 0; i--)
		result[i] = d_bar[i] - c_bar[i] * result[i + 1];
}

std::vector<double> element_prod(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a[i] * b[i];
	return c;
}

std::vector<double> element_divide(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a[i] / b[i];
	return c;
}

std::vector<double> abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = std::abs(a[i] - b[i]);
	return c;
}

std::vector<double> vec_multiplyadd(double a, double d, const std::vector<double>& b)
{
	int n = b.size();
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a * b[i] + d;
	return c;
}

std::vector<double> operator*(double a, const std::vector<double>& b)
{
	int n = b.size();
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a * b[i];
	return c;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a[i] + b[i];
	return c;
}

std::vector<double> operator+(double a, const std::vector<double>& b)
{
	int m = b.size();
	std::vector<double> c(m);
	for (int i = 0; i < m; i++)
		c[i] = a + b[i];
	return c;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = a[i] - b[i];
	return c;
}

bool operator==(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	bool flag = true;
	for (int i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			flag = false;
			break;
		}
	}
	return flag;
}

std::vector<double> compare(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		if (a[i] < b[i]) {
			c[i] = 1.0;
		}
		else {
			c[i] = 0.0;
		}
	return c;
}

std::vector<double> compare_equal(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		if (a[i] <= b[i]) {
			c[i] = 1.0;
		}
		else {
			c[i] = 0.0;
		}
		return c;
}

std::vector<double> numvecCompare(double a, const std::vector<double>& b)
{
	int m = b.size();
	std::vector<double> c(m);
	for (int i = 0; i < m; i++)
		if (a < b[i]) {
			c[i] = 1.0;
		}
		else {
			c[i] = 0.0;
		}
		return c;
}

std::vector<double> absvec(std::vector<double>& a)
{
	int n = a.size();
	for (int i = 0; i < n; i++)
		a[i] = std::abs(a[i]);
	return a;
}

double maxinvec(const std::vector<double>& a)
{
	int n = a.size();
	double temp = a[0];
	for (int i = 1; i < n; i++) {
		if (a[i] > temp) {
			temp = a[i];
		}
	}
	return temp;
}

std::vector<double> vecMAX(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = _MAX(a[i], b[i]);
	return c;
}

std::vector<double> vecMIN(const std::vector<double>& a, const std::vector<double>& b)
{
	int n = a.size();
	int m = b.size();
	assert(n == m);
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = _MIN(a[i], b[i]);
	return c;
}

std::vector<double> numvecMAX(double b, const std::vector<double>& a)
{
	int n = a.size();
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = _MAX(a[i], b);
	return c;
}

std::vector<double> numvecMIN(double b, const std::vector<double>& a)
{
	int n = a.size();
	std::vector<double> c(n);
	for (int i = 0; i < n; i++)
		c[i] = _MIN(a[i], b);
	return c;
}

double sumVec(const std::vector<double>& a)
{
	int n = a.size();
	double result;
	for (int i = 0; i < n; i++)
		result += a[i];
	return result;
}

// place total Ns grid and approximately ensure same ds level
GridInfo generate_grid_conv(double Smin, double Smax, int Ns,const std::vector<double>& targetPoint)
{
	assert(Ns <= 0);
	GridInfo result;
	int n = targetPoint.size(), index = 0;
	std::vector<int> Nss(n + 1);
	for (int i = 0; i <= n; i++) {
		if (i == 0) {
			Nss[i] = ceil((double)Ns * (targetPoint[i] - Smin) / (Smax - Smin));
		}
		else if (i == n) {
			Nss[i] = ceil((double)Ns * (Smax - targetPoint[i - 1]) / (Smax - Smin));
		}
		else {
			Nss[i] = ceil((double)Ns * (targetPoint[i] - targetPoint[i - 1]) / (Smax - Smin));
		}
	}
	// new length of the vector
	int Ns_new = 0;
	for (int i = 0; i < n + 1; i++)
		Ns_new += Nss[i];
	// grid vector
	std::vector<double> c(Ns_new + 1); // include Smin, so we need to plus one here
	std::vector<double> _ds(n + 1);
	c[0] = Smin;
	double ds, temp;
	int index_length = 0;
	while (index <= n) {
		if (index == 0) { // first interval
			ds = (targetPoint[0] - Smin) / Nss[index];
			temp = Smin;
			_ds[index] = ds;
		}
		else if (index == n) { // last interval
			ds = (Smax - targetPoint[n - 1]) / Nss[index];
			temp = targetPoint[n - 1];
			_ds[index] = ds;
		} 
		else { // intermediate interval
			ds = (targetPoint[index] - targetPoint[index - 1]) / Nss[index];
			temp = targetPoint[index - 1];
			_ds[index] = ds;
		}
		for (int i = 1; i <= Nss[index]; i++) {
			c[i + index_length] = ds * i + temp;
		}
		index_length += Nss[index];
		index++;
	}
	result.Grid = c;
	result.ds = _ds;
	return result;
}

void printGridInfo(const GridInfo& a) {
	for (int i = 0; i < a.Grid.size(); i++) {
		std::cout << a.Grid[i] << " ";
	}
	std::cout << " " << std::endl;
	std::cout << "ds vector:" << std::endl;
	for (int i = 0; i < a.ds.size(); i++) {
		std::cout << a.ds[i] << " ";
	}
}

