#pragma once
/*
** This module is modified by Wu Chen 2017/07/03
*/

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

#define _MAX(a,b) (((a)>(b))?(a):(b))
#define _MIN(a,b) (((a)<(b))?(a):(b))

struct GridInfo {
	std::vector<double> Grid;
	std::vector<double> ds;
};

void Thomas_Algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c,
	const std::vector<double>& d, std::vector<double>& result);

void my_dot_product(const std::vector<double>& a, const std::vector<double>& b,
	const std::vector<double>& c, const std::vector<double>& v, std::vector<double>& result);

std::vector<double> element_prod(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> element_divide(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> abs_diff(const std::vector<double>& a, const std::vector<double>& b);
double maxinvec(const std::vector<double>& a);

std::vector<double> operator*(double a, const std::vector<double>& b);
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator+(double a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
bool operator==(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> compare(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> compare_equal(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> numvecCompare(double a, const std::vector<double>& b);
std::vector<double> vec_multiplyadd(double a, double d, const std::vector<double>& b);
std::vector<double> vecMAX(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> vecMIN(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> numvecMAX(double b, const std::vector<double>& a);
std::vector<double> numvecMIN(double b, const std::vector<double>& a);
void printvector(const std::vector<double> a);

std::vector<double> absvec(std::vector<double>& a);
double sumVec(const std::vector<double>& a);

GridInfo generate_grid_conv(double Smin, double Smax, int Ns, const std::vector<double>& targetPoint);
void printGridInfo(const GridInfo& a);