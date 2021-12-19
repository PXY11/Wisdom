/*
** Author: Wu Chen
** 2017/06/27
*/
#ifndef _ROOT_FINDING
#define _ROOT_FINDING

#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <limits>


bool RootBracketing(std::function<double(double)> f, double &a, double &b, double min, double max);
bool RootBracketingFP(double(*fp)(double), double &a, double &b);
template<typename T> 
bool RootBracketingT(T f, double &a, double &b);
double rfbisect(std::function<double(double)> f, double a, double b, double tol);
double rfsecant(std::function<double(double)> f, double a, double b, double tol);
double rffalsi(std::function<double(double)> f, double a, double b, double tol);
double rfbrent(std::function<double(double)> f, double a, double b, double tol);

#endif
