//
//
//                          PayOff2.cpp
//
//


#include <PayOff2.h>
#include <minmax.h>

PayOffCall::PayOffCall(double Strike_) : Strike(Strike_)
{
}

/*
有 const 修饰的成员函数
（指 const 放在函数参数表的后面，
而不是在函数前面或者参数表内），
只能读取数据成员，不能改变数据成员；
没有 const 修饰的成员函数，对数据成员则是可读可写的。
*/
double PayOffCall::operator () (double Spot) const
{
    return max(Spot-Strike,0.0);
}


PayOffPut::PayOffPut(double Strike_) : Strike(Strike_)
{
}


double PayOffPut::operator () (double Spot) const
{
    return max(Strike-Spot,0.0);
}



/*
 *
 * Copyright (c) 2002
 * Mark Joshi
 *
 * Permission to use, copy, modify, distribute and sell this
 * software for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * Mark Joshi makes no representations about the
 * suitability of this software for any purpose. It is provided
 * "as is" without express or implied warranty.
*/

