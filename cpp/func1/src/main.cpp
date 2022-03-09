/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.cpp
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "func1.h"
#include "func1_terminate.h"
#include<iostream>
using namespace std;
/* Function Declarations */
static void argInit_1x4_real_T(double result[4]);
static double argInit_real_T();
static void main_func1();

/* Function Definitions */
static void argInit_1x4_real_T(double result[4])
{
  double result_tmp_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp_tmp = argInit_real_T();
  result[0] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[2] = result_tmp_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[3] = argInit_real_T();
}

static double argInit_real_T()
{
  return 0.0;
}

static void main_func1()
{
  double m_tmp_tmp_tmp_tmp;
  double Q_tmp_tmp[4];
  double dv[4];
  double b_Q_tmp_tmp[4];
  double CC[16];

  /* Initialize function 'func1' input arguments. */
  m_tmp_tmp_tmp_tmp = argInit_real_T();

  /* Initialize function input argument 'Q'. */
  argInit_1x4_real_T(Q_tmp_tmp);

  /* Initialize function input argument 'E'. */
  /* Initialize function input argument 'u'. */
  /* Initialize function input argument 'C0'. */
  /* Call the entry-point 'func1'. */
  argInit_1x4_real_T(dv);
  b_Q_tmp_tmp[0] = Q_tmp_tmp[0];
  b_Q_tmp_tmp[1] = Q_tmp_tmp[1];
  b_Q_tmp_tmp[2] = Q_tmp_tmp[2];
  b_Q_tmp_tmp[3] = Q_tmp_tmp[3];
  func1(m_tmp_tmp_tmp_tmp, m_tmp_tmp_tmp_tmp, m_tmp_tmp_tmp_tmp,
        m_tmp_tmp_tmp_tmp, m_tmp_tmp_tmp_tmp, b_Q_tmp_tmp, Q_tmp_tmp, Q_tmp_tmp,
        dv, CC);
   cout<<"in main_func1"<<endl;
   for(int i=0;i<16;i++)
   {
      // cout<<CC[i]<<endl;
   }
   cout<<"main_func1 end"<<endl;
}

int main(int, const char * const [])
{
  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_func1();

  /* Terminate the application.
     You do not need to do this more than one time. */
  func1_terminate();
  cout<<endl<<"over";
  return 0;
}

/* End of code generation (main.cpp) */
