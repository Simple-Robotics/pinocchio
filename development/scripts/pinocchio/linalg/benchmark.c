#include <stdio.h>
#include "inv_A_fun3_codegen.h"
#include "A_inv_chol_fun3_codegen.h"
#include "U_fun3_codegen.h"
#include <time.h>

#define casadi_real double
#define casadi_int long long int

int f(const casadi_real ** arg, casadi_real ** res, casadi_int * iw, casadi_real * w, int mem);

int main()
{ //
  double x[50] = {
    11.97623514, 0.5582902,   -2.39658057, -3.75451444, 2.84949994,  -4.9783341,
    0.5582902,   5.93341693,  -3.76665704, 0.18569052,  3.08882618,  -2.88752941,
    -2.39658057, -3.76665704, 6.02263329,  -0.95414839, -3.62364168, 2.17486679,
    -3.75451444, 0.18569052,  -0.95414839, 4.01772341,  1.63855904,  0.04778753,
    2.84949994,  3.08882618,  -3.62364168, 1.63855904,  5.36739647,  -4.06595882,
    -4.9783341,  -2.88752941, 2.17486679,  0.04778753,  -4.06595882, 6.69418877}; //= {9, 16};
  double y[50];                                                                   // = {1};
  double z[50];
  double z2[50];
  double z3[50];
  const double * arg[2] = {x, y};
  double * res[1] = {z};
  double * res2[1] = {z2};
  double * res3[1] = {z3};
  double w[60];

  int N_steps = 10000000;
  // Time the loop below
  clock_t start = clock();
  for (int i = 0; i < N_steps; i++)
  {
    inv_A_fun(arg, res, 0, w, 0);
    x[0] += 1e-9;
  }
  clock_t end = clock();
  double time_spent = (double)(end - start) / CLOCKS_PER_SEC / N_steps * 1e9;
  printf("Time spent in CasADi inverse: %f nanoseconds\n", time_spent);
  // Print the result

  start = clock();
  for (int i = 0; i < N_steps; i++)
  {
    A_inv_chol_fun(arg, res2, 0, w, 0);
    x[0] += 1e-9;
  }
  end = clock();
  time_spent = (double)(end - start) / CLOCKS_PER_SEC / N_steps * 1e9;
  printf("Time spent in LTL inverse: %f nanoseconds\n", time_spent);

  start = clock();
  for (int i = 0; i < N_steps; i++)
  {
    U_fun(arg, res, 0, w, 0);
    x[0] += 1e-9;
  }
  end = clock();
  time_spent = (double)(end - start) / CLOCKS_PER_SEC / N_steps * 1e9;
  printf("Time spent in LTL: %f nanoseconds\n", time_spent);

  return 0;
}
