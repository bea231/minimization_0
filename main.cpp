/*
* Sergeev Artemiy, 33601/2
* 1D minimization methods.
* Golden section search and Fibonacci search methods.
*/

#include <stdio.h>
#include <math.h>

/* (sqrt(5) - 1) / 2 constant value */
#define PHI1 0.61803398874989484820458683436564
/* (3 - sqrt(5)) / 2 constant value */
#define PHI2 0.38196601125010515179541316563436
/* Exponent constant */
#define EXP 2.718281828459045235360287471352663

/* Some segment type definition (a < b) */
struct segment_t
{
  double a, b;

  segment_t( void ) : a(0), b(0)
  {
  }

  segment_t( double _a, double _b) : a(_a), b(_b)
  {
  }

  segment_t( const segment_t &seg ) : a(seg.a), b(seg.b)
  {
  }
};

/* Pointer to 1D function type definition */
typedef double (*func_ptr)( double x );

double objFunc( double x )
{
  return x * x + 2 * (x * log10(x / EXP) - 2);
}

double getFirstLambda( segment_t &endpoints, double epsilon )
{
  double value = (endpoints.b - endpoints.a) / epsilon;
  unsigned int fib_1 = 1, fib_2 = 1, tmp;

  while (fib_2 <= value)
  {
    tmp = fib_2;
    fib_2 += fib_1;
    fib_1 = tmp;
  }
  return (double)fib_1 / (double)fib_2;
}

segment_t fibonacciSearch( func_ptr f, segment_t &endpoints, double epsilon, unsigned int *callCount = NULL )
{
  double delta_k = endpoints.b - endpoints.a, 
         delta_k1 = getFirstLambda(endpoints, epsilon) * delta_k,
         delta_k2 = delta_k - delta_k1,
         y = endpoints.a + delta_k2,
         z = endpoints.b - delta_k2,
         fy = f(y),
         fz = f(z);
  segment_t cur_seg(endpoints);
  unsigned int calls = 2;

  while (delta_k >= epsilon)
  {
    delta_k = delta_k1;
    delta_k1 = delta_k2;
    delta_k2 = delta_k - delta_k1;
    if (fy <= fz)
    {
      cur_seg.b = z;
      z = y;
      y = cur_seg.a + delta_k2;
      fy = f(y);
    }
    else
    {
      cur_seg.a = y;
      y = z;
      z = cur_seg.b - delta_k2;
      fz = f(z);
    }
    ++calls;
  }
  if (callCount)
    *callCount = calls;

  return cur_seg;
}

segment_t goldenSectionSearch( func_ptr f, segment_t &endpoints, double epsilon, unsigned int *callCount = NULL )
{
  double y = endpoints.a + PHI2 * (endpoints.b - endpoints.a),
         z = endpoints.b - PHI2 * (endpoints.b - endpoints.a), 
         fy = f(y), fz = f(z);
  segment_t cur_seg(endpoints);
  unsigned int calls = 0;

  while (cur_seg.b - cur_seg.a >= epsilon)
  {
    if (fy <= fz)
    {
      cur_seg.b = z;
      z = y;
      fz = fy;
      y = cur_seg.a + PHI2 * (cur_seg.b - cur_seg.a);
      fy = f(y);
    }
    else
    {
      cur_seg.a = y;
      y = z;
      fy = fz;
      z = cur_seg.b - PHI2 * (cur_seg.b - cur_seg.a);
      fz = f(z);
    }
    ++calls;
  }
  if (callCount)
    *callCount = calls;
  return cur_seg;
}

int main( void )
{
  segment_t endpoints(1.5, 2);
  unsigned int calls;

  /* Golden section search */
  printf("Golden section search method\n");
  printf("---------------------------------\n");
  segment_t result = goldenSectionSearch(objFunc, endpoints, 0.1, &calls);
  printf("epsilon      | 0.1\n");
  printf("[a; b]       | [%.2g; %.2g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  printf("---------------------------------\n");
  result = goldenSectionSearch(objFunc, endpoints, 0.01, &calls);
  printf("epsilon      | 0.01\n");
  printf("[a; b]       | [%.3g; %.3g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  printf("---------------------------------\n");
  result = goldenSectionSearch(objFunc, endpoints, 0.001, &calls);
  printf("epsilon      | 0.001\n");
  printf("[a; b]       | [%.4g; %.4g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  /* Fibonnaci search */
  printf("\n\nFibonacci search method\n");
  printf("---------------------------------\n");
  result = fibonacciSearch(objFunc, endpoints, 0.1, &calls);
  printf("epsilon      | 0.1\n");
  printf("[a; b]       | [%.2g; %.2g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  printf("---------------------------------\n");
  result = fibonacciSearch(objFunc, endpoints, 0.01, &calls);
  printf("epsilon      | 0.01\n");
  printf("[a; b]       | [%.3g; %.3g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  printf("---------------------------------\n");
  result = fibonacciSearch(objFunc, endpoints, 0.001, &calls);
  printf("epsilon      | 0.001\n");
  printf("[a; b]       | [%.4g; %.4g]\n", result.a, result.b);
  printf("[f(a); f(b)] | [%g; %g]\n", objFunc(result.a), objFunc(result.b));
  printf("f(x) calls   | %u\n", calls);

  return 0;
}