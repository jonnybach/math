#ifndef ROOTFINDING_H_
#define ROOTFINDING_H_

#include <stdarg.h>

// Brent function pointer
typedef double(*BrentFunc)(double, ...);

// Brent root finding function
double Brent(const double x1, const double x2, const double tol
      , const BrentFunc Func, const double funcArgs[]);

// Brent helper function
double BrentSign(const double x, const double y);

#endif /*ROOTFINDING_H_*/
