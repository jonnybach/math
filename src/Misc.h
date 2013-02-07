#ifndef MISC_
#define MISC_

extern void QuadraticRoots(const double A, const double B, const double C 
      , double *root1, double *root2, int *divideByZero, int *rootsAreImag);
      
extern void GetArrayMin(const double x[], int *iMin, double *xMin);
extern void GetArrayMax(const double x[], int *iMax, double *xMax);
extern double GetArrayAve(const double x[]);

extern double CalcSign(const double x, const double y);

extern int AreEqual(const double val1, const double val2, const double tol);

extern double PolyFit(const double x, const double A[]);

extern double Factorial(const int n);

extern double LinearZero(const double y1, const double y2, const double x1, const double x2);

#endif /*MISC_*/
