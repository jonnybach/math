
#include <math.h>
#include "Misc.h"

/*    ''' <summary>
    ''' Gets roots of quadratic equation Ax^2 + Bx + C = 0.
    ''' Roots are sorted in descending order.
    ''' </summary>
    ''' <param name="A"></param>
    ''' <param name="B"></param>
    ''' <param name="C"></param>
    ''' <param name="Root1"></param>
    ''' <param name="Root2"></param>
    ''' <param name="DivideByZero"></param>
    ''' <param name="ImaginaryRoots"></param>
    ''' <remarks>
    ''' Returns real roots only; no imaginary numbers.
    ''' </remarks>*/
 void QuadraticRoots(const double A, const double B, const double C 
      , double *root1, double *root2, int *divideByZero, int *rootsAreImag) {

		*divideByZero = 0;
		*rootsAreImag = 0;
		
        if (A == 0) { 
            *divideByZero = 1;
            return;
        } else if ((pow(B,2) - 4 * A * C) < 0) {
            *rootsAreImag = 1;
            return;
        } else {
            *root1 = (-B + sqrt(pow(B,2) - 4 * A * C)) / (2 * A);
            *root2 = (-B - sqrt(pow(B,2) - 4 * A * C)) / (2 * A);
        }

        //Now sort the roots in descending order.
        if (root2 > root1) {
            double TempVal = *root1;
            *root1 = *root2;
            *root2 = TempVal;
        }
}

void GetArrayMin(const double x[], int *iMin, double *xMin) {

        //Reads in array X and returns the min value and its location within the array.
        //Handles zero-based arrays as well as one-based arrays.
        
        int Nx = sizeof (x) / sizeof *(x);
       
       	double xMinTmp = 99999999999999.9;
       	int i, iMinTmp;
        for (i = 0; i < Nx; i++) {
            if (x[i] < xMinTmp) {
                xMinTmp = x[i];
                iMinTmp = i;
        	}
        }
        *xMin = xMinTmp;
        *iMin = iMinTmp;  
}

void GetArrayMax(const double x[], int *iMax, double *xMax) {

        //Reads in array X and returns the max value and its location within the array.
        //Handles zero-based arrays as well as one-based arrays.
        int Nx = sizeof (x) / sizeof *(x);

       	double xMaxTmp = -99999999999999.9;
       	int i, iMaxTmp;
        for (i = 0; i < Nx; i++) {
            if (x[i] > xMaxTmp) {
                xMaxTmp = x[i];
                iMaxTmp = i;
        	}
        }
        *xMax = xMaxTmp;
        *iMax = iMaxTmp;
}

double GetArrayAve(const double x[]) {
	    
	    int Nx = sizeof (x) / sizeof *(x);
	    
	    double ave = 0;

		int i;
        for (i = 0; i < Nx; i++) {
            ave += x[i];
        }
		
        return ave / Nx;
}

double CalcSign(const double x, const double y) {
	// this routine is convoluted but is used in the Numerical recipies sorting functions
	// consider moving into there and making it static
    double s = 0;
    if (y >= 0) {
        s = fabs(x);
	} else if (y < 0) {
        s = -fabs(x);
    }
    return s;
}

int AreEqual(const double val1, const double val2, const double tol) {
	int isEqual = 0;
    if (fabs(val1 - val2) <= tol) {
        isEqual = 1;
    } else {
        isEqual = 0;
    }
    return isEqual;
}

double PolyFit(const double x, const double A[]) {
    //Evaluates a polynomial of the form:
    //   y=A(0) + A(1)*x + A(2)*x^2 + A(3)*x^3 ... + A[N]*x^(N)
    //Inputs:
    //x = independent variable
    //A[] = Polynomial coefficients
	int Na = sizeof (A) / sizeof *(A);
    double yEval = 0;
    int i;
    for (i = 0; i < Na; i++) {
        yEval += A[i] * pow(x,i);
    }
    return yEval;
}

double Factorial(const int n) {
	//calculate the factorial of an integer
	double fact = 1;
	int nTmp = n;
	while (nTmp > 1) { //NOTE: no need to multiply fact by 1, just adds extra loop call
		fact *= nTmp;
		nTmp --;
	}
	return fact;
}

double LinearZero(const double y1, const double y2, const double x1, const double x2) {
        //eqn. of form y=m*x+b
        double m = (y2 - y1) / (x2 - x1);
        double b = y1 - m * x1;
        //zero at x=-b/m
        return -b / m;
}
