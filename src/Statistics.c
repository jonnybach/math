
#include <math.h>
#include "Statistics.h"
#include "Misc.h"

double PercentError(const double Experimental, const double Theoretical) {
	//Calculates percent error between experimental and theoretical values
	double Err;
	if ((Theoretical == 0) && (Experimental == 0)) {
		Err = 0.0;
	} else if ((Theoretical == 0) && (Experimental != 0)) {
	    Err = 99999999999999999999999999.9; 
	} else {
	    Err = (Experimental - Theoretical) / Theoretical;
	}
	return Err;
}

double PercentDiff(const double x1, const double x2) {
    //Calculates percent difference between two experimental values
    double Diff;
    if ((x1 + x2) == 0) {
        Diff = 99999999999999999999999999.9; 
    } else {
        Diff = 2 * fabs(x1 - x2) / (x1 + x2);
    }
    return Diff;
}

double Variance(const double A[]) {
	int Na = sizeof (A) / sizeof *(A);
    double Avg = GetArrayAve(A);
    double Var = 0;
   	int i;
    for (i = 0; i < Na; i++) {
        Var += pow((A[i] - Avg),2);
    }
    Var /= Na;
    return Var;
}

double StdDeviation(const double A[]) {
        double StdDev = sqrt(Variance(A));
        return StdDev;
}
