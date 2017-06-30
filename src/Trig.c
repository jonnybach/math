#include <math.h>
#include "Trig.h"
#include "MathConstants.h"

double Sec(const double x_rads) {
        //Secant
        if (cos(x_rads) != 0) {
            return 1 / cos(x_rads);
        } else {
            return -9999999999999;
        }
}

double Cosec(const double x_rads) {
        //Cosecant
        if (sin(x_rads) != 0) {
            return 1 / sin(x_rads);
        } else {
            return -9999999999999;
        }
}

double Cotan(const double x_rads) {
        //Cotangent
        if (tan(x_rads) != 0) {
            return 1 / tan(x_rads);
        } else {
            return -9999999999999;
        }
}

double ArcSec(const double x_rads) {
        //Inverse Secant
        return atan(x_rads / sqrt(x_rads * x_rads - 1)) + Sign(x_rads - 1) * (PI / 2);
}

double ArcCosec(const double x_rads) {
        //Inverse Cosecant
        return atan(x_rads / sqrt(x_rads * x_rads - 1)) + Sign(x_rads - 1) * (PI / 2);
}

double ArcCotan(const double x_rads) {
        //Inverse Cotangent
        return atan(x_rads) + (PI / 2);
}

double Csch(const double x_rads) {
        //Calculates hyperbolic cosecant of x_rads
        return 1 / sinh(x_rads);
}

double Sech(const double x_rads) {
        //Calculates hyperbolic secant of x_rads
        return 1 / cosh(x_rads);
}

double Coth(const double x_rads) {
        //Calculates hyperbolic cotangent of x_rads
        return 1 / tanh(x_rads);
}

int Sign(const double x) {
	return x / fabs(x);
}
