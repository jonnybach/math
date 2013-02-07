
#include <math.h>
#include "Bessel.h"
#include "MathConstants.h"
#include "Misc.h"

double ModBessFirstKind_Inx(const double x, const int n, const double tol) {
        //computes the modified bessel function of the first kind using series and conv tol
        //x is value for function(funct arg), and n is the order of the function (subscript)
        //using formulas from wolframs math world
        //function has been checked against the excel version to be accurate within fractions of %

        //start summing the series term until the change is less than the tolerence
        //also limit this fuction to 200 terms for speed
        double Inx, sum, Kval;
        int k = 0;
        sum = 0;
        do {
            Kval = pow((0.25 * pow(x,2)), k) / (Factorial(k) * Factorial(n + k));
            if (Kval < tol) {
                //if this term is less than conv you are done
                break;
            } else { //need to add this term to the sum and inc k
                sum += Kval;
                k += 1;
            }
        } while (k <= 200);
        //now multiply by the sum by the first term
        Inx = sum * pow((0.5 * x),n);

        return Inx;
}

double ModBessSecondKind_Knx(const double x, const int n_in) {
        //computes the modified bessel function of the second kind using numerical method
        //x is value for function(funct arg), and n is the order of the function (subscript)
        //based on open source numerical fortran code found on the internet
        //function has been checked against the excel version to be accurate within fractions of %
        
        double eps = 0.000000000001;

        double bignum = 1000000000000.0;
        double eul = 0.577215664901533;
        
       	double k, kf, nk1f, nkf, zn, t, s, z0, z, ans, fn, pn, pk, zmn, tlg, tox;
		double Knx;
        int i;
                
        //positive only (symetric function)
        int n = n_in;
       	if (n < 0) n = -n;
        
        if (x <= 9.55) {
            ans = 0;
            z0 = 0.25 * pow(x,2);
            fn = 1;
            pn = 0;
            zmn = 1;
            tox = 2 / x;
            if (n > 0) {
                pn = -eul;
                k = 1;
                for (i = 0; i < n; i++) { 
                	pn = pn + 1 / k;
                    k = k + 1;
                    fn = fn * k;
                }
                zmn = tox;
                if (n == 1) {
                    ans = 1 / x;
                } else {
                    nk1f = fn / n;
                    kf = 1;
                    s = nk1f;
                    z = -z0;
                    zn = 1;
                    for (i = 0; i < n; i++) {
                        nk1f = nk1f / (n - i);
                        kf = kf * i;
                        zn = zn * z;
                        t = nk1f * zn / kf;
                        s = s + t;
                        zmn = zmn * tox;
                    }
                    s = s * 0.5;
                    t = fabs(s);
                    ans = s * zmn;
                }
            }
            tlg = 2 * log(0.5 * x);
            pk = -eul;
            if (n == 0) {
                pn = pk;
                t = 1;
            } else {
                pn = pn + 1 / n;
                t = 1 / fn;
            }
            s = (pk + pn - tlg) * t;
            k = 1;
            do {
                t = t * (z0 / (k * (k + n)));
                pk = pk + 1 / k;
                pn = pn + 1 / (k + n);
                s = s + (pk + pn - tlg) * t;
                k = k + 1;
            } while (fabs(t / s) > eps);
            s = 0.5 * s / zmn;
            if (n%2 != 0) {
                s = -s;
            }
            Knx = ans + s;
		} else if (x > log(bignum)) {
		            Knx = 0;
		} else
            k = n;
            pn = 4 * k * k;
            pk = 1;
            z0 = 8 * x;
            fn = 1;
            t = 1;
            s = t;
            nkf = bignum;
            i = 0;
            do {
                z = pn - pk * pk;
                t = t * z / (fn * z0);
                nk1f = fabs(t);
                if ((i >= n) && (nk1f > nkf)) {
                    break;
                }
                nkf = nk1f;
                s = s + t;
                fn = fn + 1;
                pk = pk + 2;
                i = i + 1;
            } while (fabs(t / s) > eps);
        Knx = exp(-x) * (PI / sqrt(2 * x)) * s;
        return Knx;
}
