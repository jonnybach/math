
#include <math.h>
#include "RootFinding.h"

/*
//"Numerical Recipies Function Pointers"

//Public Delegate Function BrentFunc(ByVal x As Double, ByRef Obj As Object) As Double


//Public Delegate Function FuncTrap(ByVal x As Double, ByVal Obj As Object) As Double
double (*FuncTrap)(double, *void);

//Public Delegate Function RK4Derivs(ByVal x As Double, ByRef y() As Double, ByRef dydx() As Double) As Double
double (*RK4Derivs)(double, double [], double dydx[]);

//Public Delegate Function GenericFunction(ByVal Obj As Object) As Double
double (*GenericFunction)(void);

//Public Delegate Sub NewtonRaphSub(ByRef Obj As Object, ByRef Xvals() As Double, ByRef Fvals() As Double, ByVal Jacob(,) As Double)
double (*NewtonRaphSub)(void, double [], double[], int, int, double[][])
*/


/*    ''' <summary>
    ''' Fourth-order Runge-Kutta program
    ''' Source: "Numerical Recipes in Fortran", 2nd edition by Press, Teukolsky, Vetterling, 
    ''' and Flannery.
    ''' Given the values for the variables y(1:N) and their derivatives dydx(1:N) known at x, 
    ''' use the fourth-order Runge-Kutta method to advance the solution over an interval h 
    ''' and return the incremented variables as yout(1:N), which need not be a distinct array 
    ''' from y.  The user supplies the subroutine derivs(x,y,dydx), which returns derivatives 
    ''' dydx at x.
    ''' </summary>
    ''' <param name="y"></param>
    ''' <param name="dydx"></param>
    ''' <param name="x"></param>
    ''' <param name="h"></param>
    ''' <param name="yout"></param>
    ''' <param name="RK4Func"></param>
    ''' <remarks></remarks>
 void RungeKutta4(ByRef y() As Double, ByRef dydx() As Double, ByVal x As Double, _
      ByVal h As Double, ByRef yout() As Double, ByVal RK4Func As RK4Derivs)

        Dim dym(), dyt(), yt() As Double
        Dim hh, h6, xh As Double
        Dim i, N As Integer

        N = UBound(y)
        ReDim dym[N], dyt[N], yt[N]

        hh = h / 2
        h6 = h / 6
        xh = x + hh
        ' First step
        For i = 1 To N
            yt(i) = y(i) + hh * dydx(i)
        Next i
        ' Second step
        Call RK4Func(xh, yt, dyt)
        For i = 1 To N
            yt(i) = y(i) + hh * dyt(i)
        Next i
        ' Third step
        Call RK4Func(xh, yt, dym)
        For i = 1 To N
            yt(i) = y(i) + h * dym(i)
            dym(i) = dyt(i) + dym(i)
        Next i
        ' Fourth step
        Call RK4Func(x + h, yt, dyt)
        For i = 1 To N 'Accumulate increments with proper weights
            yout(i) = y(i) + h6 * (dydx(i) + dyt(i) + 2 * dym(i))
        Next i

}*/

/*    ''' <summary>
    ''' Uses trapezoidal rule to integrate a function between start point and end point.
    ''' Adjust NumSteps if more accuracy is required
    ''' </summary>
    ''' <param name="Obj"></param>
    ''' <param name="Func"></param>
    ''' <param name="EndPoint"></param>
    ''' <param name="StartPoint"></param>
    ''' <param name="NumSteps"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
double IntTrapRule(ByRef Obj As Object, ByRef Func As FuncTrap, _
      ByVal EndPoint As Double, Optional ByVal StartPoint As Double = 0, _
      Optional ByVal NumSteps As Integer = 1000) {

        Dim Delta As Double, x1 As Double, x2 As Double, f1 As Double, f2 As Double, Area As Double
        Dim i As Integer

        Delta = (EndPoint - StartPoint) / NumSteps
        x1 = StartPoint
        x2 = x1 + Delta
        f1 = Func(x1, Obj)
        f2 = Func(x2, Obj)
        Area = Delta * (f1 + f2) / 2
        For i = 1 To NumSteps - 1
            x1 = x2
            x2 = x1 + Delta
            f1 = Func(x1, Obj)
            f2 = Func(x2, Obj)
            Area = Area + Delta * (f1 + f2) / 2
        Next i

        return Area

}*/

/*    ''' <summary>
    ''' Using Brent's method, find the root of a function func known to lie between x1 and x2
    ''' </summary>
    ''' <param name="x1">Minimum x value</param>
    ''' <param name="x2">Maximum x value</param>
    ''' <param name="tol">Convergence tolerance</param>
    ''' <param name="Func">Function for which to obatin root</param>
    ''' <param name="Obj">Arraylist of parameters to be used in analysis.</param>
    ''' <Returns></Returns>
    ''' <remarks>
    ''' The Obj arraylist, if used, must be constructed prior to calling the Brent method.
    ''' In the function Func, the arraylist must be deconstructed into a list of objects.
    ''' </remarks>*/

inline int max(double a, double b)
{ return a > b ? a : b; }

inline int min(double a, double b)
{ return a < b ? a : b; }

double BrentSign(const double x, const double y) {
    double rslt = 0;
    if (y >= 0) {
        rslt = fabs(x);
    } else if (y < 0) {
        rslt = -fabs(x);
    }
    return rslt;
}

double Brent(const double x1, const double x2, const double tol
      , const BrentFunc Func, const double funcArgs[]) {

        int j;
        int jMax = 100;
        double zEps = 0.00000003;
        double a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, min1, min2;

        c = d = e = 0;

        a = x1;
        b = x2;
        fa = Func(a, funcArgs);
        if (fabs(fa) < tol) return a;  //a bound was the root
        fb = Func(b, funcArgs);
        if (fabs(fb) < tol) return b; //b bound was the root
        //are these bounds ok
        if ((fa * fb) >=1) {
        	//root not bounded because product of function evals using guesss is positive
        	// (either 2 positive or 2 negative numbers
        	return -9999999999999999999999.999;
        }
        
        fc = fb;
        for (j=0; j <= jMax; j++) {
            
            if (((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0))) {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
            
            if (fabs(fc) < fabs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            
            tol1 = 2 * zEps * fabs(b) + tol / 2; //convergence check
            xm = (c - b) / 2;
            
            if ((fabs(xm) <= tol1) || (fb == 0)) {
            	//within tolerance, found root!!
               	return b;
            }
            
            if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb))) {
                s = fb / fa;  //attempt inverse quadratic interpolation
                if (a == c) {
                    p = 2 * xm * s;
                    q = 1 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
                    q = (q - 1) * (r - 1) * (s - 1);
                }
                
                if (p > 0) q = -q; //check whether in bounds
                
                p = fabs(p);
                min1 = 3 * xm * q - fabs(tol1 * q);
                min2 = fabs(e * q);
                if ((2 * p) < min(min1, min2)) {
                    e = d;			//accept interpolation
                    d = p / q;
                } else {
                    d = xm;			//interpolation failed use bisection
                    e = d;
                }
            } else {
                d = xm;				//bounds decreasing too slowly, use bisection
                e = d;
            }
            
            a = b;					//move last best guess to a
            fa = fb;				
            
            if (fabs(d) > tol1) {   //evaluate new trial root
                b += d;
            } else {
                b += BrentSign(tol1, xm);
            }

            fb = Func(b, funcArgs);
            
            //ensure that rountine stops as soon as a good value is found
            if (fabs(fb) < tol) return b;
            
        } //next iteration in for loop of j
        return b;
}



/*    ''' <summary>
    ''' Given a bracketing triplet of ax, bx, and cx (such that bx is between ax and cx, and 
    ''' f(bx) is less than both f(ax) and f(ax), this routine isolates the minimum to a fractional
    ''' precision of about tol using Brent's method.
    ''' </summary>
    ''' <param name="ax"></param>
    ''' <param name="bx"></param>
    ''' <param name="cx"></param>
    ''' <param name="tol"></param>
    ''' <param name="Func"></param>
    ''' <param name="Obj"></param>
    ''' <param name="xMin"></param>
    ''' <param name="yMin"></param>
    ''' <remarks></remarks>
    Public Shared Sub BrentMinimum(ByVal ax As Double, ByVal bx As Double, _
      ByVal cx As Double, ByVal tol As Double, ByVal Func As BrentFunc, _
      ByRef Obj As Object, ByRef xMin As Double, ByVal yMin As Double)

        Dim ItMax As Integer = 100
        Dim Cgold As Double = 0.381966
        Dim zEps As Double = 0.0000000001
        Dim a, b, d, e, eTemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm As Double
        Dim iter As Integer

        a = Min(ax, cx)
        b = Max(ax, cx)
        v = bx
        w = v
        x = v
        e = 0
        fx = Func(a, Obj)
        fv = fx
        fw = fx

        For iter = 1 To ItMax
            xm = (a + b) / 2
            tol1 = tol * Abs(x) + zEps
            tol2 = 2 * tol1
            If Abs(x - xm) <= tol2 - (b - 1) / 2 Then GoTo 3
            If Abs(e) > tol1 Then
                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = w * (q - r)
                If q > 0 Then p = -p
                q = Abs(q)
                eTemp = e
                e = d
                If Abs(p) >= q * eTemp / 2 Or p <= 1 * (a - x) Or p >= q * (b - x) Then
                    GoTo 1
                End If
                d = p / q
                u = x + d
                If u - a < tol2 Or b - u < tol2 Then d = BrentSign(tol1, xm - x)
                GoTo 2
            End If

1:          If x >= xm Then
                e = a - x
            Else
                e = b - x
            End If
            d = Cgold * e
2:          If Abs(d) >= tol1 Then
                u = x + d
            Else
                u = x + BrentSign(tol1, d)
            End If
            fu = Func(u, Obj)

            If fu <= fx Then
                If u >= x Then
                    a = x
                Else
                    b = x
                End If
                v = w
                fv = fw
                w = x
                fw = fx
                x = u
                fx = fu
            Else
                If u < x Then
                    a = u
                Else
                    b = u
                End If
                If fu <= fw Or w = x Then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                ElseIf fu <= fv Or v = x Or v = w Then
                    v = u
                    fv = fu
                End If
            End If
        Next iter

3:      xMin = x
        yMin = fx

}*/

/*    ''' <summary>
    ''' Find bounds for Brent solver. This occurs where the function changes sign.
    ''' </summary>
    ''' <param name="Func">Function to evaluate</param>
    ''' <param name="xMin">Minimum bound</param>
    ''' <param name="xMax">Maximum bound</param>
    ''' <param name="NumSteps">Number of steps between xMin and xMax</param>
    ''' <param name="xLow">Lower bound, just before function crosses root</param>
    ''' <param name="xHigh">Upper bound, just after function crosses root</param>
    ''' <param name="AL">Arraylist of objects used by function Func</param>
    ''' <param name="RootExists"></param>
    ''' <param name="StartFromXmin">Set to true if starting from xMin, false if starting from xMax</param>
    ''' <remarks></remarks>
    Public Shared Sub FindBrentBounds(ByVal Func As BrentFunc, ByVal xMin As Double, _
      ByVal xMax As Double, ByVal NumSteps As Integer, ByRef xLow As Double, ByRef xHigh As Double, _
      ByRef AL As Object, Optional ByRef RootExists As Boolean = True, _
      Optional ByVal StartFromXmin As Boolean = True, Optional ByVal Tol As Double = 0, _
      Optional ByRef Root As Double = 0)

        Dim Spacing As Double = (xMax - xMin) / NumSteps
        Dim Delta, LastDelta As Double
        Dim x1, x2 As Double
        Dim i As Integer = 0

        If StartFromXmin Then
            x2 = xMin
        Else
            Spacing *= -1
            x2 = xMax
        End If

        Delta = Func(x2, AL)
        If Abs(Delta) <= Tol Then
            x1 = x2
            Root = x2
            RootExists = True
            xLow = Min(x1, x2)
            xHigh = Max(x1, x2)
            Exit Sub
        End If

        Do
            x1 = x2
            x2 += Spacing
            LastDelta = Delta
            Delta = Func(x2, AL)
            If Abs(Delta) <= Tol Then
                x1 = x2
                Root = x2
                RootExists = True
                xLow = Min(x1, x2)
                xHigh = Max(x1, x2)
                Exit Sub
            End If
            i += 1
        Loop Until Delta / LastDelta < 0 Or i = NumSteps

        xLow = Min(x1, x2)
        xHigh = Max(x1, x2)

        If Delta / LastDelta < 0 Then
            RootExists = True
        Else
            RootExists = False
        End If

}*/

/*    ''' <summary>
    ''' Root finding solver using the Ridder Method.
    ''' </summary>
    ''' <param name="LowerXBound">Negative bound of x.</param>
    ''' <param name="UpperXBound">Positive bound of x.</param>
    ''' <param name="RidderFunction">Delegate of the function that is being evaluted.</param>
    ''' <returns>Returns the x value that corresponds to the root of the function.</returns>
    ''' <remarks></remarks>
double SolveRootRidder(ByVal LowerXBound As Double, _
      ByVal UpperXBound As Double, ByVal RidderFunction As GenericFunction) {

        Dim i, Imax As Integer
        Dim Xlow, Xup, Xmid, Xroot, Xnew As Double
        Dim fLow, fUp, fMid, fRoot As Double
        Dim s, sign As Double
        Dim tol As Double = 0.000000001
        Imax = 60

        Xlow = LowerXBound
        Xup = UpperXBound
        fLow = RidderFunction(Xlow)
        fUp = RidderFunction(Xup)

        'initially set Xroot to a highly unlikely value
        Xroot = 1.0E+20

        If (fLow > 0 And fUp < 0) Or (fLow < 0 And fUp > 0) Then
            For i = 1 To Imax

                'calculate the function value at the midpoint
                Xmid = (Xlow + Xup) / 2
                fMid = RidderFunction(Xmid)

                s = Sqrt(fMid ^ 2 - (fLow * fUp))
                If s = 0 Then Return -9999999

                'determine a new x value
                sign = (Xlow - Xup) / Abs(Xlow - Xup)
                Xnew = Xmid + (Xmid - Xlow) * sign * Xmid / s

                'check for convergence on the x value at the root (y=0)
                If (Xnew - Xroot) < tol Then
                    Return Xroot
                Else
                    Xroot = Xnew
                    fRoot = RidderFunction(Xroot)
                    If fRoot = 0 Then Return Xroot
                End If

                If fMid * fRoot < 0 Then 'new x value has crossed zero between the midpoint and upper bound
                    Xlow = Xmid
                    fLow = fMid
                    Xup = Xroot
                    fUp = fRoot
                ElseIf fLow * fRoot < 0 Then 'new x value has crossed zero from the lower bound side
                    Xup = Xroot
                    fUp = fRoot
                ElseIf fUp * fRoot < 0 Then 'new x value has crossed zero from the upper bound side
                    Xlow = Xroot
                    fLow = fRoot
                End If

                If Abs(Xup - Xlow) < tol Then Return Xroot

            Next i
        ElseIf fLow = 0 Then
            Return Xlow
        ElseIf fUp = 0 Then
            Return Xup
        Else 'root is not bracketed by the lower and upper x bounds
            Return -9999999
        End If

}*/

/*    ''' <summary>
    ''' Given an initial guess of X(), take MaxIter Newton-Raphson steps to improve the root. Stop if the
    ''' root converges in either summed absolute increments TolX or summed absolute function values TolF.
    ''' </summary>
    ''' <remarks></remarks>
    Public Shared Sub NewtonRaphsonMult(ByVal GenSub As NewtonRaphSub, ByRef Obj As Object, _
      ByVal MaxIter As Integer, ByRef X() {, ByVal TolX As Double, ByVal TolF As Double)

        Dim N As Integer = UBound(X)
        Dim ErrX, ErrF As Double
        Dim Fjac(N, N), Fvec[N], P[N] As Double
        Dim iter As Integer = 1
        Do
            GenSub(Obj, X, Fvec, Fjac)
            ErrF = 0
            For i As Integer = 1 To N
                ErrF += Abs(Fvec(i))
            Next i

            For i As Integer = 1 To N 'Right-hand side of linear equations
                P(i) = -Fvec(i)
            Next i

            SolveMatrix(Fjac, P)

            ErrX = 0
            For i As Integer = 1 To N
                ErrX += Abs(P(i))
                X(i) += P(i)
            Next i
            iter += 1
        Loop Until ErrF <= TolF Or ErrX <= TolX Or iter = MaxIter

}*/

/*    ''' <summary>
    ''' Computes forward-difference approximation to Jacobian. On input, x() is the point at which the
    ''' Jacobian is to be evaluated, and Fvec() is the vector of function values at the point.
    ''' dF is the Jacobian array that is returned.
    ''' </summary>
    ''' <remarks></remarks>
    Public Shared Sub Jacobian(ByVal GenSub As NewtonRaphSub, ByRef Obj As Object, _
      ByRef X() As Double, ByRef Fvec() As Double, ByRef dF(,) As Double)

        Dim N As Integer = UBound(X)
        Const Eps As Double = 0.000001
        ReDim dF(N, N)
        Dim Temp, h As Double
        Dim f[N] As Double

        For j As Integer = 1 To N
            Temp = X(j)
            h = Eps * Abs(Temp)
            If h = 0 Then h = Eps
            X(j) = Temp + h
            h = X(j) - Temp
            GenSub(Obj, X, f, Nothing)
            X(j) = Temp
            For i As Integer = 1 To N
                dF(i, j) = (f(i) - Fvec(i)) / h
            Next i
        Next j

}*/

/*double InitialGuessRootFinder(ByVal Funct As BrentFunc, ByRef obj As Object
	, ByVal x1 As Double, ByVal x2 As Double, ByVal Tol As Double
	, ByRef FailedToConv As Boolean, Optional ByVal AllowBigDel As Boolean = True
	, Optional ByVal GuessOneIsBound As Boolean = False
	, Optional ByRef RootNotWithinTol As Boolean = False) {
    
        Dim x1Sign As Integer
        Dim maxiter As Integer = 100
        Dim iter As Integer = 1
        Dim Lb, Ub, xLast, yLast, xNow, yNow, newx, Bnd As Double
        Dim whileLoopiter As Integer

        FailedToConv = False

        If GuessOneIsBound Then
            'one input is bound other indicates direction
            'assume guess is never equal (bad input)
            Bnd = x1
        End If

        'run first point
        Dim y1 As Double = Funct(x1, obj)
        If Abs(y1) <= Tol Then Return x1 'converged with first guess
        If Double.IsNaN(y1) Or Double.IsInfinity(y1) Then
            FailedToConv = True
            Return x1 'prevent run time error divide by infinity
        End If

        x1Sign = y1 / Abs(y1) 'record initial sign
        Dim y2 As Double = Funct(x2, obj)
        If Abs(y2) <= Tol Then Return x2 'converged with second guess

        'Set New Variables for iterations
        xLast = x1
        yLast = y1
        xNow = x2
        yNow = y2

        'check for sign change between guess 1 and guess 2
        If x1Sign = y2 / Abs(y2) Then 'if did not cross zero
            'aproach root using linear approximation
            'exit this if: 
            '1 within tol of root
            '2 have crossed the zero axis(will use bisection with these bounds)
            '3 have exceed max number of iteration
            Do
                'update guess
                newx = LinearZero(yLast, yNow, xLast, xNow)
                If GuessOneIsBound Then
                    'need to make sure that the bound has not been crossed
                    '1.Half this change predicted until greater than LB
                    '2.Half this change predicted until less than UB
                    whileLoopiter = 1
                    While (newx - xNow) / (Bnd - xNow) >= 1
                        newx = xNow + (newx - xNow) / (Bnd - xNow) * 0.5
                        whileLoopiter += 1
                        If whileLoopiter > maxiter Then
                            FailedToConv = True
                            Exit Function
                        End If
                    End While
                End If
                If Not AllowBigDel Then
                    'dont allow more than a 100% change across 1 iteration
                    If Abs(newx - xNow) / xNow >= 1 Then
                        'Half this change predicted until less than 1
                        whileLoopiter = 1
                        While Abs(newx - xNow) / xNow >= 1
                            newx = xNow + (newx - xNow) / xNow * 0.5 * xNow
                            whileLoopiter += 1
                            If whileLoopiter > maxiter Then
                                FailedToConv = True
                                Exit Function
                            End If
                        End While
                    End If
                End If


                'update last to now
                xLast = xNow
                yLast = yNow
                xNow = newx
                'calculate function val
                yNow = Funct(xNow, obj)

                iter += 1
                If Abs(yNow) <= Tol Then Return xNow 'converged with this guess
                If x1Sign <> yNow / Abs(yNow) Then GoTo brent 'sign changes switch to brent with these bounds
            Loop Until iter >= maxiter
            If iter >= maxiter Then FailedToConv = True

        Else 'use brent solver to find the root (with the discovered bounds from sign change)
brent:      Lb = Min(xLast, xNow)
            Ub = Max(xLast, xNow)
            Return Brent(Lb, Ub, Tol, Funct, obj, FailedToConv, RootNotWithinTol)
        End If
}*/
