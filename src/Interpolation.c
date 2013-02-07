
#include <math.h>
#include "Interpolation.h"

/*''' <summary>
''' Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
''' x1 (less than) x2 (less than) .. .  (less than) xN, and given values yp1 and ypn for the first derivative of the interpolating
''' function at points 1 and n, respectively, this routine returns an array y2(1:n) of
''' length n which contains the second derivatives of the interpolating function at the tabulated
''' points xi. If yp1 and/or ypn are equal to 1 � 1030 or larger, the routine is signaled to set
''' the corresponding boundary condition for a natural spline, with zero second derivative on
''' that boundary.
''' </summary>
''' <param name="x"></param>
''' <param name="y"></param>
''' <param name="N"></param>
''' <param name="yp1"></param>
''' <param name="ypn"></param>
''' <param name="y2"></param>
''' <remarks></remarks>*/
//void spline(const int maxIndx, const double x[], const double y[], const double yp1, const double ypn, double y2[]) {
//
//		//FUNCTION TAKEN FROM NUMERICAL RECIPES BOOK, THEY HAVE STRICT LICENSE REQUIREMENTS DO NOT USE
//
//        double P, Qn, sig, Un;
//        int i, kx, k;
//
//		//int N = (sizeof (x) / sizeof *(x));
//		int N = maxIndx;
//		double u[N];
//
//        if (yp1 > 9.9E+29) { //The lower boundary condition is set to be natural
//            y2[0] = 0;
//            u[0] = 0;
//        } else { //The lower boundary condition is set to a specified first derivative.
//            y2[0] = -0.5;
//            u[0] = (3 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
//        }
//
//        for (i = 1; i < N; i++) { //This is the decomposition loop of the tridiagonal algorithm. y2 and u are used
//            //for temporary storage of the decomposed factors.
//            sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
//            P = sig * y2[i - 1] + 2;
//            y2[i] = (sig - 1) / P;
//            u[i] = (6 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
//                / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / P;
//        }
//
//        if (ypn > 9.9E+29) { //The upper boundary condition is set to be �natural�
//            Qn = 0;
//            Un = 0;
//        } else {   //The upper boundary condition is set to a specified first derivative.
//            Qn = 0.5;
//            Un = (3 / (x[N] - x[N-1])) * (ypn - (y[N] - y[N-1]) / (x[N] - x[N-1]));
//        }
//
//        y2[N] = (Un - Qn * u[N-1]) / (Qn * y2[N-1] + 1);
//
//        for (kx = 1; kx < N; kx++) { //This is the backsubstitution loop of the tridiagonal algorithm
//            k = (N-1) - (kx - 1);
//            y2[k] = y2[k] * y2[k + 1] + u[k];
//        }
//
//}
//
///*  ''' <summary>
//''' Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
//''' xai's in order), and given the array y2a(1:n), which is the output from spline above,
//''' and given a value of x, this routine returns a cubic-spline interpolated value y.
//''' </summary>
//''' <param name="xa"></param>
//''' <param name="ya"></param>
//''' <param name="y2a"></param>
//''' <param name="N"></param>
//''' <param name="x"></param>
//''' <param name="y"></param>
//''' <remarks></remarks>*/
//void splint(const int maxIndx, const double xa[], const double ya[], const double y2a[], const double x, double *y) {
//
//		//FUNCTION TAKEN FROM NUMERICAL RECIPES BOOK, THEY HAVE STRICT LICENSE REQUIREMENTS DO NOT USE
//
//		//int N = (sizeof (xa) / sizeof *(xa));
//        int N = maxIndx;
//        int k, khi, klo;
//        double A, B, h;
//
//        klo = 0;     /*We will find the right place in the table by means of bisection.
//        'This is optimal if sequential calls to this routine are at random
//        'values of x. If sequential calls are in order, and closely
//        'spaced, one would do better to store previous values of
//        'klo and khi and test if they remain appropriate on the
//        'next call.*/
//
//        khi = N;
//
//      	while ((khi - klo) > 1) {
//            k = (khi + klo) / 2;
//            if (xa[k] > x) {
//                khi = k;
//            } else {
//                klo = k;
//            }
//		}  //klo and khi now bracket the input value of x.
//
//        h = xa[khi] - xa[klo];
//        if (h == 0) {
//            *y = -1000000000000.0;
//            return;
//        }
//        A = (xa[khi] - x) / h; // 'Cubic spline polynomial is now evaluated.
//        B = (x - xa[klo]) / h;
//        *y = A * ya[klo] + B * ya[khi] + (pow(A,3) - A * y2a[klo] + (pow(B,3) - B) * y2a[khi]) * pow(h,2) / 6;
//}
//
//void Derivative(const int maxIndx, const double x[], const double y[], double *yp1, double *ypn) {
//
//		//FUNCTION TAKEN FROM NUMERICAL RECIPES BOOK, THEY HAVE STRICT LICENSE REQUIREMENTS DO NOT USE
//
//		int Nx, Ny;
//		Nx = maxIndx;
//		Ny = Nx;
//		//Nx = sizeof (x) / sizeof *(x);
//		//Ny = sizeof (y) / sizeof *(y);
//
//        //Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
//        //x1 < x2 < .. . < xN, This routine calculates values yp1 and ypn for the first derivative
//        //of the interpolating function at points 1 and n, respectively
//        if ((x[1] == x[0]) || (x[Nx] == x[Nx-1])) {
//            if (x[1] == x[0]) *yp1 = 1000000000000.0;
//            if (x[Nx] == x[Nx-1]) *ypn = 1000000000000.0;
//        } else {
//            *yp1 = (y[1] - y[0]) / (x[1] - x[0]);
//            *ypn = (y[Nx] - y[Nx-1]) / (x[Nx] - x[Nx-1]);
//        }
//}
//
///*	''' <summary>
//''' Given arrays xa(1:n) and ya(1:n) containing a tabulated function, i.e., yai = f(xai), with
//''' xa1 (less than) xa2 (less than) .. . (less than) xaN, This routine performs a cubic spline interpolation to calculate
//''' yy=f(xx)
//''' </summary>
//''' <param name="Xa"></param>
//''' <param name="Ya"></param>
//''' <param name="Xin"></param>
//''' <returns></returns>
//''' <remarks></remarks>*/
//double CubicInterp(const int numElems, const double Xa[], const double Ya[], const double Xin) {
//
//		//FUNCTION TAKEN FROM NUMERICAL RECIPES BOOK, THEY HAVE STRICT LICENSE REQUIREMENTS DO NOT USE
//
//		//int N = (sizeof Xa / sizeof Xa[0]); //doesn't work b/c Xa is a pointer to array and we know nothing of the type the pointer points to
//		int N = numElems-1;
//
//        double yp1 = 0, ypn = 0, yy = 0;
//        double ya2[N];
//        int i = 0;
//        for (i=0; i<=N; i++) {
//        	ya2[i]=0;
//        }
//
//        Derivative(N, Xa, Ya, &yp1, &ypn);
//        spline(N, Xa, Ya, yp1, ypn, &ya2[0]);
//        splint(N, Xa, Ya, ya2, Xin, &yy);
//
//        return yy;
//}

double LinearInterpUnsorted(const int numElems, const double Xaray[], const double Yaray[]
	, const double x, const int extrap) {

        //This routine interpolates within the input array
        //Note that the number of points in the array do NOT need to be known in advance, 
        //nor do the arrays have to be in a sequential order
		//NOTE that a significant performance impact could occur for this routine
		// if it is called in a loop since a lot of effort is expended in
		// searching for bounding values in the arrays.
		// Consider using LinearInterp for better performance if
		// it can be guaranteed that the arrays are sorted

        double xx, xLow, xHigh, yLow, yHigh, xMin, xMax;
        int i, iMin, iMax;

		int Nx, Ny;
		Nx = (sizeof (Xaray) / sizeof (Xaray[0]));  //doesn't work b/c Xaray and Yaray are pointers to arrays and we know nothing of the type the pointer points to
		Ny = (sizeof (Yaray) / sizeof (Yaray[0]));

        if (Nx != Ny) {
           return -99999999999;
        }

        iMin = 0;
        iMax = numElems-1;
   
        xx = x;

        //initializing
        xLow = 1000000000000.0;
        xHigh = -1000000000000.0;
        xMin = 1000000000000.0;
        xMax = -1000000000000.0;

        //Find two closest points in x array to the value xx
        for (i = iMin; i <= iMax; i++) {
            if (Xaray[i] == xx) { //No need to interpolate if we find an array value that matches x
                return Yaray[i];
            } else if ((Xaray[i] >= xx) && (Xaray[i] < xLow)) {
                xLow = Xaray[i];
                yLow = Yaray[i];
            } else if ((Xaray[i] <= xx) && (Xaray[i] > xHigh)) {
                xHigh = Xaray[i];
                yHigh = Yaray[i];
            }

            if (Xaray[i] > xMax) xMax = Xaray[i];
            if (Xaray[i] < xMin) xMin = Xaray[i];
        }

       if (xLow > xMax) { //value of x is higher than maximum point in array
            xLow = -1000000000000.0; //reset to low value --> now we'll find second-highest point in array
            xHigh = xMax;
            for (i=iMin; i<=iMax; i++) {
                if ((Xaray[i] < xHigh) && (Xaray[i] > xLow)) {
                    xLow = Xaray[i];
                    yLow = Yaray[i];
                }
            }
       } else if (xHigh < xMin) { //value of x is lower than minimum point in array
            xHigh = 1000000000000.0; //reset to high value --> now we'll find second-lowest point in array
            xLow = xMin;
            for (i = iMin; i<=iMax; i++) {
                if ((Xaray[i] < xHigh) && (Xaray[i] > xLow)) { 
                    xHigh = Xaray[i];
                    yHigh = Yaray[i];
            	}
            }
       }

		if (!extrap) {
            if (x < xMin) xx = xMin;
            if (x > xMax) xx = xMax;
		}	
		
		double yinterp;
        yinterp = Interpolate(xLow, xHigh, yLow, yHigh, xx);
		return yinterp;
}

double Interpolate(const double xlow, const double xup
		, const double ylow, const double yup, const double xinterp) {
        
        double yinterp;
        
        //Calculates value y2 at x2, using linear interpolation through points (x1,y1) and (x3,y3)
        if (xup == xlow) {
            yinterp = (yup - ylow) / 2.0;
        } else {
            yinterp = ylow + (xinterp - xlow) / (xup - xlow) * (yup - ylow);
        }
		return yinterp;
}


/*    ''' <summary>
    ''' Double linear interpolation function
    ''' Linearly interpolates on the first array of values, then interpolates on the second array
    ''' Example:
    ''' Let Var1=158 and Var2=1760
    ''' Let Var1Array = (50, 100, 150, 200, 250, 300)
    ''' Let Var2Array = (1700, 2000, 2300, 2600)
    ''' Interpolate table values between 150 and 200 on Var1, then interpolate table values
    ''' between 1700 and 2000 on Var2
    ''' </summary>
    ''' <param name="Var1"></param>
    ''' <param name="Var2"></param>
    ''' <param name="Var1Array"></param>
    ''' <param name="Var2Array"></param>
    ''' <param name="TableArray"></param>
    ''' <param name="Extrap"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>*/
double DoubleLinearInterp(double Var1, double Var2
	, double Var1Array[], double Var2Array[]
	, double **TableArray
	, int extrap) {

        int i, j;
        
        int Nvar1, Nvar2;
        Nvar1 = sizeof (Var1Array) / sizeof *(Var1Array);
        Nvar2 = sizeof (Var2Array) / sizeof *(Var2Array);

        double X2[Nvar2], X1[Nvar1];

        for (j = 0; j < Nvar2; j++) {
            for (i = 0; i < Nvar1; i++) {
                X1[i] = TableArray[i][j];
            }
            X2[j] = LinearInterpUnsorted(Nvar1, Var1Array, X1, Var1, extrap);
		}
 
        double yinterp = LinearInterpUnsorted(Nvar2, Var2Array, X2, Var2, extrap);
        return yinterp;
}

/*    ''' <summary>
    ''' Double linear interpolation function
    ''' Linearly interpolates on 2-D array
    ''' Does not extrapolate (if x or y outside range then Eval at Range Limit
    ''' Arrays must be sorted in increasing sequential order
    ''' Table is 2-D array (x,y)
    ''' </summary>
    ''' <param name="x"></param>
    ''' <param name="y"></param>
    ''' <param name="XArray"></param>
    ''' <param name="YArray"></param>
    ''' <param name="ZTable"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>*/
double Table2DLinearInterp(double x, double y
     , double XArray[], double YArray[]
     , double **ZTable) {
     
       int Nx, Ny;
       Nx = sizeof (XArray) / sizeof *(XArray);
       Ny = sizeof (YArray) / sizeof *(YArray);
     	
       int i;
       int Xlow, Xhigh, Ylow, Yhigh;
       
        //Find Nearest X index
        Xhigh = Nx; //Prevent extrap Upper is limit
        Xlow = Xhigh;
        
        for (i=0; i<=Nx; i++) {
            if (XArray[i] > x) {
                Xhigh = i;
                if (i!=0) {
                    Xlow = i - 1;
                } else {
                    Xlow = Xhigh; //Prevent extrap Lower is limit
                }
                break;
            }
        }
        
        //Find Nearest Y index
        Yhigh = Ny;  //Prevent extrap Upper is limit
        Ylow = Yhigh;
        for (i=0; i<=Ny; i++) {
            if (YArray[i] > y) {
                Yhigh = i;
                if (i!=1) {
                    Ylow = i - 1;
            	} else {
                    Ylow = Yhigh; //Prevent extrap Lower is limit
            	}
                break;
            }
        }

 		double InterpVal1, InterpVal2, interpVal;
        if (Xlow != Xhigh) {
            InterpVal1 = Interpolate(XArray[Xlow], XArray[Xhigh], ZTable[Xlow][Ylow], ZTable[Xhigh][Ylow], x);
            if (Ylow != Yhigh) {
                InterpVal2 = Interpolate(XArray[Xlow], XArray[Xhigh], ZTable[Xlow][Yhigh], ZTable[Xhigh][Yhigh], x);
        	} else {
                interpVal = InterpVal1;
        	}
            interpVal = Interpolate(YArray[Ylow], YArray[Yhigh], InterpVal1, InterpVal2, y);
     	} else {
            interpVal = Interpolate(YArray[Ylow], YArray[Yhigh], ZTable[Xlow][Ylow], ZTable[Xlow][Yhigh], y);
     	}
     	return interpVal;
}

/*    ''' <summary>
    ''' Routine to interpolate from table using a cubic spline interpolation
    ''' </summary>
    ''' <param name="Xin">desired x value for interpolation</param>
    ''' <param name="Yin">desired y value for interpolation</param>
    ''' <param name="Xarry">array of x values from table</param>
    ''' <param name="Yarry">array of y values from table</param>
    ''' <param name="Ztable">array of z values from table;  Ztable = f(Xarry, Yarry)</param>
    ''' <returns></returns>
    ''' <remarks>
    ''' Does not extrapolate values, will set x and y to upper or lower limits
    ''' </remarks>*/
//double TableCubicInterp(double Xin, double Yin
//	, double Xarry[], double Yarry[]
//	, double **ZTable) {
//
//       int Nx, Ny;
//       Nx = sizeof (Xarry) / sizeof *(Xarry);
//       Ny = sizeof (Yarry) / sizeof *(Yarry);
//
//        //Check limits -- Do not extrapolate
//        double X;
//        if (Xin < Xarry[0]) {
//        	X = Xarry[0];
//        } else if (Xin > Xarry[Nx]) {
//        	X = Xarry[Nx];
//        } else {
//        	X = Xin;
//        }
//
//		double Y;
//        if (Yin < Yarry[0]) {
//        	Y = Yarry[0];
//        } else if (Yin > Yarry[Ny]) {
//        	Y = Yarry[Ny];
//       	} else {
//			Y = Yin;
//		}
//
//        //Interpolate Ztable=f(Xarry) at the desired xin for each Yarry
//        //Then assemble each point into a curve of Ztable=f(Yarry)
//        double X1aray[Ny], X2aray[Nx];
//        int i, j;
//        for (j = 1; j < Nx; j++) {
//            for (i = 1; i < Ny; i++) {
//                X1aray[i] = ZTable[i][j];  //create 1-D array of columnwise values to carry into cubic interpolator
//            }
//            X2aray[j] = CubicInterp(Ny, Yarry, X1aray, Yin);
//        }
//
//        //Finally, interpolate Ztable=f(Yarry) at the desired yin
//        double yinterp = CubicInterp(Nx, Xarry, X2aray, Xin);
//		return yinterp;
//}

/*    ''' <summary>
    ''' This routine performs a double interpolation given z and x, and determines the corresponding 
    '''  for Yarray
    '''  example of when this routine would be useful is the following:
    '''  A matrix of data has been compiled for enthalpy, given temperature and pressure.
    '''  If the user knows the H and P, the would then want to know what the T is
    '''  without having to create a matrix of T, given H and P.
    '''  In this example, z represents the enthalpy value, x is the pressure, Xarray is the pressure array, 
    '''  Yarray is the temperature array, and Zvals is the H values given T and P. 
    ''' </summary>
    ''' <param name="z"></param>
    ''' <param name="x"></param>
    ''' <param name="Xarray"></param>
    ''' <param name="Yarray"></param>
    ''' <param name="Zvals"></param>
    ''' <param name="OrderYX"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
 double DoubleInterpAlternate(ByVal z As Double, ByVal x As Double, ByRef Xarray() As Double, _
      ByRef Yarray() As Double, ByRef Zvals(,) As Double, Optional ByVal OrderYX As Boolean = False) As Double

        Dim i, j, k, iLow, iHigh, kLow, kHigh As Integer
        Dim xLow, xHigh, zLow, zHigh As Double
        Dim Xvals(), Yvals() As Double

        //Account for the possibility that the z matrix is formatted with y, then x
        If OrderYX Then
            ReDim Xvals(UBound(Yarray))
            ReDim Yvals(UBound(Xarray))
            For i = 1 To UBound(Xarray)
                Yvals(i) = Xarray(i)
            Next i
            For j = 1 To UBound(Yarray)
                Xvals(j) = Yarray(j)
            Next j
        Else
            ReDim Xvals(UBound(Xarray))
            ReDim Yvals(UBound(Yarray))
            For i = 1 To UBound(Xarray)
                Xvals(i) = Xarray(i)
            Next i
            For j = 1 To UBound(Yarray)
                Yvals(j) = Yarray(j)
            Next j
        End If

        //Find closest points within x array
        //Initializing
        If UBound(Xarray) = 1 Then
            iLow = 1
            iHigh = 1
            xLow = Xvals(1)
            xHigh = Xvals(1)
        Else
            xLow = -1000000000000.0
            xHigh = 1000000000000.0
            For i = 1 To UBound(Xarray)
                If Xvals(i) > xLow And Xvals(i) <= x Then
                    xLow = Xvals(i)
                    iLow = i
                End If
                If Xvals(i) < xHigh And Xvals(i) > x Then
                    xHigh = Xvals(i)
                    iHigh = i
                End If
            Next i
        End If

        //Now let's create a new z array at the specified x location along the y vector
        Dim Zarray() As Double
        ReDim Zarray(UBound(Yvals))
        For j = 1 To UBound(Zarray)
            If OrderYX Then
                If iLow = iHigh Then
                    Zarray(j) = Zvals(iLow, j)
                Else
                    Zarray(j) = Interpolate(Xvals(iLow), x, Xvals(iHigh), Zvals(iLow, j), Zvals(iHigh, j))
                End If
            Else
                If iLow = iHigh Then
                    Zarray(j) = Zvals(j, iLow)
                Else
                    Zarray(j) = Interpolate(Xvals(iLow), x, Xvals(iHigh), Zvals(j, iLow), Zvals(j, iHigh))
                End If
            End If
        Next j

        //Finally, find closest points along y vector and interpolate
        //Initializing
        zLow = -1000000000000.0
        zHigh = 1000000000000.0
        For k = 1 To UBound(Zarray)
            If Zarray(k) > zLow And Zarray(k) <= z Then
                zLow = Zarray(k)
                kLow = k
            End If
            If Zarray(k) < zHigh And Zarray(k) > z Then
                zHigh = Zarray(k)
                kHigh = k
            End If
        Next k

        Return Interpolate(Zarray(kLow), z, Zarray(kHigh), Yvals(kLow), Yvals(kHigh))
}*/


/*    ''' <summary>
    ''' Reverse Double linear interpolation function
    ''' Linearly interpolates on 2-D array
    ''' Does not extrapolate (if x or y outside range then Eval at Range Limit
    ''' Arrays must be sorted in increasing sequential order
    ''' Table is 2-D array (x,y)
    ''' </summary>
    ''' <param name="Zval"></param>
    ''' <param name="y"></param>
    ''' <param name="XArray"></param>
    ''' <param name="YArray"></param>
    ''' <param name="ZTable"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
 double Table2DReverseLinearInterp(ByVal Zval As Double, ByVal y As Double, _
      ByRef XArray() As Double, ByRef YArray() As Double, ByRef ZTable(,) As Double) As Double

        Dim i As Integer
        Dim InterpVal1, InterpVal2 As Double
        Dim Ylow, Yhigh As Integer
        Dim FoundYlowZ, FoundYhighZ As Boolean
        'Find Nearest Y index
        Yhigh = UBound(YArray)  'Prevent extrap Upper is limit
        Ylow = Yhigh
        For i = 1 To UBound(YArray)
            If YArray(i) > y Then
                Yhigh = i
                If i <> 1 Then
                    Ylow = i - 1
                Else
                    Ylow = Yhigh 'Prevent extrap Lower is limit
                    FoundYhighZ = True
                End If
                Exit For
            End If
        Next i


        'Find the nearest Z index and interpolate
        For i = 1 To UBound(XArray)
            'Ylow
            If Not FoundYlowZ And ZTable(i, Ylow) > Zval Then
                FoundYlowZ = True
                If i <> 1 Then
                    InterpVal1 = Interpolate(ZTable(i - 1, Ylow), Zval, ZTable(i, Ylow), XArray(i - 1), XArray(i))
                Else
                    InterpVal1 = XArray(1) 'ZTable(1, Ylow) 'Prevent extrap Lower is limit
                End If
            End If
            'Yhigh
            If Not FoundYhighZ And ZTable(i, Yhigh) > Zval Then
                FoundYhighZ = True
                If i <> 1 Then
                    InterpVal2 = Interpolate(ZTable(i - 1, Yhigh), Zval, ZTable(i, Yhigh), XArray(i - 1), XArray(i))
                Else
                    InterpVal2 = XArray(1) 'ZTable(1, Yhigh) 'Prevent extrap Lower is limit
                End If
            End If
            If FoundYlowZ And FoundYhighZ Then Exit For
        Next i
        'Prevent extrap Upper is limit
        If Not FoundYlowZ Then InterpVal1 = ZTable(UBound(XArray), Ylow)
        If Not FoundYhighZ Then InterpVal2 = ZTable(UBound(XArray), Yhigh)


        'INterpolate in Y direction
        If Ylow <> Yhigh Then
            Return Interpolate(YArray(Ylow), y, YArray(Yhigh), InterpVal1, InterpVal2)
        Else
            Return InterpVal1
        End If
 }*/


/*    ''' <summary>
    ''' This routine interpolates within the input arrays and find the intercept of the two 
    ''' Extrapolation is not done 
    ''' Arrays are 1 based and must be sorted in sequential increasing order
    ''' Developed by David Giel 1-16-08
    ''' Note that routine has been validated for the cases which it has been developed
    ''' given the numerious potential applications of this routine 
    ''' it most likely has some robustness bugs.
    ''' If you plan to use this routine please test it for your specific application
    ''' </summary>
    ''' <param name="Xarry1"></param>
    ''' <param name="Yarry1"></param>
    ''' <param name="Xarry2"></param>
    ''' <param name="Yarry2"></param>
    ''' <param name="xInt"></param>
    ''' <param name="Exists"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
double FindInterSectOFTwoDataSets(ByRef Xarry1() As Double, ByRef Yarry1() As Double, _
      ByRef Xarry2() As Double, ByRef Yarry2() As Double, ByRef xInt As Double, _
      ByRef Exists As Boolean) As Double

        Dim i, j, iMin1, iMax1, jMin2, jMax2 As Integer //markers for global range index
        Dim x1L, x2L, x1H, x2H As Integer //markers for cross over range
        Dim dY, dYlast As Double //indicators of y intercept
        //marks the 4 points which seround the intercept
        Dim xi1, xi2, Interp1Low, Interp1High, Interp2Low, Interp2High As Integer

        iMin1 = 1
        iMax1 = UBound(Xarry1)
        jMin2 = 1
        jMax2 = UBound(Xarry2)

        Exists = False
        //x arrays must overlap
        //start at min of Xarray of 1
        //march up the Xarray of 1 until the min bound of 2 is crossed 
        //this is the start of overlap
        For i = iMin1 To iMax1
            If Xarry1(i) >= Xarry2(jMin2) Then
                //crossed or reached X2 min
                If Xarry1(i) = Xarry2(jMin2) Then
                    x1L = i
                Else //cross in between use prev
                    x1L = Max(i - 1, iMin1)
                End If

                If i = iMin1 Then
                    //find the lowest of 2
                    //step down until xl1 is crossed
                    j = jMax2
                    Do
                        x2L = j
                        If Xarry1(i) >= Xarry2(j) Then Exit Do
                        If j = jMin2 Then Exists = False
                        j -= 1
                    Loop Until j = jMin2 - 1
                Else
                    x2L = jMin2
                End If
                Exists = True
                GoTo upper
            End If
        Next i

upper:  If Exists Then
            Exists = False
            //step down from X1 max until X2 max is crossed
            i = iMax1
            Do
                j = jMax2
                Do
                    If Xarry1(i) <= Xarry2(j) Then
                        //'crossed or reached the upper bound of array 2
                        x1H = i
                        If Xarry1(i) = Xarry2(j) Then
                            x2H = j //2 does not extend past 1
                        Else
                            //cross inbetween use next j if not max
                            x2H = Min(j + 1, jMax2)
                        End If
                        Exists = True
                        GoTo Ycross
                    End If
                Loop
            Loop

Ycross:     //now using xL and xH the non-overlapping parts of arrays are chopped off
            If Exists Then
                //starting at the lowest x val march in increasing x until y cross
                xi1 = x1L
                xi2 = x2L
                Do
                    //check if x vals are equal, 1 is less, or 2 is less
                    If Xarry1(xi1) = Xarry2(xi2) Then
                        //can check dY driectly
                        dYlast = dY
                        dY = Yarry1(xi1) - Yarry2(xi2)

                        //check for a change in sign if not first time
                        If dYlast <> 0 And dY / Abs(dY) <> dYlast / Abs(dYlast) Then
                            //sign change just occured so the upper interp val is known
                            Interp1High = xi1
                            Interp2High = xi2
                            Exit Do
                        Else
                            //no sign chnage yet so update the lower interp val
                            Interp1Low = xi1
                            Interp2Low = xi2
                        End If




                    ElseIf Xarry1(xi1) < Xarry2(xi2) Then
                        //check y2 directly
                        //need to interpolate for y1 at x2
                        dYlast = dY
                        Dim NextX1, NextY1 As Double
                        If xi1 + 1 > x1H Then
                            //no crossing exists
                            Exists = False
                            Exit Do
                        Else
                            NextX1 = Xarry1(xi1 + 1)
                            NextY1 = Yarry1(xi1 + 1)
                        End If
                        dY = Interpolate(Xarry1(xi1), Xarry2(xi2), NextX1, Yarry1(xi1), NextY1) - Yarry2(xi2)


                        //check for a change in sign if not first time
                        If dYlast <> 0 And dY / Abs(dY) <> dYlast / Abs(dYlast) Then
                            //sign change just occured so the upper interp val is known
                            Interp1High = xi1 + 1
                            Interp2High = xi2
                            Exit Do
                        Else
                            //no sign chnage yet so update the lower interp val
                            Interp1Low = xi1
                            Interp2Low = xi2
                        End If



                    ElseIf Xarry1(xi1) > Xarry2(xi2) Then
                        //check y1 directly
                        //need to interpolate for y2 at x1
                        dYlast = dY
                        Dim NextX2, NextY2 As Double
                        If xi2 + 1 > x2H Then
                            //no crossing exists
                            Exists = False
                            Exit Do
                        Else
                            NextX2 = Xarry2(xi2 + 1)
                            NextY2 = Yarry2(xi2 + 1)
                        End If
                        dY = Yarry1(xi1) - Interpolate(Xarry2(xi2), Xarry1(xi1), NextX2, Yarry2(xi2), NextY2)


                        //check for a change in sign if not first time
                        If dYlast <> 0 And dY / Abs(dY) <> dYlast / Abs(dYlast) Then
                            //sign change just occured so the upper interp val is known
                            Interp1High = xi1
                            Interp2High = xi2 + 1
                            Exit Do
                        Else
                            //no sign chnage yet so update the lower interp val
                            Interp1Low = xi1
                            Interp2Low = xi2
                        End If

                    End If

                    //increase index (can increase both if next are equal)
                    If Xarry1(Min(xi1 + 1, x1H)) = Xarry2(Min(xi2 + 1, x2H)) Then
                        xi1 = Min(xi1 + 1, x1H)
                        xi2 = Min(xi2 + 1, x2H)
                    ElseIf Xarry1(Min(xi1 + 1, x1H)) < Xarry2(Min(xi2 + 1, x2H)) Then
                        //1 is less increase only 1
                        xi1 = Min(xi1 + 1, x1H)
                    ElseIf Xarry2(Min(xi2 + 1, x2H)) < Xarry1(Min(xi1 + 1, x1H)) Then
                        //2 is less increase only 2
                        xi2 = Min(xi2 + 1, x2H)
                    End If

                    If xi1 = x1H And xi2 = x2H Then
                        Exists = False
                        Exit Do
                    End If

                Loop
            End If
        End If

        If Exists Then
            //use interp values to calculate the intersection point
            Dim m1, m2, b1, b2 As Double
            m1 = (Yarry1(Interp1High) - Yarry1(Interp1Low)) / (Xarry1(Interp1High) - Xarry1(Interp1Low))
            m2 = (Yarry2(Interp2High) - Yarry2(Interp2Low)) / (Xarry2(Interp2High) - Xarry2(Interp2Low))
            b1 = Yarry1(Interp1Low) - m1 * Xarry1(Interp1Low)
            b2 = Yarry2(Interp2Low) - m2 * Xarry2(Interp2Low)
            xInt = (b2 - b1) / (m1 - m2)

            Return Interpolate(Xarry1(Interp1Low), xInt, Xarry1(Interp1High), Yarry1(Interp1Low), Yarry1(Interp1High))
        Else
            Return Nothing
        End If
}*/
