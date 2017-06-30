#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

//void spline(const int numElems, const double x[], const double y[], const double yp1, const double ypn, double y2[]);
//
//void splint(const int numElems, const double xa[], const double ya[], const double y2a[], const double x, double *y);
//
//void Derivative(const int numElems, const double x[], const double y[], double *yp1, double *ypn);
//
//double CubicInterp(const int numElems, const double Xa[], const double Ya[], const double Xin);

double LinearInterpUnsorted(const int numElems, const double Xaray[], const double Yaray[]
	, const double x, const int extrap);

double Interpolate(const double xlow, const double xup, const double ylow
	, const double yup, const double xinterp);

double DoubleLinearInterp(double Var1, double Var2
	, double Var1Array[], double Var2Array[]
	, double **TableArray
	, int extrap);
	
double Table2DLinearInterp(double x, double y
     , double XArray[], double YArray[]
     , double **ZTable);
	
double TableCubicInterp(double Xin, double Yin
	, double Xarry[], double Yarry[]
	, double **ZTable);

#endif 
/*INTERPOLATION_H_*/
