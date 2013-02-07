#ifndef GEOMETRY_H_
#define GEOMETRY_H_

double HydDiaRect(const double L, const double W);

double CalcDist2D(const double x1, const double y1, const double x2, const double y2);
double CalcDist3D(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);

double CalcSlope(const double x1, const double y1, const double x2,const double y2);

double AreaTriangle(const double LengthA, const double LengthB, const double LengthC);

#endif /*GEOMETRY_H_*/
