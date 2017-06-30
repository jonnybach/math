
#include <math.h>
#include "Geometry.h"

double HydDiaRect(const double L, const double W) {
        //Calculates hydraulic diameter of a rectangle of length L and width W
        return (2 * L * W / (L + W));
}

double CalcDist2D(const double x1, const double y1, const double x2, const double y2) {
	return sqrt(pow((x2 - x1),2) + pow((y2 - y1),2));
}

double CalcDist3D(const double x1, const double y1, const double z1
	, const double x2, const double y2, const double z2) {
	return sqrt(pow((x2 - x1),2) + pow((y2 - y1),2) + pow((z2 - z1),2));
}

double CalcSlope(const double x1, const double y1, const double x2,const double y2) {
    double m;
    if (x2 == x1) { //avoid divide by zero for vertical lines
        m = 99999999999999999999999999999999999999.999;
    } else {
        m = (y2 - y1) / (x2 - x1);
    }
    return m;
}


double AreaTriangle(const double LengthA, const double LengthB, const double LengthC) {
        double s = (LengthA + LengthB + LengthC) / 2;
        double Area = sqrt(s * (s - LengthA) * (LengthB) * (LengthC));
        return Area;
}

/*
void GetPtAlongLine(ByVal x1 As Double, ByVal y1 As Double, ByVal z1 As Double, _
      ByVal x2 As Double, ByVal y2 As Double, ByVal z2 As Double, ByVal D As Double, ByRef x3 As Double, _
      ByRef y3 As Double, ByRef z3 As Double)
        'Calculates coordinates for a point (x3,y3,z3) that is a distance D from point (x1,y1,z1) along the line
        'connecting points (x1,y1,z1) and (x2,y2,z2).
        'Point (x3,y3,z3) is between (x1,y1,z1) and (x2,y2,z2)

        Dim Dist_1to2 As Double

        Dist_1to2 = CalcDist3D(x1, y1, z1, x2, y2, z2)

        x3 = D / Dist_1to2 * (x2 - x1)
        y3 = D / Dist_1to2 * (y2 - y1)
        z3 = D / Dist_1to2 * (z2 - z1)
}

void FindClosestPoint(ByVal RefPt As DataPoint, ByRef OppPts() As DataPoint, _
      ByRef ClosestPoint As DataPoint)

        Dim i As Integer
        Dim Dist, MinDist As Double

        MinDist = 1000000000000.0
        For i = LBound(OppPts) To UBound(OppPts)
            Dist = CalcDist(RefPt.X, RefPt.Y, OppPts(i).X, OppPts(i).Y)
            If Dist < MinDist Then
                MinDist = Dist
                ClosestPoint = OppPts(i)
            End If
        Next i
}
*/

/*double CalcDistPnt(DataPoint Pt1, DataPoint Pt2) {
        Dim Dist As Double = Sqrt((Pt2.X - Pt1.X) ^ 2 + (Pt2.Y - Pt1.Y) ^ 2)
        Return Dist
}*/

/*
Public Shared Sub CentroidQuad(ByVal Pt1 As DataPoint, ByVal Pt2 As DataPoint, ByVal Pt3 As DataPoint, _
      ByVal Pt4 As DataPoint, ByRef PtC As DataPoint)
        'Calculates centroid of a quadrilateral, given the four vertices
        'This routine assumes that the vertices have ben arranged in a clockwise or counter-clockwise manner;
        'ie. there are no intersecting lines

        Dim CentTriangle(4) As DataPoint
        For i As Integer = 1 To 4
            CentTriangle(i) = New DataPoint
        Next i

        Dim m1, m2, b1, b2 As Double

        'Splitting the quadrilateral into four triangles
        '1) (x1,y1), (x2,y2), and (x3,y3)
        '2) (x1,y1), (x3,y3), and (x4,y4)
        '3) (x1,y1), (x2,y2), and (x4,y4)
        '4) (x2,y2), (x3,y3), and (x4,y4)

        'Now get centroids of the four triangles
        CentroidTriangle(Pt1, Pt2, Pt4, CentTriangle(1))
        CentroidTriangle(Pt1, Pt4, Pt3, CentTriangle(2))
        CentroidTriangle(Pt1, Pt2, Pt3, CentTriangle(3))
        CentroidTriangle(Pt2, Pt4, Pt3, CentTriangle(4))

        'The centroid of the quadrilateral is the intersection point between the lines connecting the centroids of the
        'triangles.
        'Calculating slopes and intercepts of the lines connecting the triangle centroids
        m1 = CalcSlope(CentTriangle(1), CentTriangle(2))
        b1 = CentTriangle(1).Y - m1 * CentTriangle(1).X
        m2 = CalcSlope(CentTriangle(3), CentTriangle(4))
        b2 = CentTriangle(3).Y - m2 * CentTriangle(3).X

        'Finally, get centroid coordinates by the intersection of these two lines
        CalcVertex(m1, b1, m2, b2, PtC)

    End Sub

    Public Shared Sub CentroidTriangle(ByVal Pt1 As DataPoint, ByVal Pt2 As DataPoint, ByVal Pt3 As DataPoint, _
      ByRef PtC As DataPoint)

        PtC.X = (Pt1.X + Pt2.X + Pt3.X) / 3
        PtC.Y = (Pt1.Y + Pt2.Y + Pt3.Y) / 3

    End Sub

    Public Shared Sub CheckQuadVertices(ByRef x1 As Double, ByRef y1 As Double, ByRef x2 As Double, _
      ByRef y2 As Double, ByRef x3 As Double, ByRef y3 As Double, ByRef x4 As Double, ByRef y4 As Double)

        'This subroutine checks the four vertices of a quadrilateral and re-orders them to avoid crossed lines

        Dim m(3) As Double, mMax As Double, mMin As Double, i As Integer, iMax As Integer, iMin As Integer
        Dim x(4) As Double, y(4) As Double, xNew(4) As Double, yNew(4) As Double

        x(1) = x1
        y(1) = y1
        x(2) = x2
        y(2) = y2
        x(3) = x3
        y(3) = y3
        x(4) = x4
        y(4) = y4

        m(1) = CalcSlope(x(1), y(1), x(2), y(2))
        m(2) = CalcSlope(x(1), y(1), x(3), y(3))
        m(3) = CalcSlope(x(1), y(1), x(4), y(4))

        mMax = -1000000000000.0#
        mMin = 1000000000000.0#

        For i = 1 To 3
            If m(i) >= mMax Then
                mMax = m(i)
                iMax = i
            End If
            If m(i) <= mMin Then
                mMin = m(i)
                iMin = i
            End If
        Next i

        'Initialize new x and y matrices
        For i = 1 To 4
            xNew(i) = x(i)
            yNew(i) = y(i)
        Next i

        If Abs(iMax - iMin) = 1 Then 'Re-order vertices
            For i = 1 To 3
                If i = iMax Then
                    xNew(4) = x(i + 1)
                    yNew(4) = y(i + 1)
                ElseIf i = iMin Then
                    xNew(2) = x(i + 1)
                    yNew(2) = y(i + 1)
                Else
                    xNew(3) = x(i + 1)
                    yNew(3) = y(i + 1)
                End If
            Next i
        End If

        x1 = xNew(1)
        y1 = yNew(1)
        x2 = xNew(2)
        y2 = yNew(2)
        x3 = xNew(3)
        y3 = yNew(3)
        x4 = xNew(4)
        y4 = yNew(4)

    End Sub

    Public Overloads Shared Sub CalcVertex(ByVal m1 As Double, ByVal b1 As Double, ByVal m2 As Double, _
      ByVal b2 As Double, ByRef x As Double, ByRef y As Double)

        If m1 <> m2 Then
            x = (b2 - b1) / (m1 - m2)
            y = m1 * x + b1
        Else
            x = Double.MaxValue
            y = Double.MaxValue
        End If

    End Sub

    Public Overloads Shared Sub CalcVertex(ByVal m1 As Double, ByVal b1 As Double, ByVal m2 As Double, _
      ByVal b2 As Double, ByRef Pt As DataPoint)

        If m1 <> m2 Then
            Pt.X = (b2 - b1) / (m1 - m2)
            Pt.Y = m1 * Pt.X + b1
        Else
            Pt.X = Double.MaxValue
            Pt.Y = Double.MaxValue
        End If

    End Sub
*/

/*
Public Overloads Shared Function CalcSlope(ByVal Pt1 As DataPoint, ByVal Pt2 As DataPoint) As Double

        Dim m As Double

        If Pt1.X = Pt2.X Then 'avoid divide by zero for vertical lines
            m = Double.MaxValue
        Else
            m = (Pt2.Y - Pt1.Y) / (Pt2.X - Pt1.X)
        End If

        Return m

End Function

Public Shared Function CalcIntercept(ByVal Slope As Double, ByVal Pt As DataPoint) As Double

        Dim b As Double = Pt.Y - Slope * Pt.X

        Return b
End Function
*/
