//==========================================================================
// file defining a 2D point
//==========================================================================
#ifndef __POINT_CLASS_HPP__
#define __POINT_CLASS_HPP__


#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
using namespace std;



//==========================================================================
// class Point : 2D point
//==========================================================================
class Point
{
public : 
    double x, y; // 2D point coordinates
    
    // constructors
    Point(); // (0,0)
    Point(double a); // (a,a)
    Point(double a, double b); // (a,b)

    Point& operator +=(const Point& Q); // P+=Q
    Point& operator -=(const Point& Q); // P-=Q
    Point& operator +=(const double& a); // P+=a (ajout de a à chaque composante)
    Point& operator -=(const double& a); // P-=a (soustraction de a à chaque composante)
    Point& operator *=(const double& a); // P*=a
    Point& operator /=(const double& a); // P/=a
};
//==========================================================================
// <<
ostream& operator <<(ostream& out, const Point& P);
// +P
Point operator +(const Point& P);
// -P
Point operator -(const Point& P);
// P+Q 
Point operator +(const Point & P, const Point & Q);
// P-Q
Point operator -(const Point & P, const Point & Q);
// P+a
Point operator +(const Point & P, double a);
// a+P
Point operator +(double a, const Point & P);
// P-a
Point operator -(const Point & P, double a);
// a-P
Point operator -(double a, const Point & P);
// P*a
Point operator *(const Point & P, double a);
// a*P
Point operator *(double a, const Point & P);
// P/a
Point operator /(const Point & P, double a);
// P==Q
bool operator ==(const Point & P, const Point & Q);
// P!=Q
bool operator !=(const Point & P, const Point & Q);
// P<Q
bool operator <(const Point & P, const Point & Q);
// P|Q
double operator |(const Point & P, const Point & Q);
// |P|
double norme(const Point & P);
// surface of a triangle
double surface(const Point & A, const Point & B, const Point & C);



#endif // __POINT_CLASS_HPP__