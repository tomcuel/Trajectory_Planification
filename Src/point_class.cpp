#include "point_class.hpp"



//==========================================================================
// class Point : 2D point
//==========================================================================

// member functions
//==========================================================================
// (0,0) constructor
Point::Point()
{
    x = 0;
    y = 0;
}
// (a,a) constructor 
Point::Point(double a) : x(a), y(a)
{

}
// (a,b) constructor 
Point::Point(double a, double b) : x(a), y(b)
{

}
// P+=Q 
Point& Point::operator +=(const Point& Q)
{
    this->x += Q.x;
    this->y += Q.y;
    return *this;
}
// P-=Q
Point& Point::operator -=(const Point& Q)
{
    this->x -= Q.x;
    this->y -= Q.y;
    return *this;
}
// P+=a (add a to each component)
Point& Point::operator +=(const double& a)
{
    this->x += a;
    this->y += a;
    return *this;
}
// P-=a (substract a to each component)
Point& Point::operator -=(const double& a)
{
    this->x -= a;
    this->y -= a;
    return *this;
}
// P*=a
Point& Point::operator *=(const double& a)
{   
    this->x *= a;
    this->y *= a;
    return *this;
}
// P/=a
Point& Point::operator /=(const double& a)
{
    if (a == 0){
        throw invalid_argument("Error: Division by zero.");
    }
    this->x /= a;
    this->y /= a;
    return *this;
}

// external functions
//==========================================================================
// <<
ostream& operator <<(ostream& out, const Point& P)
{   
    out<<"("<<P.x<<","<<P.y<<")";
    return out;
}
// +P
Point operator +(const Point& P)
{
    Point result = P;
    return result;
}
// -P
Point operator -(const Point& P)
{
    Point result = P;
    result *= -1;
    return result;
}
// P+Q 
Point operator +(const Point & P, const Point & Q)
{
    Point result = P;
    result += Q;
    return result;
}
// P-Q
Point operator -(const Point & P, const Point & Q)
{
    Point result = P;
    result -= Q;
    return result;
}
// P+a
Point operator +(const Point & P, double a)
{
    Point result = P;
    result += a;
    return result;
}
// a+P
Point operator +(double a, const Point & P)
{
    return P + a;
}
// P-a
Point operator -(const Point & P, double a)
{
    Point result = P;
    result -= a;
    return result;
}
// a-P
Point operator -(double a, const Point & P)
{
    return P - a;
}
// P*a
Point operator *(const Point & P, double a)
{
    Point result = P;
    result *= a;
    return result;
}
// a*P
Point operator *(double a, const Point & P)
{
    return P * a;
}
// P/a
Point operator /(const Point & P, double a)
{
    Point result = P;
    result /= a;
    return result;
}
// P==Q
bool operator ==(const Point & P, const Point & Q)
{
    return P.x == Q.x && P.y == Q.y;
}
// P!=Q
bool operator !=(const Point & P, const Point & Q)
{
    return !(P == Q);
}
// P<Q
bool operator <(const Point & P, const Point & Q)
{
    return P.x < Q.x || (P.x == Q.x && P.y < Q.y);
}
// P|Q
double operator |(const Point & P, const Point & Q)
{
    return P.x * Q.x + P.y * Q.y;
}
// |P|
double norme(const Point & P)
{ 
    return sqrt(P|P);
}
// surface of a triangle
double surface(const Point & A, const Point & B, const Point & C){
    return fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y)) / 2.0;
}
