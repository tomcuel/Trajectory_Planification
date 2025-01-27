#include <iostream>
#include <cmath>
#include <functional>
#include <string>
using namespace std;

struct Point
{
    double x, y, z;

    friend ostream& operator<<(ostream& out, const Point& p) 
    {
        out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
        return out;
    }
};

class Param_Fun 
{
public:
    function<Point(double)> func;
    string formula;

    Param_Fun(function<Point(double)> f, const string& formula) : func(f), formula(formula) 
    {
    }
    Point operator()(double t) const 
    {
        return func(t);
    }
};

int main() {
    // Define a symbolic parametric function
    Param_Fun param(
        [](double t) { return Point{cos(t), sin(t), t}; },
        "x(t) = cos(t), y(t) = sin(t), z(t) = t"
    );

    double t = 1.0;
    Point p = param(t); // Evaluate at t = 1.0
    cout << "Point: " << p << "\n";
    cout << "Formula: " << param.formula << "\n";

    return 0;
}
