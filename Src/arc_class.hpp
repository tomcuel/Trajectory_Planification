//==========================================================================
// file defining the Arcs that links the points, and the Param_Fun that defines the parametric function and its symbolic representation
//==========================================================================
#ifndef __ARC_CLASS_HPP__
#define __ARC_CLASS_HPP__


#include "point_class.hpp"
#include <string>
#include <functional>
#include <fstream>



//==========================================================================
// class Param_Fun defining what a parametric function is, along with its symbolic representation
//==========================================================================
class Param_Fun 
{
public:
    function<Point(double)> func;
    string formula;

    Param_Fun(function<Point(double)> f, const string& formula);
    Point operator()(double t) const;
};



//==========================================================================
// class Arc 
//==========================================================================
class Arc
{
public :
    string type=""; // arc type
    vector<Point> points; // arc extremities
    virtual double length() const = 0; // arc length
    virtual Arc& reverse() = 0; // reverse the arc (inversion of extremities)
    virtual void exporte(ostream& out) const = 0; // export the arc in a file
};
//==========================================================================
ostream& operator <<(ostream& out, const Arc& A); // << operator for Arc



//==========================================================================
// class Segment heriting from Arc
//==========================================================================
class Segment : public Arc 
{
public : 
    Segment(const Point& A, const Point& B); // constructor
    double length() const override; // function length 
    Arc& reverse() override; // reverse the arc (inversion of extremities of the segment)
    void exporte(ostream &out) const override; // export the arc in a file
};



//==========================================================================
// class ArcCircle heriting from Arc 
//==========================================================================
class ArcCircle : public Arc 
{
public :
    Point center; // center of the circle                                    
    double radius, angle1, angle2; // radius and angles of the arc extremities
    ArcCircle(const Point& C, double R, double a1, double a2); // constructor
    double length() const override; // function length
    Arc& reverse() override; // reverse the arc (inversion of extremities and the corresponding angles)
    void exporte(ostream &out) const override; // export the arc in a file
};



//==========================================================================
// class ArcParameter heriting from Arc
//==========================================================================
class ArcParameter : public Arc 
{
public : 
    Param_Fun parametrisation; // parametrisation function
    double t0, tf; // interval of the parametrisation
    int nbt; // number of intervals to calculate the approximate value
    ArcParameter(Param_Fun f, double t_debut, double t_fin, int n); // constructor
    double length() const override; // function length
    Arc& reverse() override; // reverse the arc (inversion of extremities and each value of the parametrisation)
    void exporte(ostream &out) const override; // export the arc in a file
};



//==========================================================================
// functions defining whether two arcs intersect
//==========================================================================
// check if one point is on a segment
bool is_point_on_segment(const Point& P, const Segment& seg);
// get the intersection point of two segments
Point get_segment_intersection(const Segment& seg1, const Segment& seg2);
// check if two segments intersect
bool do_segments_intersect(const Segment& seg1, const Segment& seg2);
// check the relative orientation of three points
int orientation(const Point& p, const Point& q, const Point& r);
// check if seg2 is contained within seg1
bool is_segment_contained_within(const Segment& seg1, const Segment& seg2);
// check if two segments are collinear
bool are_segments_collinear(const Segment& seg1, const Segment& seg2);
// function to check if two segments are overlapping
bool are_overlapping(const Segment& seg1, const Segment& seg2);
// check if a segment intersects with an arc circle
bool do_segment_intersect_arc_circle(const Segment& segment, const ArcCircle& arc_circle);
// check if a segment intersects with a parametric arc
bool do_segment_intersect_parametric(const Segment& segment, const ArcParameter& param_arc);
// check if two arc circles intersect
bool do_arc_circle_intersect(const ArcCircle& arc1, const ArcCircle& arc2);
// check if an arc circle intersects with a parametric arc
bool do_arc_circle_intersect_parametric(const ArcCircle& arc_circle, const ArcParameter& param_arc);
// check if two parametric arcs intersect
bool do_parametric_intersect(const ArcParameter& arc1, const ArcParameter& arc2);



#endif // __ARC_CLASS_HPP__