#include "arc_class.hpp"



//==========================================================================
// class Param_Fun defining what a parametric function is 
//==========================================================================

// member functions
//==========================================================================
// constructor
Param_Fun::Param_Fun(function<Point(double)> f, const string& formula) : func(f), formula(formula) 
{

}
// operator() to evaluate the function at a given value
Point Param_Fun::operator()(double t) const 
{
    return func(t);
}



//==========================================================================
// class Arc
//==========================================================================

// external functions
//==========================================================================
//<<operator for Arc
ostream& operator <<(ostream& out, const Arc& A)
{
    out<<A.type<<" : "<<A.points[0]<<" -> "<<A.points[1];
    return out;
}



//==========================================================================
// class Segment
//==========================================================================

// member functions
//==========================================================================
// constructor
Segment::Segment(const Point& A, const Point& B)
{
    type = "Segment";
    points.push_back(A);
    points.push_back(B);
}
// length of the segment 
double Segment::length() const 
{
    return norme(points[1] - points[0]);
}
// reverse the arc (inversion of extremities of the segment)
Arc& Segment::reverse() 
{
    swap(points[0], points[1]);
    return *this;
}
// export the arc in a file
void Segment::exporte(ostream &out) const 
{
    out<<"Segment : "<<points[0]<<" -> "<<points[1]<<endl;
}



//==========================================================================
// class ArcCircle
//==========================================================================

// member functions
//==========================================================================
// constructor 
ArcCircle::ArcCircle(const Point& C, double R, double a1, double a2) : center(C), radius(R), angle1(a1), angle2(a2) 
{
    type = "ArcCircle";
    points.push_back(C + Point(radius * cos(angle1), radius * sin(angle1)));
    points.push_back(C + Point(radius * cos(angle2), radius * sin(angle2)));
}
// length of the arc circle
double ArcCircle::length() const 
{
    return abs(angle1-angle2)*radius;
}
// reverse the arc (inversion of extremities and the corresponding angles)
Arc& ArcCircle::reverse() 
{
    swap(points[0], points[1]);
    swap(angle1, angle2);
    return *this;
}
// export the arc in a file
void ArcCircle::exporte(ostream &out) const 
{
    out<<"ArcCircle : "<<points[0]<<" -> "<<points[1]<<" (center : "<<center<<", radius : "<<radius<<", angles : "<<angle1<<" -> "<<angle2<<")"<<endl;
}



//==========================================================================
//class ArcParameter
//==========================================================================

// member functions
//==========================================================================
// constructor 
ArcParameter::ArcParameter(Param_Fun f, double t_debut, double t_fin, int n) : parametrisation(f), t0(t_debut), tf(t_fin), nbt(n) 
{
        type = "ArcParameter";
        points.push_back(f(t0));
        points.push_back(f(tf));
}
// length of the arc parameter
double ArcParameter::length() const
{
    double result = 0;
    int n = nbt -1;
    double dt = (tf-t0) / n;
    for (int k=0; k<n; k++){
        result += norme(parametrisation(t0 + (k+1) * dt) - parametrisation(t0 + k * dt));
    }
    return result;
}                             
// reverse the arc (inversion of extremities and each value of the parametrisation)
Arc& ArcParameter::reverse()
{
    swap(points[0], points[1]);
    swap(t0, tf);
    return *this;
}                                   
// export the arc in a file
void ArcParameter::exporte(ostream &out) const
{
    out<<"ArcParameter : "<<points[0]<<" -> "<<points[1]<<" (parametrisation from "<<t0<<" to "<<tf<<" in "<<nbt<<" intervals with the formula : "<<parametrisation.formula<<")"<<endl;
} 



//==========================================================================
// functions defining whether two arcs intersect
//==========================================================================
// check if one point is on a segment
bool is_point_on_segment(const Point& P, const Segment& seg)
{
    // check if the point is on the segment
    double d1 = norme(seg.points[0] - P);
    double d2 = norme(seg.points[1] - P);
    double d3 = norme(seg.points[0] - seg.points[1]);
    return abs(d1 + d2 - d3) < 1e-4;
}
// get the intersection point of two segments
Point get_segment_intersection(const Segment& seg1, const Segment& seg2)
{
    double x1 = seg1.points[0].x, y1 = seg1.points[0].y;
    double x2 = seg1.points[1].x, y2 = seg1.points[1].y;
    double x3 = seg2.points[0].x, y3 = seg2.points[0].y;
    double x4 = seg2.points[1].x, y4 = seg2.points[1].y;

    // calculate the denominators
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (fabs(denom) < 1e-9) {
        return false; // segments are parallel or coincident
    }

    // calculate the intersection point using parameter t and u
    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
    double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

    Point intersection(x1 + t * (x2 - x1), y1 + t * (y2 - y1));
    if (t < 0 || t > 1 || u < 0 || u > 1) {
        return Point(); // intersection point is not within the segments
    }
    return intersection; // segments intersect
}
// check if two segments intersect
bool do_segments_intersect(const Segment& seg1, const Segment& seg2) 
{
    double x1 = seg1.points[0].x, y1 = seg1.points[0].y;
    double x2 = seg1.points[1].x, y2 = seg1.points[1].y;
    double x3 = seg2.points[0].x, y3 = seg2.points[0].y;
    double x4 = seg2.points[1].x, y4 = seg2.points[1].y;

    // calculate the denominators
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (fabs(denom) < 1e-9) {
        return false; // Segments are parallel or coincident
    }

    // calculate the intersection point using parameter t and u
    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    double u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom;

    // check if t and u are within [0, 1]
    if (t < 0 || t > 1 || u < 0 || u > 1){
        return false; // intersection point is not within the segments
    }
    if ((t == 0 && u == 0) || (t == 1 && u == 1) || (t == 0 && u == 1) || (t == 1 && u == 0)){
        return false; // segments have a common extremity, not considered as an intersection
    }

    return true; // segments intersect
}
// check the relative orientation of three points
int orientation(const Point& p, const Point& q, const Point& r) 
{
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0){
        return 0;  // collinear
    }
    return (val > 0) ? 1 : 2;  // clockwise or counterclockwise
}
// check if two segments are collinear
bool are_segments_collinear(const Segment& seg1, const Segment& seg2) {
    return orientation(seg1.points[0], seg1.points[1], seg2.points[0]) == 0 && orientation(seg1.points[0], seg1.points[1], seg2.points[1]) == 0;
}
// check if seg2 is contained within seg1
bool is_segment_contained_within(const Segment& seg1, const Segment& seg2) {
    return is_point_on_segment(seg2.points[0], seg1) && is_point_on_segment(seg2.points[1], seg1);
}
// check if two segments are overlapping
bool are_overlapping(const Segment& seg1, const Segment& seg2)
{
    Point p1 = seg1.points[0];
    Point p2 = seg1.points[1];
    Point q1 = seg2.points[0];
    Point q2 = seg2.points[1];

    // check if orientation of (p1, p2, q1) and (p1, p2, q2) is collinear
    int o1 = orientation(p1, p2, q1);
    int o2 = orientation(p1, p2, q2);
    int o3 = orientation(q1, q2, p1);
    int o4 = orientation(q1, q2, p2);

    // if o1 == o2 == o3 == o4 == 0, we check if the segments overlap
    if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0){
        return is_point_on_segment(p1, seg2) || is_point_on_segment(p2, seg2) || is_point_on_segment(q1, seg1) || is_point_on_segment(q2, seg1);
    }
    
    // otherwise, they are not overlapping
    return false;
}
// check if a segment intersects with an arc circle
bool do_segment_intersect_arc_circle(const Segment& segment, const ArcCircle& arc_circle) 
{
    // parametric equations for the segment
    Point p1 = segment.points[0];
    Point p2 = segment.points[1];

    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;

    // Circle equation: (x - Cx)^2 + (y - Cy)^2 = R^2
    // Arc center and radius
    Point center = arc_circle.center;
    double radius = arc_circle.radius;
    
    double Cx = center.x, Cy = center.y;

    // Parametric equation for the segment: P(t) = (x1 + t * (x2 - x1), y1 + t * (y2 - y1)), t ∈ [0, 1]
    // Solve for t where the segment intersects the circle

    double dx = x2 - x1;
    double dy = y2 - y1;
    double fx = x1 - Cx;
    double fy = y1 - Cy;

    // Solve the quadratic equation to find the intersection points
    double a = dx * dx + dy * dy;
    double b = 2 * (fx * dx + fy * dy);
    double c = fx * fx + fy * fy - radius * radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0){
        return false; // no intersection
    }

    // find the intersection points (t values)
    discriminant = sqrt(discriminant);
    double t1 = (-b - discriminant) / (2 * a);
    double t2 = (-b + discriminant) / (2 * a);

    // Check if any intersection points are within the segment (t ∈ [0, 1])
    vector<Point> intersection_points;
    if (t1 >= 0 && t1 <= 1){
        intersection_points.push_back(Point(x1 + t1 * dx, y1 + t1 * dy));
    }
    if (t2 >= 0 && t2 <= 1){
        intersection_points.push_back(Point(x1 + t2 * dx, y1 + t2 * dy));
    }

    // check if the intersection points are within the arc's angular range
    for (const Point& p : intersection_points){
        double angle = atan2(p.y - Cy, p.x - Cx);
        if (angle < 0) angle += 2 * M_PI; // Normalize to [0, 2π]
        double start_angle = fmod(arc_circle.angle1, 2 * M_PI);
        double end_angle = fmod(arc_circle.angle2, 2 * M_PI);

        // check if the angle lies within [angle1, angle2]
        bool within_arc = false;
        if (start_angle <= end_angle){
            within_arc = (angle >= start_angle && angle <= end_angle);
        }
        else { // arc wraps around 0 radians
            within_arc = (angle >= start_angle || angle <= end_angle);
        }
        // valid intersection point
        if (within_arc){
            return true; 
        }
    }

    return false; // no valid intersection within the arc bounds
}
// check if a segment intersects with a parametric arc
bool do_segment_intersect_parametric(const Segment& segment, const ArcParameter& param_arc)
{
    // sample the parametric arc at multiple points
    const int num_samples = 100; // 100 should be enough for most cases
    for (int i = 0; i <= num_samples; ++i) {
        double t = i / (double)num_samples;
        Point p = param_arc.parametrisation(t);  // get the point on the parametric arc
        
        // check if the point on the parametric arc lies on the segment
        if (do_segments_intersect(segment, Segment(segment.points[0], p))) {
            return true;
        }
    }
    return false;
}
// check if two arc circles intersect
bool do_arc_circle_intersect(const ArcCircle& arc1, const ArcCircle& arc2)
{
    // calculate the distance between the centers of the two circles
    double dist = sqrt((arc1.center.x - arc2.center.x) * (arc1.center.x - arc2.center.x) + (arc1.center.y - arc2.center.y) * (arc1.center.y - arc2.center.y));

    //iIf the distance between the centers is greater than the sum of their radius, no intersection
    if (dist > arc1.radius + arc2.radius){
        return false;
    }

    // if the distance is less than the absolute difference of their radius, one circle is inside the other
    if (dist < abs(arc1.radius - arc2.radius)){
        return false;
    }

    // calculate the angle range for both arcs
    // normalize the angles within the range [0, 2*pi)
    double angle1_start = fmod(arc1.angle1, 2 * M_PI);
    double angle1_end = fmod(arc1.angle2, 2 * M_PI);
    double angle2_start = fmod(arc2.angle1, 2 * M_PI);
    double angle2_end = fmod(arc2.angle2, 2 * M_PI);

    // normalize the angular ranges if the arc is clockwise
    if (angle1_start > angle1_end){
        std::swap(angle1_start, angle1_end);  // swap the start and end for clockwise arc
    }

    if (angle2_start > angle2_end){
        std::swap(angle2_start, angle2_end);  // swap the start and end for clockwise arc
    }

    // check if the arcs' angular ranges overlap
    if (angle1_start < angle2_end && angle1_end > angle2_start){
        // the arcs overlap in angular space
        return true;
    }

    return false;
}
// check if an arc circle intersects with a parametric arc
bool do_arc_circle_intersect_parametric(const ArcCircle& arc_circle, const ArcParameter& param_arc) 
{
    const int num_samples = 100; // 100 should be enough for most cases
    const double margin = 0.5;  // margin for floating-point comparison
    
    // First, check if the entire parametric arc lies within the arc circle
    for (int i = 0; i <= num_samples; ++i) {
        double t = i / (double)num_samples;
        Point p = param_arc.parametrisation(t);  // Get the point on the parametric arc
        
        // Check if the point lies within the arc circle's boundary
        double distance = sqrt((p.x - arc_circle.center.x) * (p.x - arc_circle.center.x) + 
                               (p.y - arc_circle.center.y) * (p.y - arc_circle.center.y));
        
        // Points must be within the radius plus margin
        if (distance > arc_circle.radius + margin) {
            return false; // If any point is outside, there's no intersection
        }
    }

    // Check for intersection with the boundary of the arc circle
    // Sample points and check if they lie on the arc circle's boundary
    for (int i = 0; i <= num_samples; ++i) {
        double t = i / (double)num_samples;
        Point p = param_arc.parametrisation(t);  // Get the point on the parametric arc
        
        double distance = sqrt((p.x - arc_circle.center.x) * (p.x - arc_circle.center.x) + 
                               (p.y - arc_circle.center.y) * (p.y - arc_circle.center.y));
        
        // Check if the point is close to the boundary of the arc circle
        if (abs(distance - arc_circle.radius) <= margin) {
            // Check if the point lies within the angular range of the arc
            double angle = atan2(p.y - arc_circle.center.y, p.x - arc_circle.center.x);
            
            // Normalize the angle to be within [0, 2 * pi)
            angle = fmod(angle, 2 * M_PI);
            if (angle < 0) angle += 2 * M_PI;
            
            if (angle >= arc_circle.angle1 && angle <= arc_circle.angle2) {
                return true; // The point is within the arc's angular range and close to the boundary
            }
        }
    }

    return false; // No intersection found
}
// check if two parametric arcs intersect
bool do_parametric_intersect(const ArcParameter& arc1, const ArcParameter& arc2) 
{
    const int num_samples = 100; // 100 should be enough for most cases
    const double margin = 0.5;  // margin for floating-point comparison
    
    // sample points along arc1
    for (int i = 0; i <= num_samples; ++i){
        double t1 = i / (double)num_samples;
        Point p1 = arc1.parametrisation(t1);  // get the point on the first parametric arc
        
        // sample points along arc2
        for (int j = 0; j <= num_samples; ++j){
            double t2 = j / (double)num_samples;
            Point p2 = arc2.parametrisation(t2);  // get the point on the second parametric arc
            
            // check if the points on the two arcs are within the margin of each other
            double distance = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
            if (distance <= margin){
                return true; // points are close enough, meaning the arcs intersect
            }
        }
    }
    return false; // no intersection found
}
